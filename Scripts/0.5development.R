
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# Define phenotypes of interest

phenotypes_of_interest = 'prs_scaled_IQ'
outdir = 'results/processedHeatmapscaled_IQ'

phenotypes_of_interest <- "prs_scaled_mem"
outdir = 'results/processedHeatmap_scaled_mem'


phenotypes_of_interest <- c(
  "Dementia",
  "prs_scaled_AD_Bellenguez",
  "prs_scaled_AD2",
  "AmygPlaquesValue",
  "AmygTanglesValue", 
  "CerebralAtrophy",
  "MidGliosisValue",
  "SupGliosisValue",
  "CDR_Memory",
  "AD",
  "Cognitive_and_Tau_Resilience"
)
outdir = 'results/processed5Heatmap'

# Load the results
results_path <- "results/drug_pheno_assoc_anal_psychAD.csv"
dt <- fread(results_path)

# Remove NA values in Estimate column
dt <- dt[!is.na(Estimate)]

# Create p-value column from either z or t test results
dt[, p_value := fcase(
  !is.na(`Pr(>|z|)`), `Pr(>|z|)`,
  !is.na(`Pr(>|t|)`), `Pr(>|t|)`,
  default = NA
)]

# Calculate FDR adjusted p-values
dt[, fdr_adjusted := p.adjust(p_value, method = 'fdr')]

# Add significance flags
dt[, significant := fdr_adjusted < 0.05]

# Check which phenotypes of interest exist in the data
existing_phenotypes <- phenotypes_of_interest[phenotypes_of_interest %in% unique(dt$phenotype)]

if(length(existing_phenotypes) == 0) {
  # If specified phenotypes don't exist, use top phenotypes with most significant compounds
  phenotype_counts <- dt[significant == TRUE, .N, by = phenotype][order(-N)][1:13]
  existing_phenotypes <- phenotype_counts$phenotype
  cat("Specified phenotypes not found. Using top phenotypes by significant compounds instead:\n")
  print(existing_phenotypes)
}

# Find compounds that are significant for at least one of these phenotypes
significant_compounds <- dt[phenotype %in% existing_phenotypes & significant == TRUE, unique(compound)]


# Check number of compounds per phenotypes and their total sum
dt[compound %in% significant_compounds & phenotype %in% existing_phenotypes & significant == TRUE][, .N, by = .(phenotype)]
sum(dt[compound %in% significant_compounds & phenotype %in% existing_phenotypes & significant == TRUE][, .N, by = .(phenotype)]$N)

if(length(significant_compounds) > 0){

  # Filter data for these compounds and phenotypes
  filtered_data <- dt[compound %in% significant_compounds & phenotype %in% existing_phenotypes]
  
  # Create a matrix for the heatmap
  matrix_data <- dcast(filtered_data, compound ~ phenotype, value.var = "Estimate", fill = 0)
  
  # Convert to matrix format
  heatmap_matrix <- as.matrix(matrix_data[, -1])
  rownames(heatmap_matrix) <- matrix_data$compound
  
  # Create a significance matrix
  sig_matrix <- dcast(filtered_data, compound ~ phenotype, value.var = "fdr_adjusted", fill = 1)
  sig_matrix_values <- as.matrix(sig_matrix[, -1])
  rownames(sig_matrix_values) <- sig_matrix$compound
  
  # If the matrix is too large, limit to top compounds
  if(nrow(heatmap_matrix) > 100) {
    # Count significant associations per compound
    compound_sig_counts <- rowSums(sig_matrix_values < 0.05)
    
    # Get the top 50 compounds by number of significant associations
    top_compounds <- names(sort(compound_sig_counts, decreasing = TRUE)[1:50])
    
    # Filter matrices
    heatmap_matrix <- heatmap_matrix[top_compounds, ]
    sig_matrix_values <- sig_matrix_values[top_compounds, ]
  }
  
  # Set up color function with better contrast
  max_abs_value <- max(abs(heatmap_matrix))
  col_fun <- colorRamp2(
    c(-max_abs_value, 0, max_abs_value), 
    c("#0072B2", "white", "#D55E00")  # Colorblind-friendly blue and orange
  )
  
  # Calculate the height based on number of compounds
  height <- max(12, nrow(heatmap_matrix) * 0.25)

  if(!dir.exists(outdir)) dir.create(outdir)  
  # Generate the combined heatmap
  pdf(file.path(outdir,"combined_phenotype_compound_heatmap.pdf"), width = 14, height = height)


  show.col.dend = ifelse(length(phenotypes_of_interest) > 1, TRUE, FALSE)
  
  if(length(phenotypes_of_interest) > 1) {
    clust.dist.columns <- "euclidean"
  } else {
    clust.dist.columns <- NULL
  }
  clust.cols  = ifelse(length(phenotypes_of_interest) > 1, TRUE, FALSE)
  if(length(phenotypes_of_interest) == 1){
    # Convert to a proper matrix with one column
    heatmap_matrix <- matrix(heatmap_matrix, ncol=1)
    rownames(heatmap_matrix) <- names(heatmap_matrix)
    colnames(heatmap_matrix) <- phenotypes_of_interest[1]

    # Do the same for significance values
    sig_matrix_values <- matrix(sig_matrix_values, ncol=1)
    rownames(sig_matrix_values) <- names(sig_matrix_values)
  }

  ht <- Heatmap(
    heatmap_matrix,
    name = "Effect Size",
    col = col_fun,
    cluster_columns = clust.cols,
    rect_gp = gpar(col = "white", lwd = 0.5),
    column_title = "Effect Sizes Across AD-Related Phenotypes",
    column_title_gp = gpar(fontsize = 16, fontface = "bold"),
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    row_title = "Compounds",
    row_title_gp = gpar(fontsize = 14, fontface = "bold"),
    column_title_side = "top",
    row_title_side = "left",
    cell_fun = function(j, i, x, y, width, height, fill) {
      # Use ifelse to handle both single and multiple phenotype cases
      sig_value <- ifelse(length(phenotypes_of_interest) == 1, 
                          sig_matrix_values[i, 1], 
                          sig_matrix_values[i, j])
      
      if(sig_value < 0.05) {
        grid.text("*", x, y, gp = gpar(fontsize = 12))
      }
      if(sig_value < 0.01) {
        grid.text("**", x, y, gp = gpar(fontsize = 12))
      }
      if(sig_value < 0.001) {
        grid.text("***", x, y, gp = gpar(fontsize = 12))
      }
    },
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = clust.dist.columns,
    show_row_dend = TRUE,
    show_column_dend = show.col.dend,
    heatmap_legend_param = list(
      title = "Effect Size",
      title_gp = gpar(fontsize = 12, fontface = "bold"),
      labels_gp = gpar(fontsize = 10)
    )
  )
  
  # Draw the heatmap
  draw(ht, heatmap_legend_side = "right")
  
  # Add a legend for significance
  pushViewport(viewport(x = 0.9, y = 0.05, width = 0.2, height = 0.1, just = c("right", "bottom")))
  grid.text("Significance:", x = 0, y = 0.8, just = "left", gp = gpar(fontsize = 10, fontface = "bold"))
  grid.text("* FDR < 0.05", x = 0, y = 0.6, just = "left", gp = gpar(fontsize = 10))
  grid.text("** FDR < 0.01", x = 0, y = 0.4, just = "left", gp = gpar(fontsize = 10))
  grid.text("*** FDR < 0.001", x = 0, y = 0.2, just = "left", gp = gpar(fontsize = 10))
  upViewport()
  
  dev.off()
  
  # Convert PDF to PNG for easier viewing
  has_convert <- system("which convert", intern = TRUE)
  if(length(has_convert) > 0) {
    system("convert -density 150 combined_phenotype_compound_heatmap.pdf combined_phenotype_compound_heatmap.png")
  }
  
  cat("Combined heatmap created successfully with", nrow(heatmap_matrix), "compounds and", ncol(heatmap_matrix), "phenotypes.\n")
  
  # Check if ADBC compounds are included
#  if(exists("adbc")) {
  if(F){
    adbc_in_heatmap <- sum(rownames(heatmap_matrix) %in% adbc)
    cat(adbc_in_heatmap, "out of", length(adbc), "ADBC compounds are included in the heatmap.\n")
  }
} else {
  cat("No significant compounds found for any of the specified phenotypes.\n")
}
