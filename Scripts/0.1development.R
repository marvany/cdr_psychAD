library(data.table)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(viridis)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

outwrite <- function(x,filename){fwrite(x,file.path("results/processed", filename))}

# Load the results
results_path <- "results/drug_pheno_assoc_anal_psychAD.csv"
dt <- fread(results_path)
outdir = "results/processed"

# ============= 1. Data Cleaning and Preparation =============

# Remove NA values in Estimate column
dt <- dt[!is.na(Estimate)]

# Create p-value column from either z or t test results
dt[, p_value := fcase(
  !is.na(`Pr(>|z|)`), `Pr(>|z|)`,
  !is.na(`Pr(>|t|)`), `Pr(>|t|)`,
  default = NA
)]

# Create analysis type column
dt[, analysis_type := fcase(
  !is.na(`Pr(>|z|)`), 'logistic',
  !is.na(`Pr(>|t|)`), 'linear',
  default = NA
)]

# Calculate FDR and Bonferroni adjusted p-values
dt[, fdr_adjusted := p.adjust(p_value, method = 'fdr')]
dt[, bonferroni_adjusted := p.adjust(p_value, method = 'bonferroni')]

# Add significance flags
dt[, significant_05 := fdr_adjusted < 0.05]
dt[, significant_01 := fdr_adjusted < 0.01]
dt[, significant_bonferroni := bonferroni_adjusted < 0.05]

# ============= 2. Summary Statistics =============

# Count of significant associations at different thresholds
sig_counts <- data.frame(
  threshold = c("p < 0.05", "FDR < 0.05", "FDR < 0.01", "Bonferroni < 0.05"),
  count = c(
    sum(dt$p_value < 0.05, na.rm = TRUE),
    sum(dt$fdr_adjusted < 0.05, na.rm = TRUE),
    sum(dt$fdr_adjusted < 0.01, na.rm = TRUE),
    sum(dt$bonferroni_adjusted < 0.05, na.rm = TRUE)
  )
)
print(sig_counts)
outwrite(sig_counts, 'sig_counts.csv')

# Top 10 most significant associations
top_assoc <- dt[order(p_value)][1:10]
print(top_assoc[, .(compound, phenotype, Estimate, p_value, fdr_adjusted)])

# ============= 3. Volcano Plot =============

# Add log10 p-value for volcano plot
dt[, neg_log10_p := -log10(p_value)]

# Create factor for coloring points
dt[, significance := factor(
  case_when(
    bonferroni_adjusted < 0.05 ~ "Bonferroni p < 0.05",
    fdr_adjusted < 0.05 ~ "FDR p < 0.05", 
    p_value < 0.05 ~ "p < 0.05",
    TRUE ~ "Not significant"
  ),
  levels = c("Bonferroni p < 0.05", "FDR p < 0.05", "p < 0.05", "Not significant")
)]

# Add labels for top significant compounds
dt[, to_label := fdr_adjusted < 0.01]

# Volcano plot
p_volcano <- ggplot(dt, aes(x = Estimate, y = neg_log10_p, color = significance)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = -log10(0.05/nrow(dt)), linetype = "dotted", color = "darkgrey") +
  geom_vline(xintercept = 0, linetype = "solid", color = "darkgrey") +
  scale_color_manual(values = c("red", "orange", "blue", "grey")) +
  labs(
    title = "Volcano Plot of Drug-Phenotype Associations",
    x = "Effect Size (Estimate)",
    y = "-log10(p-value)",
    color = "Significance"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Add labels for top hits
p_volcano_labeled <- p_volcano +
  geom_text_repel(
    data = dt[to_label == TRUE],
    aes(label = compound),
    size = 3,
    max.overlaps = 15,
    box.padding = 0.5
  )

ggsave(file.path(outdir,"volcano_plot.pdf"), p_volcano_labeled, width = 10, height = 8)

# ============= 4. Manhattan-style Plot =============

# Get top phenotypes with significant associations
top_phenotypes <- dt[fdr_adjusted < 0.1, .N, by = phenotype][order(-N)][1:20]$phenotype

# Filter data for top phenotypes
dt_top_pheno <- dt[phenotype %in% top_phenotypes]

# Arrange compounds by phenotype for plotting
dt_top_pheno[, compound_index := .GRP, by = compound]

# Manhattan-style plot
p_manhattan <- ggplot(dt_top_pheno, aes(x = compound_index, y = neg_log10_p, color = phenotype)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = -log10(0.05/nrow(dt)), linetype = "dotted", color = "darkgrey") +
  labs(
    title = "Manhattan Plot of Drug Associations Across Phenotypes",
    x = "Compounds (sorted)",
    y = "-log10(p-value)",
    color = "Phenotype"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(file.path(outdir,"manhattan_plot.pdf"), p_manhattan, width = 12, height = 8)

# ============= 5. Heatmap of Top Associations =============

# Get top drugs with significant associations
top_drugs <- dt[fdr_adjusted < 0.05, .N, by = compound][order(-N)][1:20]$compound

# Filter for top drugs and phenotypes
dt_heatmap <- dt[compound %in% top_drugs & phenotype %in% top_phenotypes]

# Create matrix for heatmap
heatmap_data <- dcast(dt_heatmap, compound ~ phenotype, value.var = "Estimate", fill = 0)
heatmap_matrix <- as.matrix(heatmap_data[, -1])
rownames(heatmap_matrix) <- heatmap_data$compound

# Create significance matrix
significance_matrix <- dcast(dt_heatmap, compound ~ phenotype, value.var = "fdr_adjusted", fill = 1)
significance_matrix <- as.matrix(significance_matrix[, -1])
rownames(significance_matrix) <- heatmap_data$compound

# Set up color functions
col_fun <- colorRamp2(c(min(heatmap_matrix), 0, max(heatmap_matrix)), c("blue", "white", "red"))
sig_fun <- function(x) ifelse(x < 0.05, "*", "")

# Generate PDF heatmap using ComplexHeatmap
pdf(file.path(outdir, "association_heatmap.pdf"), width = 12, height = 10)
Heatmap(
  heatmap_matrix,
  name = "Effect Size",
  col = col_fun,
  rect_gp = gpar(col = "white", lwd = 1),
  column_title = "Drug-Phenotype Association Heatmap",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(significance_matrix[i, j] < 0.05) {
      grid.text("*", x, y, gp = gpar(fontsize = 14))
    }
    if(significance_matrix[i, j] < 0.01) {
      grid.text("**", x, y, gp = gpar(fontsize = 14))
    }
    if(significance_matrix[i, j] < 0.001) {
      grid.text("***", x, y, gp = gpar(fontsize = 14))
    }
  },
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  show_row_dend = TRUE,
  show_column_dend = TRUE
)
dev.off()

# ============= 6. Barplot of Top Drugs =============

# Count significant associations per drug
drug_sig_counts <- dt[fdr_adjusted < 0.05, .N, by = compound][order(-N)]
top_n_drugs <- drug_sig_counts[1:15]

# Create barplot
p_drug_counts <- ggplot(top_n_drugs, aes(x = reorder(compound, N), y = N)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top Drugs by Number of Significant Associations",
    x = NULL,
    y = "Number of Significant Associations (FDR < 0.05)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10)
  )

ggsave(file.path(outdir, "top_drugs_barplot.pdf"), p_drug_counts, width = 10, height = 8)

# ============= 7. Network Visualization =============

# Filter for significant associations
dt_network <- dt[fdr_adjusted < 0.05]

# Limit to top associations for readability
if(nrow(dt_network) > 100) {
  dt_network <- dt_network[order(fdr_adjusted)][1:100]
}

# Export network data for visualization in external tools like Cytoscape or Gephi
outwrite(dt_network[, .(compound, phenotype, Estimate, p_value, fdr_adjusted)], 
          "network_data.csv")

# ============= 8. Effect Size Distribution =============

# Histogram of effect sizes
p_effect_dist <- ggplot(dt, aes(x = Estimate, fill = significance)) +
  geom_histogram(bins = 50, position = "identity", alpha = 0.7) +
  scale_fill_manual(values = c("red", "orange", "blue", "grey")) +
  labs(
    title = "Distribution of Effect Sizes",
    x = "Effect Size (Estimate)",
    y = "Count",
    fill = "Significance"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(file.path(outdir, "effect_size_distribution.pdf"), p_effect_dist, width = 10, height = 6)

# ============= 9. Summary Table of Significant Results =============

# Create a summary table of significant results
sig_results <- dt[fdr_adjusted < 0.05][order(fdr_adjusted)]
sig_results <- sig_results[, .(
  compound, 
  phenotype, 
  analysis_type,
  Estimate, 
  `Std. Error`, 
  p_value, 
  fdr_adjusted, 
  bonferroni_adjusted
)]

# Export to CSV
outwrite(sig_results, "significant_associations.csv")

# Print summary
cat(sprintf("Total associations analyzed: %d\n", nrow(dt)))
cat(sprintf("Significant associations (FDR < 0.05): %d\n", nrow(sig_results)))
cat(sprintf("Number of unique compounds with significant associations: %d\n", 
            length(unique(sig_results$compound))))
cat(sprintf("Number of unique phenotypes with significant associations: %d\n", 
            length(unique(sig_results$phenotype))))

# ============= 10. Phenotype Analysis =============

# Count significant associations per phenotype
phenotype_sig_counts <- dt[fdr_adjusted < 0.05, .N, by = phenotype][order(-N)]
top_n_phenotypes <- phenotype_sig_counts[1:15]

# Create barplot for phenotypes
p_phenotype_counts <- ggplot(top_n_phenotypes, aes(x = reorder(phenotype, N), y = N)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  coord_flip() +
  labs(
    title = "Top Phenotypes by Number of Significant Drug Associations",
    x = NULL,
    y = "Number of Significant Associations (FDR < 0.05)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10)
  )

ggsave(file.path(outdir, "top_phenotypes_barplot.pdf"), p_phenotype_counts, width = 10, height = 8)

# Create an index HTML file to view all visualizations
cat('<!DOCTYPE html>
<html>
<head>
    <title>Drug-Phenotype Association Analysis Results</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1, h2 { color: #333366; }
        .container { display: flex; flex-wrap: wrap; }
        .viz { margin: 10px; padding: 10px; border: 1px solid #ddd; }
        img { max-width: 100%; height: auto; }
        table { border-collapse: collapse; width: 100%; }
        th, td { text-align: left; padding: 8px; border-bottom: 1px solid #ddd; }
        th { background-color: #f2f2f2; }
    </style>
</head>
<body>
    <h1>Drug-Phenotype Association Analysis Results</h1>
    
    <h2>Summary Statistics</h2>
    <p>Total associations analyzed: ', nrow(dt), '</p>
    <p>Significant associations (FDR < 0.05): ', nrow(sig_results), '</p>
    <p>Number of unique compounds with significant associations: ', length(unique(sig_results$compound)), '</p>
    <p>Number of unique phenotypes with significant associations: ', length(unique(sig_results$phenotype)), '</p>
    
    <div class="container">
        <div class="viz">
            <h2>Volcano Plot</h2>
            <img src="volcano_plot.png" alt="Volcano Plot">
            <p>This plot shows effect sizes vs. statistical significance, highlighting the most significant drug-phenotype associations.</p>
        </div>
        
        <div class="viz">
            <h2>Manhattan Plot</h2>
            <img src="manhattan_plot.png" alt="Manhattan Plot">
            <p>Drug associations across different phenotypes, highlighting significant hits.</p>
        </div>
        
        <div class="viz">
            <h2>Association Heatmap</h2>
            <img src="association_heatmap.png" alt="Heatmap">
            <p>Heatmap showing the effect sizes of top drug-phenotype associations. Asterisks indicate statistical significance.</p>
        </div>
        
        <div class="viz">
            <h2>Top Drugs</h2>
            <img src="top_drugs_barplot.png" alt="Top Drugs">
            <p>Drugs with the highest number of significant phenotype associations.</p>
        </div>
        
        <div class="viz">
            <h2>Top Phenotypes</h2>
            <img src="top_phenotypes_barplot.png" alt="Top Phenotypes">
            <p>Phenotypes with the highest number of significant drug associations.</p>
        </div>
        
        <div class="viz">
            <h2>Effect Size Distribution</h2>
            <img src="effect_size_distribution.png" alt="Effect Size Distribution">
            <p>Distribution of effect sizes across all analyzed drug-phenotype pairs.</p>
        </div>
    </div>
    
    <h2>Top Significant Associations</h2>
    <table>
        <tr>
            <th>Compound</th>
            <th>Phenotype</th>
            <th>Effect Size</th>
            <th>P-value</th>
            <th>FDR-adjusted P</th>
        </tr>', paste0(apply(top_assoc[1:10, .(compound, phenotype, Estimate, p_value, fdr_adjusted)], 1, function(row) {
        sprintf("<tr><td>%s</td><td>%s</td><td>%.4f</td><td>%.2e</td><td>%.2e</td></tr>", 
                row[1], row[2], as.numeric(row[3]), as.numeric(row[4]), as.numeric(row[5]))
    }), collapse=""), '
    </table>
    
    <p><a href="significant_associations.csv">Download complete list of significant associations</a></p>
    <p><a href="network_data.csv">Download network data for visualization in Cytoscape/Gephi</a></p>
    
    <footer>
        <p>Generated on ', format(Sys.time(), "%Y-%m-%d"), '</p>
    </footer>
</body>
</html>', file = "association_results.html")

# Convert PDFs to PNGs for HTML display
pdf_files <- c("volcano_plot.pdf", "manhattan_plot.pdf", "top_drugs_barplot.pdf", 
               "top_phenotypes_barplot.pdf", "effect_size_distribution.pdf")

# Check if ImageMagick is available
has_convert <- system("which convert", intern = TRUE)
if(length(has_convert) > 0) {
  for(pdf_file in pdf_files) {
    png_file <- sub("\\.pdf$", ".png", pdf_file)
    system(paste("convert -density 150", pdf_file, png_file))
  }
}

# For ComplexHeatmap PDF
if(length(has_convert) > 0) {
  system("convert -density 150 association_heatmap.pdf association_heatmap.png")
}

cat("Analysis complete. Results are available in the current directory.\n")
cat("Open 'association_results.html' to view a summary of all visualizations.\n")