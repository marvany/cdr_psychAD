library(data.table)
library(ggplot2)
library(dplyr)
library(forcats)
library(RColorBrewer)
library(viridis)
library(scales)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(patchwork)

bcmarks = read_excel('/sc/arion/projects/va-biobank/PROJECTS/cdr.comparative.efficacy.marios/Resources/benchmark_drugs.xlsx')
bc = as.data.table(bcmarks)
ad.bc = bc[comb_indication == 'AD']
known_ad_drugs = ad.bc$pert_iname

# Load the results
results_path <- "results/drug_pheno_assoc_anal_psychAD.csv"
dt <- fread(results_path)
outdir3 = "results/processed3"

outwrite3 <- function(x,filename){fwrite(x,file.path("results/processed3", filename))}

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

# Define the ADBC compounds (this vector should already exist in your environment)
adbc <- known_ad_drugs


# Filter data for ADBC compounds
adbc_data <- dt[compound %in% adbc]

# ============= 1. Number of significant associations per ADBC compound =============
adbc_data[significant == TRUE,]
# Count significant associations per compound
adbc_compound_counts <- adbc_data[significant == TRUE, .N, by = compound][order(-N)]
adbc_compound_counts[, compound := fct_reorder(compound, N)]

# Create bar plot
p1 <- ggplot(adbc_compound_counts, aes(x = compound, y = N)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = N), hjust = -0.2, size = 3) +
  coord_flip() +
  labs(
    title = "Number of Significant Associations per ADBC Compound",
    subtitle = "FDR-adjusted p-value < 0.05",
    x = NULL,
    y = "Number of Significant Associations"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 10)
  )

# Save plot
ggsave("adbc_compound_significant_associations.pdf", p1, width = 10, 
       height = max(8, nrow(adbc_compound_counts) * 0.25))

# ============= 2. Top phenotypes with significant associations to ADBC compounds =============

# Count significant associations per phenotype for ADBC compounds
adbc_phenotype_counts <- adbc_data[significant == TRUE, .N, by = phenotype][order(-N)]
adbc_phenotype_counts[, phenotype := fct_reorder(phenotype, N)]

# Filter to phenotypes with at least one significant association
adbc_phenotype_counts <- adbc_phenotype_counts[N > 0]

# Create bar plot
p2 <- ggplot(adbc_phenotype_counts, aes(x = phenotype, y = N)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  geom_text(aes(label = N), hjust = -0.2, size = 3) +
  coord_flip() +
  labs(
    title = "Phenotypes with Significant Associations to ADBC Compounds",
    subtitle = "FDR-adjusted p-value < 0.05",
    x = NULL,
    y = "Number of Significant Associations"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 10)
  )

# Save plot
ggsave("adbc_phenotype_significant_associations.pdf", p2, width = 10, 
       height = max(8, nrow(adbc_phenotype_counts) * 0.25))

# ============= 3. Heatmap of ADBC compounds by top phenotypes =============

# Get top phenotypes
top_n_phenotypes <- min(20, nrow(adbc_phenotype_counts))
top_phenotypes <- adbc_phenotype_counts$phenotype[1:top_n_phenotypes]

# Create a matrix of effect sizes
heatmap_data <- adbc_data[phenotype %in% top_phenotypes]
heatmap_matrix <- dcast(heatmap_data, compound ~ phenotype, value.var = "Estimate", fill = 0)
heatmap_matrix_values <- as.matrix(heatmap_matrix[, -1])
rownames(heatmap_matrix_values) <- heatmap_matrix$compound

# Create a matrix of p-values for significance markers
significance_matrix <- dcast(heatmap_data, compound ~ phenotype, value.var = "fdr_adjusted", fill = 1)
significance_matrix_values <- as.matrix(significance_matrix[, -1])
rownames(significance_matrix_values) <- significance_matrix$compound

# Set up color function
col_fun <- colorRamp2(
  c(min(heatmap_matrix_values), 0, max(heatmap_matrix_values)), 
  c("blue", "white", "red")
)

# Generate PDF heatmap
pdf("adbc_compound_phenotype_heatmap.pdf", width = 12, height = 10)
Heatmap(
  heatmap_matrix_values,
  name = "Effect Size",
  col = col_fun,
  rect_gp = gpar(col = "white", lwd = 1),
  column_title = "ADBC Compounds - Phenotype Association Heatmap",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(significance_matrix_values[i, j] < 0.05) {
      grid.text("*", x, y, gp = gpar(fontsize = 14))
    }
    if(significance_matrix_values[i, j] < 0.01) {
      grid.text("**", x, y, gp = gpar(fontsize = 14))
    }
    if(significance_matrix_values[i, j] < 0.001) {
      grid.text("***", x, y, gp = gpar(fontsize = 14))
    }
  },
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  show_row_dend = TRUE,
  show_column_dend = TRUE
)
dev.off()

# ============= 4. Volcano plot of ADBC compound associations =============

# Add log10 p-value for volcano plot
adbc_data[, neg_log10_p := -log10(p_value)]

# Create factor for coloring points
adbc_data[, significance := factor(
  case_when(
    fdr_adjusted < 0.01 ~ "FDR p < 0.01",
    fdr_adjusted < 0.05 ~ "FDR p < 0.05", 
    p_value < 0.05 ~ "p < 0.05",
    TRUE ~ "Not significant"
  ),
  levels = c("FDR p < 0.01", "FDR p < 0.05", "p < 0.05", "Not significant")
)]

# Add compound information
adbc_data[, label := paste0(compound, " - ", phenotype)]

# Add labels for top hits
adbc_data[, to_label := fdr_adjusted < 0.01]

# Create volcano plot
p3 <- ggplot(adbc_data, aes(x = Estimate, y = neg_log10_p, color = significance)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = -log10(0.05/nrow(adbc_data)), linetype = "dotted", color = "darkgrey") +
  geom_vline(xintercept = 0, linetype = "solid", color = "darkgrey") +
  scale_color_manual(values = c("red", "orange", "blue", "grey")) +
  labs(
    title = "Volcano Plot of ADBC Compound Associations",
    x = "Effect Size (Estimate)",
    y = "-log10(p-value)",
    color = "Significance"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Add labels for top hits using ggrepel
if(requireNamespace("ggrepel", quietly = TRUE)) {
  library(ggrepel)
  p3 <- p3 + 
    geom_text_repel(
      data = adbc_data[to_label == TRUE],
      aes(label = label),
      size = 3,
      max.overlaps = 15,
      box.padding = 0.5
    )
}

# Save plot
ggsave("adbc_compound_volcano_plot.pdf", p3, width = 12, height = 10)

# ============= 5. Effect size by compound and phenotype category =============

# Define phenotype categories if possible
# This is a simple example - you might want to customize based on your specific phenotypes
adbc_data[, phenotype_category := case_when(
  grepl("cognit|memory|CDR|BRAAK|Dementia|AD", phenotype, ignore.case = TRUE) ~ "Cognitive",
  grepl("mood|depress|anxiety|MDD|BD|psych", phenotype, ignore.case = TRUE) ~ "Psychiatric",
  grepl("nps_|insomnia|sleep", phenotype, ignore.case = TRUE) ~ "Neuropsychiatric",
  grepl("motor|PD|move", phenotype, ignore.case = TRUE) ~ "Motor",
  TRUE ~ "Other"
)]

# Create a faceted box plot by compound and phenotype category
p4 <- ggplot(adbc_data, aes(x = compound, y = Estimate, fill = phenotype_category)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
  facet_wrap(~ phenotype_category, scales = "free_y") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Effect Sizes by ADBC Compound and Phenotype Category",
    x = NULL,
    y = "Effect Size (Estimate)",
    fill = "Phenotype Category"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

# Save plot
ggsave("adbc_compound_effect_by_category.pdf", p4, width = 12, height = 10)

# ============= 6. Top significant associations table =============

# Extract top 20 most significant associations
top_associations <- adbc_data[order(p_value)][1:min(20, nrow(adbc_data)), .(
  compound, 
  phenotype, 
  Estimate, 
  `Std. Error`, 
  p_value, 
  fdr_adjusted
)]

# Export to CSV
write.csv(top_associations, "top_adbc_compound_associations.csv", row.names = FALSE)

# ============= 7. Dot plot of significant associations =============

# Create a dot plot showing significance and effect sizes
significant_assocs <- adbc_data[significant == TRUE]

# Reorder compounds and phenotypes by number of significant associations
compound_order <- adbc_compound_counts$compound
phenotype_order <- adbc_phenotype_counts$phenotype

# Convert to factors with desired order
significant_assocs[, compound := factor(compound, levels = compound_order)]
significant_assocs[, phenotype := factor(phenotype, levels = phenotype_order)]

# Create dot plot
p5 <- ggplot(significant_assocs, 
            aes(x = phenotype, y = compound, size = abs(Estimate), color = Estimate)) +
  geom_point() +
  scale_size_continuous(name = "Effect Size (abs)", range = c(2, 10)) +
  scale_color_gradient2(
    name = "Effect Size",
    low = "blue", 
    mid = "white", 
    high = "red", 
    midpoint = 0
  ) +
  labs(
    title = "Significant Associations Between ADBC Compounds and Phenotypes",
    subtitle = "FDR-adjusted p-value < 0.05",
    x = "Phenotype",
    y = "ADBC Compound"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )

# Save plot
ggsave("adbc_significant_associations_dotplot.pdf", p5, width = 14, height = 10)

# ============= 8. Interactive HTML report =============

# Create an HTML report summarizing all the visualizations
html_output <- '<!DOCTYPE html>
<html>
<head>
<title>ADBC Compounds Phenotype Association Analysis</title>
<style>
body {font-family: Arial, sans-serif; margin: 20px; line-height: 1.6;}
h1, h2 {color: #2c3e50;}
.container {max-width: 1200px; margin: 0 auto;}
.summary {background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin-bottom: 20px;}
.viz {margin-bottom: 30px;}
table {border-collapse: collapse; width: 100%; margin: 20px 0;}
th, td {text-align: left; padding: 12px; border-bottom: 1px solid #ddd;}
th {background-color: #f2f2f2;}
img {max-width: 100%; height: auto; box-shadow: 0 4px 8px rgba(0,0,0,0.1);}
</style>
</head>
<body>
<div class="container">
<h1>ADBC Compounds Phenotype Association Analysis</h1>

<div class="summary">
<h2>Summary Statistics</h2>
<table>
<tr><th>Metric</th><th>Value</th></tr>
<tr><td>Number of ADBC compounds analyzed</td><td>'

html_output <- paste0(html_output, length(adbc), '</td></tr>
<tr><td>Total associations tested for ADBC compounds</td><td>', nrow(adbc_data), '</td></tr>
<tr><td>Significant associations (FDR < 0.05)</td><td>', nrow(adbc_data[significant == TRUE]), '</td></tr>
<tr><td>Percentage of significant associations</td><td>', 
                     sprintf("%.2f%%", nrow(adbc_data[significant == TRUE])/nrow(adbc_data)*100), '</td></tr>')

if(nrow(adbc_compound_counts) > 0) {
  html_output <- paste0(html_output, 
                       '<tr><td>Top ADBC compound with most significant associations</td><td>', 
                       as.character(adbc_compound_counts$compound[1]), ' (', 
                       adbc_compound_counts$N[1], ')</td></tr>')
}

if(nrow(adbc_phenotype_counts) > 0) {
  html_output <- paste0(html_output, 
                       '<tr><td>Top phenotype for ADBC compounds</td><td>', 
                       as.character(adbc_phenotype_counts$phenotype[1]), ' (', 
                       adbc_phenotype_counts$N[1], ')</td></tr>')
}

html_output <- paste0(html_output, '</table>
</div>

<div class="viz">
<h2>Number of Significant Associations per ADBC Compound</h2>
<p>This bar chart shows the count of phenotypes that have a significant association (FDR-adjusted p-value < 0.05) with each ADBC compound.</p>
<img src="adbc_compound_significant_associations.png" alt="Bar chart of significant associations per ADBC compound">
</div>

<div class="viz">
<h2>Phenotypes with Significant Associations to ADBC Compounds</h2>
<p>This bar chart shows the phenotypes with the most significant associations to ADBC compounds.</p>
<img src="adbc_phenotype_significant_associations.png" alt="Bar chart of phenotypes with significant associations to ADBC compounds">
</div>

<div class="viz">
<h2>ADBC Compound-Phenotype Association Heatmap</h2>
<p>This heatmap shows the effect sizes of associations between ADBC compounds and phenotypes. Asterisks indicate statistical significance (* p<0.05, ** p<0.01, *** p<0.001).</p>
<img src="adbc_compound_phenotype_heatmap.png" alt="Heatmap of ADBC compound-phenotype associations">
</div>

<div class="viz">
<h2>Volcano Plot of ADBC Compound Associations</h2>
<p>This plot shows effect sizes vs. statistical significance for all ADBC compound-phenotype associations.</p>
<img src="adbc_compound_volcano_plot.png" alt="Volcano plot of ADBC compound associations">
</div>

<div class="viz">
<h2>Effect Sizes by ADBC Compound and Phenotype Category</h2>
<p>These box plots show the distribution of effect sizes for each ADBC compound, categorized by phenotype type.</p>
<img src="adbc_compound_effect_by_category.png" alt="Box plots of effect sizes by compound and phenotype category">
</div>

<div class="viz">
<h2>Significant Associations Dot Plot</h2>
<p>This dot plot shows significant associations between ADBC compounds and phenotypes. The size represents the absolute effect size, and the color indicates the direction and magnitude of the effect.</p>
<img src="adbc_significant_associations_dotplot.png" alt="Dot plot of significant associations">
</div>

<h2>Top Significant Associations</h2>
<table>
<tr><th>Compound</th><th>Phenotype</th><th>Effect Size</th><th>P-value</th><th>FDR-adjusted P</th></tr>')

if(nrow(top_associations) > 0) {
  for(i in 1:nrow(top_associations)) {
    html_output <- paste0(html_output, 
                         '<tr>',
                         '<td>', top_associations$compound[i], '</td>',
                         '<td>', top_associations$phenotype[i], '</td>',
                         '<td>', sprintf("%.4f", top_associations$Estimate[i]), '</td>',
                         '<td>', format(top_associations$p_value[i], scientific = TRUE, digits = 3), '</td>',
                         '<td>', format(top_associations$fdr_adjusted[i], scientific = TRUE, digits = 3), '</td>',
                         '</tr>')
  }
}

html_output <- paste0(html_output, '</table>
<p><a href="top_adbc_compound_associations.csv">Download complete table (CSV)</a></p>

<footer>
<p>Generated on ', format(Sys.time(), "%Y-%m-%d"), '</p>
</footer>
</div>
</body>
</html>')

# Write HTML file
writeLines(html_output, "adbc_phenotype_associations.html")

# Convert PDFs to PNGs for HTML display
pdf_files <- c(
  "adbc_compound_significant_associations.pdf",
  "adbc_phenotype_significant_associations.pdf",
  "adbc_compound_phenotype_heatmap.pdf",
  "adbc_compound_volcano_plot.pdf",
  "adbc_compound_effect_by_category.pdf",
  "adbc_significant_associations_dotplot.pdf"
)

# Check if ImageMagick is available
has_convert <- system("which convert", intern = TRUE)
if(length(has_convert) > 0) {
  for(pdf_file in pdf_files) {
    png_file <- sub("\\.pdf$", ".png", pdf_file)
    system(paste("convert -density 150", pdf_file, png_file))
  }
}

cat("Analysis complete. Visualizations of ADBC compound associations across phenotypes are available.\n")
cat("Open 'adbc_phenotype_associations.html' to view a summary of all visualizations.\n")