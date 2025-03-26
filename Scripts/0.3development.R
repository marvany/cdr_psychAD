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

# Define AD-related keywords and compounds
ad_keywords <- c(
  'alzheimer', 'dementia', 'amyloid', 'tau', 'cholinesterase', 
  'acetylcholine', 'nmda', 'memantine', 'donepezil', 'rivastigmine', 
  'galantamine', 'aducanumab', 'cognit', 'memory'
)



dt[compound %in% known_ad_drugs & fdr.adjusted < 0.05, .(phenotype, compound)]
dt[compound %in% known_ad_drugs & bonferroni.adjusted < 0.05, .(phenotype, compound)]

# AD-related phenotypes
ad_phenotypes <- c('AD', 'Dementia', 'BRAAK_AD', 'CDRScore', 'Cognitive_Resilience')

# Load the results
results_path <- "results/drug_pheno_assoc_anal_psychAD.csv"
dt <- fread(results_path)

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

# Identify AD-related compounds by name
is_ad_compound <- function(compound) {
  compound_lower <- tolower(as.character(compound))
  return(any(sapply(ad_keywords, function(kw) grepl(kw, compound_lower))))
}

# Find AD compounds in our dataset
all_compounds <- unique(dt$compound)
ad_compounds <- all_compounds[sapply(all_compounds, is_ad_compound)]

# If we don't have enough compounds by keyword, add known AD drugs found in dataset
if(length(ad_compounds) < 10) {
  known_drugs_found <- all_compounds[sapply(all_compounds, function(compound) {
    compound_lower <- tolower(as.character(compound))
    return(any(sapply(known_ad_drugs, function(drug) grepl(drug, compound_lower))))
  })]
  
  ad_compounds <- unique(c(ad_compounds, known_drugs_found))
}

# If still not enough, add compounds that have significant associations with AD phenotypes
if(length(ad_compounds) < 15) {
  ad_phenotype_associations <- dt[phenotype %in% ad_phenotypes & significant == TRUE]
  ad_significant_compounds <- unique(ad_phenotype_associations$compound)
  
  ad_compounds <- unique(c(ad_compounds, ad_significant_compounds))
}

# Take top 30 compounds if we have more than that
if(length(ad_compounds) > 30) {
  # Count significant associations per compound
  ad_compound_counts <- dt[compound %in% ad_compounds & significant == TRUE, 
                         .N, by = compound][order(-N)]
  ad_compounds <- ad_compound_counts$compound[1:30]
}

# Extract data for AD compounds
ad_data <- dt[compound %in% ad_compounds]

# ============= 1. Number of significant associations per AD compound =============

# Count significant associations per compound
ad_compound_counts <- ad_data[significant == TRUE, .N, by = compound][order(-N)]
ad_compound_counts[, compound := fct_reorder(compound, N)]

# Create bar plot
p1 <- ggplot(ad_compound_counts, aes(x = compound, y = N)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = N), hjust = -0.2, size = 3) +
  coord_flip() +
  labs(
    title = "Number of Significant Associations per AD Compound",
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
ggsave("ad_compound_significant_associations.pdf", p1, width = 10, height = max(8, nrow(ad_compound_counts) * 0.25))

# ============= 2. Top phenotypes with significant associations to AD compounds =============

# Count significant associations per phenotype for AD compounds
ad_phenotype_counts <- ad_data[significant == TRUE, .N, by = phenotype][order(-N)]
ad_phenotype_counts[, phenotype := fct_reorder(phenotype, N)]

# Filter to phenotypes with at least one significant association
ad_phenotype_counts <- ad_phenotype_counts[N > 0]

# Create bar plot
p2 <- ggplot(ad_phenotype_counts, aes(x = phenotype, y = N)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  geom_text(aes(label = N), hjust = -0.2, size = 3) +
  coord_flip() +
  labs(
    title = "Phenotypes with Significant Associations to AD Compounds",
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
ggsave("ad_phenotype_significant_associations.pdf", p2, width = 10, height = max(8, nrow(ad_phenotype_counts) * 0.25))

# ============= 3. Heatmap of top AD compounds by top phenotypes =============

# Get top compounds and phenotypes
top_n_compounds <- min(15, nrow(ad_compound_counts))
top_n_phenotypes <- min(15, nrow(ad_phenotype_counts))

top_compounds <- ad_compound_counts$compound[1:top_n_compounds]
top_phenotypes <- ad_phenotype_counts$phenotype[1:top_n_phenotypes]

# Create a matrix of effect sizes
heatmap_data <- ad_data[compound %in% top_compounds & phenotype %in% top_phenotypes]
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
pdf("ad_compound_phenotype_heatmap.pdf", width = 12, height = 10)
Heatmap(
  heatmap_matrix_values,
  name = "Effect Size",
  col = col_fun,
  rect_gp = gpar(col = "white", lwd = 1),
  column_title = "AD Compounds - Phenotype Association Heatmap",
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

# ============= 4. Volcano plot of AD compound associations =============

# Add log10 p-value for volcano plot
ad_data[, neg_log10_p := -log10(p_value)]

# Create factor for coloring points
ad_data[, significance := factor(
  case_when(
    fdr_adjusted < 0.01 ~ "FDR p < 0.01",
    fdr_adjusted < 0.05 ~ "FDR p < 0.05", 
    p_value < 0.05 ~ "p < 0.05",
    TRUE ~ "Not significant"
  ),
  levels = c("FDR p < 0.01", "FDR p < 0.05", "p < 0.05", "Not significant")
)]

# Add labels for top hits
ad_data[, to_label := fdr_adjusted < 0.01]

# Create volcano plot
p3 <- ggplot(ad_data, aes(x = Estimate, y = neg_log10_p, color = significance)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = -log10(0.05/nrow(ad_data)), linetype = "dotted", color = "darkgrey") +
  geom_vline(xintercept = 0, linetype = "solid", color = "darkgrey") +
  scale_color_manual(values = c("red", "orange", "blue", "grey")) +
  labs(
    title = "Volcano Plot of AD Compound Associations",
    x = "Effect Size (Estimate)",
    y = "-log10(p-value)",
    color = "Significance"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Save plot
ggsave("ad_compound_volcano_plot.pdf", p3, width = 10, height = 8)

# ============= 5. Effect size distribution by compound =============

# Box plot of effect sizes for top compounds
top_compounds_data <- ad_data[compound %in% top_compounds]
top_compounds_data[, compound := factor(compound, levels = levels(ad_compound_counts$compound))]

p4 <- ggplot(top_compounds_data, aes(x = compound, y = Estimate, fill = compound)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Distribution of Effect Sizes for Top AD Compounds",
    x = NULL,
    y = "Effect Size (Estimate)"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Save plot
ggsave("ad_compound_effect_sizes.pdf", p4, width = 10, height = 8)

# ============= 6. Top significant associations table =============

# Extract top 20 most significant associations
top_associations <- ad_data[order(p_value)][1:20, .(
  compound, 
  phenotype, 
  Estimate, 
  `Std. Error`, 
  p_value, 
  fdr_adjusted
)]

# Export to CSV
write.csv(top_associations, "top_ad_compound_associations.csv", row.names = FALSE)

# ============= 7. Comparison of AD vs non-AD compounds =============

# Calculate percentage of significant associations
ad_sig_percentage <- nrow(ad_data[significant == TRUE]) / nrow(ad_data) * 100
non_ad_sig_percentage <- nrow(dt[!compound %in% ad_compounds & significant == TRUE]) / 
                         nrow(dt[!compound %in% ad_compounds]) * 100

comparison_data <- data.frame(
  compound_type = c("AD Compounds", "Non-AD Compounds"),
  significant_percentage = c(ad_sig_percentage, non_ad_sig_percentage)
)

p5 <- ggplot(comparison_data, aes(x = compound_type, y = significant_percentage, fill = compound_type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2f%%", significant_percentage)), vjust = -0.5, size = 4) +
  labs(
    title = "Percentage of Significant Associations",
    subtitle = "AD Compounds vs. Non-AD Compounds",
    x = NULL,
    y = "Percentage of Significant Associations (%)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Save plot
ggsave("ad_vs_nonad_comparison.pdf", p5, width = 8, height = 6)

# ============= 8. Analysis of AD phenotypes specifically =============

# Focus on the associations between AD compounds and AD phenotypes
ad_phenotype_specific <- ad_data[phenotype %in% ad_phenotypes]

# Count significant associations for each AD phenotype
ad_phenotype_specific_counts <- ad_phenotype_specific[significant == TRUE, .N, by = phenotype][order(-N)]

# Create bar plot
p6 <- ggplot(ad_phenotype_specific_counts, aes(x = reorder(phenotype, -N), y = N, fill = phenotype)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = N), vjust = -0.5, size = 4) +
  labs(
    title = "Significant Associations Between AD Compounds and AD Phenotypes",
    x = NULL,
    y = "Number of Significant Associations"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save plot
ggsave("ad_compound_ad_phenotype_associations.pdf", p6, width = 10, height = 7)

# ============= 9. Combined visualizations =============

# Combine plots using patchwork
combined_plot1 <- p1 / p2
combined_plot2 <- p3 | p5
combined_plot3 <- p4 / p6

ggsave("ad_compounds_analysis_1.pdf", combined_plot1, width = 12, height = 14)
ggsave("ad_compounds_analysis_2.pdf", combined_plot2, width = 14, height = 8)
ggsave("ad_compounds_analysis_3.pdf", combined_plot3, width = 12, height = 14)

# ============= 10. Create HTML report =============

html_output <- '<!DOCTYPE html>
<html>
<head>
<title>Alzheimer\'s Disease Compounds Association Analysis</title>
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
<h1>Alzheimer\'s Disease Compounds Association Analysis</h1>

<div class="summary">
<h2>Summary Statistics</h2>
<table>
<tr><th>Metric</th><th>Value</th></tr>
<tr><td>Number of AD compounds analyzed</td><td>', length(ad_compounds), '</td></tr>
<tr><td>Total associations for AD compounds</td><td>', nrow(ad_data), '</td></tr>
<tr><td>Significant associations (FDR < 0.05)</td><td>', nrow(ad_data[significant == TRUE]), '</td></tr>
<tr><td>Percentage of significant associations</td><td>', sprintf("%.2f%%", ad_sig_percentage), '</td></tr>
<tr><td>Percentage of significant associations (non-AD compounds)</td><td>', sprintf("%.2f%%", non_ad_sig_percentage), '</td></tr>
<tr><td>Top AD compound with most significant associations</td><td>', as.character(ad_compound_counts$compound[1]), ' (', ad_compound_counts$N[1], ')</td></tr>
<tr><td>Top phenotype for AD compounds</td><td>', as.character(ad_phenotype_counts$phenotype[1]), ' (', ad_phenotype_counts$N[1], ')</td></tr>
</table>
</div>

<div class="viz">
<h2>Number of Significant Associations per AD Compound</h2>
<p>This bar chart shows the count of phenotypes that have a significant association (FDR-adjusted p-value < 0.05) with each AD-related compound.</p>
<img src="ad_compound_significant_associations.png" alt="Bar chart of significant associations per AD compound">
</div>

<div class="viz">
<h2>Phenotypes with Significant Associations to AD Compounds</h2>
<p>This bar chart shows the phenotypes with the most significant associations to AD compounds.</p>
<img src="ad_phenotype_significant_associations.png" alt="Bar chart of phenotypes with significant associations to AD compounds">
</div>

<div class="viz">
<h2>AD Compound-Phenotype Association Heatmap</h2>
<p>This heatmap shows the effect sizes of associations between top AD compounds and phenotypes. Asterisks indicate statistical significance (* p<0.05, ** p<0.01, *** p<0.001).</p>
<img src="ad_compound_phenotype_heatmap.png" alt="Heatmap of AD compound-phenotype associations">
</div>

<div class="viz">
<h2>Volcano Plot of AD Compound Associations</h2>
<p>This plot shows effect sizes vs. statistical significance for all AD compound-