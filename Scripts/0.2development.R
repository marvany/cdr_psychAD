library(data.table)
library(openxlsx)
library(readxl)
library(data.table)
library(ggplot2)
library(dplyr)
library(forcats)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(scales)
library(tidyr)

outwrite2 <- function(x,filename){fwrite(x,file.path("results/processed2", filename))}



bcmarks = read_excel('/sc/arion/projects/va-biobank/PROJECTS/cdr.comparative.efficacy.marios/Resources/benchmark_drugs.xlsx')
bc = as.data.table(bcmarks)
ad.bc = bc[comb_indication == 'AD']
adbc = ad.bc$pert_iname


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

# Create analysis type column
dt[, analysis_type := fcase(
  !is.na(`Pr(>|z|)`), 'logistic',
  !is.na(`Pr(>|t|)`), 'linear',
  default = NA
)]

# Calculate FDR adjusted p-values
dt[, fdr_adjusted := p.adjust(p_value, method = 'fdr')]

# Add significance flags
dt[, significant := fdr_adjusted < 0.05]

# ============= 1. Count FDR-significant compounds per phenotype =============

# Count number of significant compounds per phenotype
phenotype_counts <- dt[significant == TRUE, .N, by = phenotype][order(-N)]
phenotype_counts[, phenotype := factor(phenotype, levels = phenotype)]

# Add a count of total tests per phenotype
total_tests_per_phenotype <- dt[, .N, by = phenotype]
setnames(total_tests_per_phenotype, "N", "total_tests")

# Merge to get both significant counts and total tests
phenotype_counts <- merge(phenotype_counts, total_tests_per_phenotype, by = "phenotype")
phenotype_counts[, percentage := 100 * N / total_tests]

# ============= 2. Create horizontal bar plot =============

# For better visualization, limit to phenotypes with at least one significant hit
phenotype_counts_filtered <- phenotype_counts[N > 0]

# Reorder by number of significant compounds
phenotype_counts_filtered[, phenotype := fct_reorder(phenotype, N)]

# Create bar plot
p1 <- ggplot(phenotype_counts_filtered, aes(x = phenotype, y = N)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = N), hjust = -0.2, size = 3) +
  coord_flip() +
  labs(
    title = "Number of FDR-Significant Compounds per Phenotype",
    subtitle = "FDR-adjusted p-value < 0.05",
    x = NULL,
    y = "Number of Significant Compounds"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 8)
  )

# Save plot
ggsave("fdr_significant_compounds_per_phenotype.pdf", p1, width = 10, height = max(8, nrow(phenotype_counts_filtered) * 0.25))

# ============= 3. Create percentage bar plot =============

p2 <- ggplot(phenotype_counts_filtered, aes(x = phenotype, y = percentage)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), hjust = -0.2, size = 3) +
  coord_flip() +
  labs(
    title = "Percentage of FDR-Significant Compounds per Phenotype",
    subtitle = "FDR-adjusted p-value < 0.05",
    x = NULL,
    y = "Percentage of Tested Compounds (%)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 8)
  )

# Save plot
ggsave("percentage_fdr_significant_compounds_per_phenotype.pdf", p2, width = 10, height = max(8, nrow(phenotype_counts_filtered) * 0.25))

# ============= 4. Create bubble plot =============

# Bubble plot combining count and percentage
p3 <- ggplot(phenotype_counts_filtered, 
            aes(x = N, y = reorder(phenotype, N), size = percentage, color = percentage)) +
  geom_point() +
  scale_size_continuous(name = "% Significant", range = c(2, 10)) +
  scale_color_viridis(name = "% Significant") +
  labs(
    title = "FDR-Significant Compounds per Phenotype",
    subtitle = "Size and color represent percentage of tested compounds",
    x = "Number of Significant Compounds",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

# Save plot
ggsave("bubble_fdr_significant_compounds_per_phenotype.pdf", p3, width = 10, height = max(8, nrow(phenotype_counts_filtered) * 0.25))

# ============= 5. Extract the top compounds for the top phenotypes =============

# Identify top phenotypes with significant compounds (top 5 or all if fewer)
top_n <- min(5, nrow(phenotype_counts_filtered))
top_phenotypes <- phenotype_counts_filtered[1:top_n]$phenotype

# Get top significant compounds for each top phenotype
top_compounds_per_phenotype <- dt[significant == TRUE & phenotype %in% top_phenotypes, 
                                 .(compound, phenotype, Estimate, p_value, fdr_adjusted)]

# For each phenotype, find top 5 compounds (or all if fewer)
top_compounds_list <- list()
for (pheno in top_phenotypes) {
  pheno_compounds <- top_compounds_per_phenotype[phenotype == pheno][order(fdr_adjusted)]
  top_n_compounds <- min(5, nrow(pheno_compounds))
  top_compounds_list[[pheno]] <- pheno_compounds[1:top_n_compounds]
}

top_compounds_df <- rbindlist(top_compounds_list)

# Create a summary table
write.csv(top_compounds_df, "top_compounds_for_top_phenotypes.csv", row.names = FALSE)

# ============= 6. Create a heat tile plot for phenotype-analysis combinations =============

# Count by phenotype and analysis type
phenotype_analysis_counts <- dt[significant == TRUE, .N, by = .(phenotype, analysis_type)][order(-N)]

# Create a wide format for the heat tile
phenotype_analysis_wide <- dcast(phenotype_analysis_counts, phenotype ~ analysis_type, value.var = "N", fill = 0)

# Convert back to long format for plotting
phenotype_analysis_long <- melt(phenotype_analysis_wide, 
                               id.vars = "phenotype", 
                               variable.name = "analysis_type", 
                               value.name = "count")

# Only include phenotypes with at least one significant hit
phenotype_analysis_long <- phenotype_analysis_long[phenotype %in% phenotype_counts_filtered$phenotype]

# Reorder phenotypes by total count
phenotype_totals <- phenotype_analysis_long[, .(total = sum(count)), by = phenotype]
phenotype_order <- phenotype_totals[order(-total)]$phenotype
phenotype_analysis_long[, phenotype := factor(phenotype, levels = phenotype_order)]

# Create heat tile plot
p4 <- ggplot(phenotype_analysis_long, aes(x = analysis_type, y = phenotype, fill = count)) +
  geom_tile() +
  geom_text(aes(label = ifelse(count > 0, count, "")), color = "black", size = 3) +
  scale_fill_viridis(name = "Count", option = "plasma") +
  labs(
    title = "FDR-Significant Compounds by Phenotype and Analysis Type",
    x = "Analysis Type",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

# Save plot
ggsave("heatmap_phenotype_analysis_type.pdf", p4, width = 10, height = max(8, nrow(phenotype_counts_filtered) * 0.25))

# ============= 7. Create a combined visualization =============

# Combine plots using patchwork
combined_plot <- p1 / p3
ggsave("combined_phenotype_significant_compounds.pdf", combined_plot, width = 12, height = 12)

# ============= 8. Create a summary stats table =============

summary_stats <- data.table(
  Metric = c(
    "Total unique phenotypes",
    "Phenotypes with at least one significant compound",
    "Total compounds tested",
    "Total significant compound-phenotype associations",
    "Phenotype with most significant compounds",
    "Highest percentage of significant compounds"
  ),
  Value = c(
    length(unique(dt$phenotype)),
    nrow(phenotype_counts_filtered),
    length(unique(dt$compound)),
    sum(dt$significant),
    as.character(phenotype_counts_filtered$phenotype[1]),
    sprintf("%.1f%% (%s)", 
            phenotype_counts_filtered[which.max(percentage)]$percentage,
            as.character(phenotype_counts_filtered[which.max(percentage)]$phenotype))
  )
)

write.csv(summary_stats, "phenotype_significant_compounds_summary.csv", row.names = FALSE)

# Generate an HTML report that shows all visualizations
html_output <- '<!DOCTYPE html>
<html>
<head>
<title>FDR-Significant Compounds per Phenotype</title>
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
<h1>FDR-Significant Compounds per Phenotype Analysis</h1>

<div class="summary">
<h2>Summary Statistics</h2>
<table>
<tr><th>Metric</th><th>Value</th></tr>'

# Add each row from summary_stats
for (i in 1:nrow(summary_stats)) {
  html_output <- paste0(html_output, 
                       '<tr><td>', summary_stats$Metric[i], '</td><td>', summary_stats$Value[i], '</td></tr>')
}

html_output <- paste0(html_output, '</table>
</div>

<div class="viz">
<h2>Number of FDR-Significant Compounds per Phenotype</h2>
<p>This bar chart shows the count of compounds that have a significant association (FDR-adjusted p-value < 0.05) with each phenotype.</p>
<img src="fdr_significant_compounds_per_phenotype.png" alt="Bar chart of significant compounds per phenotype">
</div>

<div class="viz">
<h2>Percentage of FDR-Significant Compounds per Phenotype</h2>
<p>This bar chart shows the percentage of tested compounds that have a significant association with each phenotype.</p>
<img src="percentage_fdr_significant_compounds_per_phenotype.png" alt="Bar chart of percentage of significant compounds per phenotype">
</div>

<div class="viz">
<h2>Bubble Plot of FDR-Significant Compounds per Phenotype</h2>
<p>This bubble plot combines the count and percentage of significant compounds. The size and color of each bubble represent the percentage of tested compounds that are significant.</p>
<img src="bubble_fdr_significant_compounds_per_phenotype.png" alt="Bubble plot of significant compounds per phenotype">
</div>

<div class="viz">
<h2>Heatmap of Significant Compounds by Phenotype and Analysis Type</h2>
<p>This heatmap shows the count of significant compounds for each phenotype, broken down by analysis type (linear or logistic).</p>
<img src="heatmap_phenotype_analysis_type.png" alt="Heatmap of phenotype and analysis type">
</div>

<div class="viz">
<h2>Top Compounds for Top Phenotypes</h2>
<p>The table below shows the top 5 most significant compounds for each of the top phenotypes.</p>
<table>
<tr><th>Phenotype</th><th>Compound</th><th>Effect Size</th><th>P-value</th><th>FDR-adjusted P</th></tr>')

# Check if top_compounds_df exists and has rows
if (exists("top_compounds_df") && nrow(top_compounds_df) > 0) {
  for (i in 1:nrow(top_compounds_df)) {
    html_output <- paste0(html_output, 
                         '<tr>',
                         '<td>', top_compounds_df$phenotype[i], '</td>',
                         '<td>', top_compounds_df$compound[i], '</td>',
                         '<td>', round(top_compounds_df$Estimate[i], 4), '</td>',
                         '<td>', format(top_compounds_df$p_value[i], scientific = TRUE, digits = 3), '</td>',
                         '<td>', format(top_compounds_df$fdr_adjusted[i], scientific = TRUE, digits = 3), '</td>',
                         '</tr>')
  }
}

html_output <- paste0(html_output, '</table>
<p><a href="top_compounds_for_top_phenotypes.csv">Download complete table (CSV)</a></p>
</div>

<footer>
<p>Generated on ', format(Sys.time(), "%Y-%m-%d"), '</p>
</footer>
</div>
</body>
</html>')

# Write HTML file
writeLines(html_output, "fdr_significant_compounds_per_phenotype.html")

# Convert PDFs to PNGs for HTML display
pdf_files <- c("fdr_significant_compounds_per_phenotype.pdf", 
               "percentage_fdr_significant_compounds_per_phenotype.pdf", 
               "bubble_fdr_significant_compounds_per_phenotype.pdf",
               "heatmap_phenotype_analysis_type.pdf")

# Check if ImageMagick is available
has_convert <- system("which convert", intern = TRUE)
if(length(has_convert) > 0) {
  for(pdf_file in pdf_files) {
    png_file <- sub("\\.pdf$", ".png", pdf_file)
    system(paste("convert -density 150", pdf_file, png_file))
  }
}

cat("Analysis complete. Visualizations of FDR-significant compounds per phenotype are available.\n")
cat("Open 'fdr_significant_compounds_per_phenotype.html' to view a summary of all visualizations.\n")