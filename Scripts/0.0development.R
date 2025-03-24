library(data.table)
library(ggplot2)
library(ggrepel)

# LOAD TABLE
all.results = fread("/sc/arion/projects/va-biobank/PROJECTS/ma_cdr_psychAD/results/drug_pheno_assoc_anal_psychAD.csv")

 # Extract adjusted p-values
if(anl == 'logistic'){
      coefficients_df$fdr.adjusted = p.adjust(coefficients_df$`Pr(>|z|)`, method = "fdr")
      coefficients_df$bonferroni.adjusted = p.adjust(coefficients_df$`Pr(>|z|)`, method = "bonferroni")
      coefficients_df$BH.adjusted = p.adjust(coefficients_df$`Pr(>|z|)`, method = "BH")
} else{
    coefficients_df$fdr.adjusted = p.adjust(coefficients_df$`Pr(>|t|)`, method = "fdr")
    coefficients_df$bonferroni.adjusted = p.adjust(coefficients_df$`Pr(>|t|)`, method = "bonferroni")
    coefficients_df$BH.adjusted = p.adjust(coefficients_df$`Pr(>|t|)`, method = "BH")
}


#fcase(
  # condition1, value1,
  # condition2, value2,
  # condition3, value3,
  #default = default_value
#)

# REMOVE NA VALUES
dt_all.cols = all.results[!is.na(Estimate),]


# CREATE PVALUE COLUMN
dt_all.cols[, p := fcase(
  !is.na(`Pr(>|z|)`), `Pr(>|z|)`,
  !is.na(`Pr(>|t|)`), `Pr(>|t|)`,
  default = NA
)]


# CREATE ANALYSIS COLUMN
dt_all.cols[, analysis := fcase(
  !is.na(`Pr(>|z|)`), 'logistic',
  !is.na(`Pr(>|t|)`), 'linear',
  default = NA
)]

# CALCULATE FDR & BONFERRONI
dt_all.cols[, fdr.adjusted := p.adjust(p, method = 'fdr')]
dt_all.cols[, bonferroni.adjusted := p.adjust(p, method = 'bonferroni')]
sum(dt_all.cols$fdr.adjusted < 0.05)
sum(dt_all.cols$bonferroni.adjusted < 0.05)

# 

