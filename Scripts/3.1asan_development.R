library(data.table)
library(ordinal)
library(MASS)
library(parallel)
library(pbmcapply)

create_coefficients_df = function(fit.lm, phen, cmp, anl) {
  # Extract and organize coefficients
  coef_summary <- coef(summary(fit.lm))[-1, , drop=FALSE]
  coefficients_df <- data.table(t(as.data.frame(coef_summary)))
  setnames(coefficients_df, names(coefficients_df), rownames(coef_summary))
  
  # Add metadata
  coefficients_df[, `:=`(
    phenotype = phen,
    compound = cmp,
    analysis = anl
  )]
  
  # Set column names based on analysis type
  if(anl == 'logistic') {
    stat_col <- "z value"
    p_col <- "Pr(>|z|)"
  } else {
    stat_col <- "t value"
    p_col <- "Pr(>|t|)"
  }
  
  # Calculate confidence intervals
  estimates <- coef_summary[, "Estimate"]
  std_errors <- coef_summary[, "Std. Error"]
  t_value <- qt(0.975, df.residual(fit.lm))
  
  # Add all statistics
  coefficients_df[, `:=`(
    conf_int_low = estimates - t_value * std_errors,
    conf_int_high = estimates + t_value * std_errors,
    fdr.adjusted = p.adjust(get(p_col), method = "fdr"),
    bonferroni.adjusted = p.adjust(get(p_col), method = "bonferroni"),
    BH.adjusted = p.adjust(get(p_col), method = "BH"),
    multiple_r_squared = summary(fit.lm)$r.squared,
    adjusted_r_squared = summary(fit.lm)$adj.r.squared,
    AIC = AIC(fit.lm)
  )]
  
  # Add deviance statistics - these are model-dependent
  model_class <- class(fit.lm)[1]
  rse <- summary(fit.lm)$sigma
  freedom_degrees <- df.residual(fit.lm)
  
  if(!is.null(rse)) {
    coefficients_df[, deviance := rse^2 * freedom_degrees]
  }
  
  if(model_class == "glm") {
    coefficients_df[, `:=`(
      glm.deviance = fit.lm$deviance,
      null_deviance = fit.lm$null.deviance
    )]
  }
  
  # Add residual summary statistics
  deviance_residuals <- residuals(fit.lm, type = "deviance")
  coefficients_df[, `:=`(
    min_residual_se = min(deviance_residuals),
    first_q_residual_se = quantile(deviance_residuals, 0.25),
    median_residual_se = median(deviance_residuals),
    mean_residual_se = mean(deviance_residuals),
    third_q_residual_se = quantile(deviance_residuals, 0.75),
    max_residual_se = max(deviance_residuals)
  )]
  
  # Calculate variance partitioning
  var_part <- tryCatch({
    calcVarPart(fit.lm)
  }, error = function(e) {
    c(NA_real_, NA_real_)
  })
  
  if(length(var_part) > 0) {
    coefficients_df[, `:=`(
      varPart = var_part[1],
      residuals = var_part[2]
    )]
  }
  
  # Add number of observations
  coefficients_df[, n_observations := nrow(fit.lm$model)]
  
  return(coefficients_df)
}

# Function to determine analysis type based on variable characteristics
determine_analysis_type <- function(data, variable_name) {
  # Remove NA values for this variable
  valid_data <- data[!is.na(get(variable_name))]
  
  # Get unique values
  unique_vals <- unique(valid_data[[variable_name]])
  n_unique <- length(unique_vals)
  
  # If only 2 unique values, use logistic regression
  if(n_unique == 2) {
    return("logistic")
  } 
  # If less than 7 unique values, use ordinal regression
  else if(n_unique < 7 && all(unique_vals %in% c(0:10))) { 
    return("ordinal")
  }
  # Otherwise, use linear regression
  else {
    return("linear")
  }
}

# Function to standardize a variable (scale)
standardize_variable <- function(x) {
  if(is.numeric(x)) {
    return(as.numeric(scale(x)))
  } else {
    return(x)
  }
}

# Source necessary scripts
  source("/sc/arion/projects/va-biobank/PROJECTS/ma_neurocog_grex/Scripts/GReX_to_pheno/accessoryScripts.R")
  p.stat_calcVarPart = "/sc/arion/projects/va-biobank/PROJECTS/ma_neurocog_grex/Scripts/GReX_to_pheno/stat_calcVarPart.R"
  p.stat_PseudoR2 = "/sc/arion/projects/va-biobank/PROJECTS/ma_neurocog_grex/Scripts/GReX_to_pheno/stat_PseudoR2.R"
  source(p.stat_calcVarPart)
  source(p.stat_PseudoR2)
  setwd("/sc/arion/projects/va-biobank/PROJECTS/ma_cdr_psychAD")

##############
# DEFINE PATHS
p.drugs.phenos = 'Resources/result_analysis/cdr_psychAD_IC_microglia_v4.csv'
p.psychAD = './Resources/psychAD/clinical_metadata.csv'
p.psychAD.IDmap = '/sc/arion/projects/roussp01a/deepika/merging_psychAD_SNParray_WGS/common_variants_psychAD/ancestry_pca_psychAD_1429_samples/psychAD_20PC_3_methods_ancestry.tsv'




##################
# MUNGE PHENOTYPES

# We will examine which columns that correspond to phenotypes need to be filtered out
# Confirm that all phenotypes are present
psychAD = fread(p.psychAD)


# We keep only numeric, logical and columns 'ApoE_gt', 'Sex', 'Ethnicity'
numeric_or_logical = sapply(psychAD, function(x) (is.numeric(x) | is.logical(x)))
psychAD = psychAD[, c(names(psychAD)[numeric_or_logical], c('ApoE_gt', 'Sex', 'Ethnicity')), with = FALSE]



# Load all the variables
source('/sc/arion/projects/va-biobank/PROJECTS/ma_cdr_psychAD/Scripts/clinical_variables.R')


# Variables with which we will associate
association.vars = c(
  alzheimer_vars,
  cdr_vars,
  plaque_vars,
  tangle_vars,
  gliosis_vars
)

# Variables with which we will exclude
exclude.based.on.vars = c(
  neuro_vars,
  psych_vars
)

# Variables for initiall filtering
all.vars = c(general, association.vars, exclude.based.on.vars)

###################
# MUNGE CDR RESULTS (attention you are munging the aggregated results from the jupyter notebooks part)
# LOAD MERGED DATAFRAME (contains combined psychAD and CDR results; created in Jupyter Notebook)
df.merged = fread(p.drugs.phenos) # everything will be sourced from here
