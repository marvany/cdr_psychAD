# Environment setup
{
  library(data.table)
  library(ordinal)
  library(MASS)
  library(parallel)
  library(pbmcapply)
  
  # Source necessary scripts
  source("/sc/arion/projects/va-biobank/PROJECTS/ma_neurocog_grex/Scripts/GReX_to_pheno/accessoryScripts.R")
  p.stat_calcVarPart = "/sc/arion/projects/va-biobank/PROJECTS/ma_neurocog_grex/Scripts/GReX_to_pheno/stat_calcVarPart.R"
  p.stat_PseudoR2 = "/sc/arion/projects/va-biobank/PROJECTS/ma_neurocog_grex/Scripts/GReX_to_pheno/stat_PseudoR2.R"
  source(p.stat_calcVarPart)
  source(p.stat_PseudoR2)
  setwd("/sc/arion/projects/va-biobank/PROJECTS/ma_cdr_psychAD")
  
  # Helper functions
  create_coefficients_df = function(fit, fit.lm, phen, cmp, anl){
    # Extract coefficients
    coefs = as.data.frame(coef(summary(fit.lm))[-1, ])
    temp2 = transpose(coefs)
    colnames(temp2) = rownames(coefs)
    coefficients_df = temp2
    
    # Add metadata
    coefficients_df$phenotype = phen
    coefficients_df$compound = cmp
    
    # Format based on analysis type
    if(anl == 'logistic'){
      coefficients_df = coefficients_df[, c("compound", "phenotype", "Estimate", "Std. Error", "z value", "Pr(>|z|)")]
    } else {
      coefficients_df = coefficients_df[, c("compound", "phenotype", "Estimate", "Std. Error", "t value", "Pr(>|t|)")]
    }
    
    # Extract coefficients and standard errors for CI calculation
    coeff_summary <- coef(summary(fit.lm))
    estimates <- coeff_summary[-1, "Estimate"]
    std_errors <- coeff_summary[-1, "Std. Error"]
    
    # Calculate confidence intervals
    t_value <- qt(0.975, df.residual(fit.lm))
    coefficients_df$conf_int_low <- estimates - t_value * std_errors
    coefficients_df$conf_int_high <- estimates + t_value * std_errors
    
    # Add adjusted p-values
    if(anl == 'logistic'){
      coefficients_df$fdr.adjusted = p.adjust(coefficients_df$`Pr(>|z|)`, method = "fdr")
      coefficients_df$bonferroni.adjusted = p.adjust(coefficients_df$`Pr(>|z|)`, method = "bonferroni")
      coefficients_df$BH.adjusted = p.adjust(coefficients_df$`Pr(>|z|)`, method = "BH")
    } else{
      coefficients_df$fdr.adjusted = p.adjust(coefficients_df$`Pr(>|t|)`, method = "fdr")
      coefficients_df$bonferroni.adjusted = p.adjust(coefficients_df$`Pr(>|t|)`, method = "bonferroni")
      coefficients_df$BH.adjusted = p.adjust(coefficients_df$`Pr(>|t|)`, method = "BH")
    }
    
    # Extract model statistics
    coefficients_df$multiple_r_squared <- summary(fit.lm)$r.squared
    coefficients_df$adjusted_r_squared <- summary(fit.lm)$adj.r.squared
    
    # Extract Residual Deviance
    rse = summary(fit.lm)$sigma
    freedom_degrees = df.residual(fit.lm)
    if(!is.null(rse)) coefficients_df$deviance = rse^2 * freedom_degrees
    coefficients_df$glm.deviance = fit$deviance
    
    # Extract Null Deviance
    coefficients_df$null_deviance <- fit$null.deviance
    
    # Extract AIC
    coefficients_df$AIC <- AIC(fit.lm)
    
    # Extract residual standard error summary
    deviance_residuals = residuals(fit.lm, type = "deviance")
    coefficients_df$min_residual_se = min(deviance_residuals) 
    coefficients_df$first_q_residual_se = quantile(deviance_residuals, 0.25)
    coefficients_df$median_residual_se = median(deviance_residuals)
    coefficients_df$mean_residual_se = mean(deviance_residuals)  
    coefficients_df$third_q_residual_se = quantile(deviance_residuals, 0.75)
    coefficients_df$max_residual_se = max(deviance_residuals)
    
    # Calculate variance partitioning
    if(anl == 'linear') var_part = calcVarPart(fit) else var_part = calcVarPart(fit.lm)
    coefficients_df$varPart <- var_part[1]
    coefficients_df$residuals <- var_part[2]
    
    # Add number of observations
    coefficients_df$n_observations <- nrow(fit.lm$model)
    
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
}

##############
# DEFINE PATHS
p.drugs.phenos = 'Resources/result_analysis/cdr_psychAD_IC_microglia_v4.csv'
p.psychAD = './Resources/psychAD/clinical_metadata.csv'
p.psychAD.IDmap = '/sc/arion/projects/roussp01a/deepika/merging_psychAD_SNParray_WGS/common_variants_psychAD/ancestry_pca_psychAD_1429_samples/psychAD_20PC_3_methods_ancestry.tsv'

# Define variables of interest
neuro_vars <- c("AD", "MCI", "Dementia", "PD", "PD_uncertain_plus_encephalitic", 
                "DLBD", "FTD", "ALS", "NormPressHydrocephalus")
psych_vars <- c("SCZ", "MDD", "BD_unspecific", "BD_I", "BD_II", "PTSD", "ADHD", "OCD")

##################
# MUNGE PHENOTYPES
message("Loading psychAD data...")
psychAD <- fread(p.psychAD)

# Keep only numeric, logical and specified categorical columns
numeric_or_logical <- sapply(psychAD, function(x) (is.numeric(x) | is.logical(x)))
psychAD <- psychAD[, c(which(numeric_or_logical), which(names(psychAD) %in% c('ApoE_gt', 'Sex', 'Ethnicity'))), with = FALSE]

# Get the phenotype column names
psychAD_cols <- names(psychAD)

# Load PCs for covariates
pcs_data <- fread(p.psychAD.IDmap)
pcs_columns <- paste0("PC", 1:10)
pcs_data <- pcs_data[, c("SubID", ..pcs_columns)]
setnames(pcs_data, "SubID", "ind_ID_clin")

###################
# LOAD MERGED DATA
message("Loading merged drug-phenotype data...")
df.merged <- fread(p.drugs.phenos)

# Format specific columns
df.merged[, ApoE_gt := as.numeric(ApoE_gt)]
df.merged[, Sex := ifelse(Sex == 'Male', 1, 0)]
df.merged[, Ethnicity := NULL]

# Add PC covariates to merged data
df.merged <- merge(df.merged, pcs_data, by = "ind_ID_clin", all.x = TRUE)

# Identify numeric and binary phenotypes for analysis type
phenotype_vars <- intersect(names(df.merged), psychAD_cols)
neuropsych_vars <- c(neuro_vars, psych_vars)

# Filter out individuals with neurological or psychiatric conditions
for(var in neuropsych_vars) {
  if(var %in% names(df.merged)) {
    df.merged <- df.merged[get(var) == 0 | is.na(get(var))]
  }
}

message(paste("After filtering neuropsych conditions:", nrow(df.merged), "individuals remain"))

# Get association variables (phenotypes we want to test)
association_vars <- setdiff(phenotype_vars, c(neuropsych_vars, "Sex", "Ethnicity", "ind_ID", "ind_ID_clin"))

# Identify drug columns
drug_cols <- setdiff(names(df.merged), c('ind_ID', 'ind_ID_clin', phenotype_vars, pcs_columns))

# Determine analysis type for each phenotype
analysis_types <- sapply(association_vars, function(var) {
  determine_analysis_type(df.merged, var)
})

# Create a mapping of variables to analysis types
var_to_analysis <- setNames(analysis_types, association_vars)

# Display analysis types
analysis_summary <- data.table(
  phenotype = names(var_to_analysis),
  analysis_type = unname(var_to_analysis)
)
print(analysis_summary)

# Create combo grid for drug-phenotype pairs
message("Creating combination grid for drug-phenotype pairs...")
combo_grid <- CJ(compound = drug_cols, phenotype = association_vars, sorted = FALSE)
combo_grid[, analysis := var_to_analysis[phenotype]]

# Scale numeric variables for analysis
message("Scaling variables for analysis...")
for(col in c(drug_cols, association_vars, pcs_columns, "age")) {
  if(col != "Sex" && is.numeric(df.merged[[col]])) {
    df.merged[, (col) := standardize_variable(get(col))]
  }
}

# Create age squared column
df.merged[, age_sq := df.merged$age^2]

###################
# PERFORM ANALYSIS
message("Running associations...")
results <- rbindlist(pbmclapply(1:nrow(combo_grid), function(i) {
  thisrow <- combo_grid[i]
  cmp <- thisrow$compound
  phen <- thisrow$phenotype
  anl <- thisrow$analysis
  
  message(sprintf("Processing %d/%d: %s ~ %s (%s)", i, nrow(combo_grid), phen, cmp, anl))
  
  # Select and prepare data
  model_vars <- c(cmp, phen, pcs_columns, "age", "age_sq", "Sex")
  temp <- na.omit(df.merged[, ..model_vars])
  
  # Skip if insufficient data
  if(nrow(temp) <= length(model_vars) + 5) {
    message("  Insufficient data points")
    return(NULL)
  }
  
  # Set up formula - use all covariates
  formula_str <- paste(phen, "~", cmp, "+", 
                      paste(c(pcs_columns, "age", "age_sq", "Sex"), collapse = " + "))
  
  # Run appropriate model
  if(anl == "linear") {
    # Linear Regression
    fit.lm <- lm(as.formula(formula_str), data = temp)
    fit <- glm(as.formula(formula_str), data = temp, family = gaussian(link = "identity"))
    dt <- create_coefficients_df(fit = fit, fit.lm = fit.lm, phen = phen, cmp = cmp, anl = anl)
    dt$n_phen <- nrow(temp)
    dt$analysis <- anl
    
  } else if(anl == "logistic") {
    # Ensure binary variable is factor
    temp[[phen]] <- as.factor(temp[[phen]])
    fit.lm <- glm(as.formula(formula_str), data = temp, family = binomial)
    dt <- create_coefficients_df(fit = fit.lm, fit.lm = fit.lm, phen = phen, cmp = cmp, anl = anl)
    dt$n_phen <- nrow(temp)
    dt$analysis <- anl
    
  } else if(anl == "ordinal") {
    # Ordinal Regression
    temp[[phen]] <- factor(temp[[phen]], ordered = TRUE)
    tryCatch({
      fit.clm <- clm(as.formula(formula_str), data = temp, link = "logit")
      # Extract coefficients (only for compound)
      coef_summ <- summary(fit.clm)$coefficients
      compound_coef <- coef_summ[grepl(cmp, rownames(coef_summ)), , drop = FALSE]
      
      dt <- data.table(
        compound = cmp,
        phenotype = phen,
        Estimate = compound_coef[1, "Estimate"],
        `Std. Error` = compound_coef[1, "Std. Error"],
        `z value` = compound_coef[1, "z value"],
        `Pr(>|z|)` = compound_coef[1, "Pr(>|z|)"]
      )
      
      # Add CIs
      dt$conf_int_low <- dt$Estimate - 1.96 * dt$`Std. Error`
      dt$conf_int_high <- dt$Estimate + 1.96 * dt$`Std. Error`
      
      # Add adjusted p-values
      dt$fdr.adjusted <- p.adjust(dt$`Pr(>|z|)`, method = "fdr")
      dt$bonferroni.adjusted <- p.adjust(dt$`Pr(>|z|)`, method = "bonferroni")
      dt$BH.adjusted <- p.adjust(dt$`Pr(>|z|)`, method = "BH")
      
      dt$n_phen <- nrow(temp)
      dt$analysis <- anl
      
    }, error = function(e) {
      message("  Error in ordinal regression: ", e$message)
      return(NULL)
    })
  }
  
  return(dt)
}, mc.cores = detectCores() - 2), fill = TRUE)

# Clean up results
results <- results[!is.na(compound)]

# Compute unified p-value column for easier interpretation
results[, p := fcase(
  !is.na(`Pr(>|z|)`), `Pr(>|z|)`,
  !is.na(`Pr(>|t|)`), `Pr(>|t|)`,
  default = NA
)]

# Export results
fwrite(results, 'results/drug_pheno_assoc_anal_psychAD_clean.csv')

# Summary statistics
message("Analysis complete!")
message(paste("Total associations tested:", nrow(combo_grid)))
message(paste("Significant associations (FDR < 0.05):", sum(results$fdr.adjusted < 0.05, na.rm = TRUE)))
message(paste("Significant associations (Bonferroni < 0.05):", sum(results$bonferroni.adjusted < 0.05, na.rm = TRUE)))

# Top significant hits
top_hits <- results[order(p)][1:10]
print(top_hits[, .(compound, phenotype, analysis, Estimate, p, fdr.adjusted)])