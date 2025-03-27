# Load required libraries
library(data.table)
library(ordinal)
library(MASS)
library(parallel)
library(pbmcapply)

# Source necessary scripts
source('/sc/arion/projects/va-biobank/PROJECTS/ma_cdr_psychAD/Scripts/clinical_variables.R')
source("/sc/arion/projects/va-biobank/PROJECTS/ma_neurocog_grex/Scripts/GReX_to_pheno/accessoryScripts.R")
p.stat_calcVarPart = "/sc/arion/projects/va-biobank/PROJECTS/ma_neurocog_grex/Scripts/GReX_to_pheno/stat_calcVarPart.R"
p.stat_PseudoR2 = "/sc/arion/projects/va-biobank/PROJECTS/ma_neurocog_grex/Scripts/GReX_to_pheno/stat_PseudoR2.R"
source(p.stat_calcVarPart)
source(p.stat_PseudoR2)
setwd("/sc/arion/projects/va-biobank/PROJECTS/ma_cdr_psychAD")

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

# Function to standardize variables
standardize_variable <- function(x) {
        if(is.numeric(x)) {
          return(as.numeric(scale(x)))
        } else {
          return(x)
        }
}

# Function to process regression results
create_coefficients_df = function(fit.lm, phen, cmp, anl){
  browser()
      # Extract coefficients based on model type
      model_type <- class(fit.lm)[1]
      is_clm <- inherits(fit.lm, "clm")
      
      # Extract coefficients with appropriate method
      coef_summary <- summary(fit.lm)$coefficients
      
      # Process differently based on model type
      if(is_clm) {
        compound_coef <- coef_summary[grepl(cmp, rownames(coef_summary)), , drop = FALSE]
        coefficients_df <- data.table(
              compound = cmp, phenotype = phen, analysis = anl,
              Estimate = compound_coef[1, "Estimate"],
              `Std. Error` = compound_coef[1, "Std. Error"],
              `z value` = compound_coef[1, "z value"],
              `Pr(>|z|)` = compound_coef[1, "Pr(>|z|)"]
        )
        
        # Add CI values and p-value adjustments
        critical_value <- 1.96 # Z distribution
        p_col <- "Pr(>|z|)"
      } else {
        # Regular glm/lm models
        coef_summary <- coef_summary[-1, , drop=FALSE]
        coefficients_df <- data.table(t(as.data.frame(coef_summary)))
        setnames(coefficients_df, names(coefficients_df), rownames(coef_summary))
        
        # Add basic metadata
        coefficients_df[, `:=`(phenotype = phen, compound = cmp, analysis = anl)]
        
        # Set analysis-specific values
        p_col <- if(anl == 'logistic') "Pr(>|z|)" else "Pr(>|t|)"
        critical_value <- qt(0.975, df.residual(fit.lm))

        coefficients_df[, `:=`(
                deviance_explained := 1 - (fit.lm$deviance / fit.lm$null.deviance)
                
                #multiple_r_squared = summary(fit.lm)$r.squared,
                #adjusted_r_squared = summary(fit.lm)$adj.r.squared
        )]

        if(!is.null(rse <- summary(fit.lm)$sigma)) {
              coefficients_df[, deviance := rse^2 * df.residual(fit.lm)]
            }
            
        if(inherits(fit.lm, "glm")) {
          coefficients_df[, `:=`(
                glm.deviance = fit.lm$deviance,
                null_deviance = fit.lm$null.deviance
          )]
        }
        
        # Add residual statistics for non-clm models
        deviance_residuals <- residuals(fit.lm, type = "deviance")
        coefficients_df[, `:=`(
              min_residual_se = min(deviance_residuals),
              first_q_residual_se = quantile(deviance_residuals, 0.25),
              median_residual_se = median(deviance_residuals),
              mean_residual_se = mean(deviance_residuals),
              third_q_residual_se = quantile(deviance_residuals, 0.75),
              max_residual_se = max(deviance_residuals)
        )]
        
        # Try to calculate variance partitioning for non-clm models
        tryCatch({
          var_part <- calcVarPart(fit.lm)
          if(length(var_part) > 0) {
            coefficients_df[, `:=`(varPart = var_part[1], residuals = var_part[2])]
          }
        }, error = function(e) {})

      }
      
      # Add common statistics for all model types
      coefficients_df[, `:=`(
              conf_int_low = Estimate - critical_value * `Std. Error`,
              conf_int_high = Estimate + critical_value * `Std. Error`,
              AIC = AIC(fit.lm),                  # <- Missing comma here
              n_observations = nrow(fit.lm$model)
      )]
      
      return(coefficients_df)
}


#########################################
# NOTES: Here you run the glm ; maybe better create diff function
# for each of the lm , lg , clm
# 

# Define paths
p.drugs.phenos = 'Resources/result_analysis/cdr_psychAD_IC_microglia_v4.csv'
p.psychAD = './Resources/psychAD/clinical_metadata.csv'
p.psychAD.IDmap = '/sc/arion/projects/roussp01a/deepika/merging_psychAD_SNParray_WGS/common_variants_psychAD/ancestry_pca_psychAD_1429_samples/psychAD_20PC_3_methods_ancestry.tsv'

outdir = paste0(getwd(), '/results/assoc_anal/v1')

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


# Load and prepare the data
message("Loading merged drug-phenotype data...")
df.merged <- fread(p.drugs.phenos)

# Load psychAD data for reference
psychAD <- fread(p.psychAD)

# Convert necessary columns
df.merged[, ApoE_gt := as.numeric(ApoE_gt)]
df.merged[, Sex := ifelse(Sex == 'Male', 1, 0)]
df.merged[, Ethnicity := NULL]

# Load PCs data for covariates
pcs_data <- fread(p.psychAD.IDmap)
pc_columns <- paste0("PC", 1:10)
pcs_data <- pcs_data[, c("SubID", ..pc_columns)]
setnames(pcs_data, "SubID", "ind_ID_clin")

# Merge PC data with main dataset
df.merged <- merge(df.merged, pcs_data, by = "ind_ID_clin", all.x = TRUE)

# Filter out individuals with neurological or psychiatric conditions
message("Filtering out individuals with neurological or psychiatric conditions...")
for(var in c(neuro_vars, psych_vars)) {
  if(var %in% names(df.merged)) {
    df.merged <- df.merged[get(var) == 0 | is.na(get(var))]
  }
}

message(paste("After filtering neuropsych conditions:", nrow(df.merged), "individuals remain"))

# Get rid of '-' to avoid errors
colnames(df.merged) = gsub('-', '_', colnames(df.merged))
psychAD_cols = gsub('-', '_', colnames(psychAD))

# Identify variable sets for analysis
phenotype_variables <- intersect(names(df.merged), psychAD_cols)
phenotype_variables <- setdiff(phenotype_variables, c(neuro_vars, psych_vars, "Sex", "Ethnicity", "ind_ID", "ind_ID_clin"))

# Identify drug columns (everything that's not a phenotype or ID or PC)
drug_columns <- setdiff(names(df.merged), c(phenotype_variables, 'ind_ID', 'ind_ID_clin', pc_columns))

# Only keep association variables that we want to test
phenotype_variables <- intersect(phenotype_variables, c('', association.vars))

# Determine analysis type for each phenotype
analysis_types <- sapply(phenotype_variables, function(var) {
  determine_analysis_type(df.merged, var)
})

# Map variables to analysis types for reference
var_to_analysis <- setNames(analysis_types, phenotype_variables)

# Display analysis summary
analysis_summary <- data.table(
  phenotype = names(var_to_analysis),
  analysis_type = unname(var_to_analysis)
)
print(analysis_summary)

# Create combo grid for drug-phenotype pairs
message("Creating combination grid for drug-phenotype pairs...")
combo_grid <- CJ(compound = drug_columns, phenotype = phenotype_variables, sorted = FALSE) # CJ creates data.table with all combinations of provided vectors
combo_grid[, analysis := var_to_analysis[phenotype]]


# Create age squared column
df.merged[, age_sq := df.merged$Age^2]

# Scale numeric variables for analysis
message("Scaling variables for analysis...")
# Always scale drug_columns, pc_columns, and Age
always_scale <- c(drug_columns, pc_columns, "Age", "age_sq")
# Only scale phenotype variables getting linear regression
linear_phenos <- names(var_to_analysis[var_to_analysis == "linear"])
scale_cols <- c(always_scale, linear_phenos)
numeric_cols <- scale_cols[scale_cols != "Sex" & 
                           sapply(df.merged[, ..scale_cols], is.numeric)]
df.merged[, (numeric_cols) := lapply(.SD, scale), .SDcols = numeric_cols]



# Create directory for results if it doesn't exist
if(!dir.exists(outdir)) {
  dir.create(outdir)
}

# Perform association analysis
message("Running associations...")

#results <- rbindlist(pbmclapply(1:nrow(combo_grid), function(i) {
results <- rbindlist(lapply(1:nrow(combo_grid), function(i) {
  thisrow <- combo_grid[i]
  cmp <- thisrow$compound
  phen <- thisrow$phenotype
  anl <- thisrow$analysis
  message(i)
  # Prepare data and formula
  model_vars <- c(cmp, phen, pc_columns, "Age", "age_sq", "Sex")
  temp <- na.omit(df.merged[, ..model_vars])
  
  if(i == 2)browser()
  if(nrow(temp) <= length(model_vars) + 5) return(NULL)
  
  formula_str <- paste(phen, "~", cmp, "+", paste(c(pc_columns, "Age", "age_sq", "Sex"), collapse = " + "))
  # Run model based on analysis type
  if(anl == "linear") {
    fit.lm <- glm(as.formula(formula_str), data = temp, family = gaussian())
    return(create_coefficients_df(fit.lm, phen, cmp, anl))
    
  } else if(anl == "logistic") {
    temp[[phen]] <- as.factor(temp[[phen]])
    fit.lm <- glm(as.formula(formula_str), data = temp, family = binomial())
    return(create_coefficients_df(fit.lm, phen, cmp, anl))
    
  } else if(anl == "ordinal") {
    temp[[phen]] <- factor(temp[[phen]], ordered = TRUE)
    tryCatch({
      fit.clm <- clm(as.formula(formula_str), data = temp, link = "logit")
      return(create_coefficients_df(fit.clm, phen, cmp, anl))
    }, error = function(e) return(NULL))
  }
  
  return(NULL)
}))
#}, mc.cores = detectCores() - 2), fill = TRUE)


# Clean up results
results <- results[!is.na(compound)]

# Compute unified p-value column for easier interpretation
results[, p_value := fcase(
  !is.na(`Pr(>|z|)`), `Pr(>|z|)`,
  !is.na(`Pr(>|t|)`), `Pr(>|t|)`,
  default = NA
)]

# Export results
fwrite(results, 'results/drug_pheno_assoc_anal_psychAD_final.csv')

# Summary statistics
message("Analysis complete!")
message(paste("Total associations tested:", nrow(combo_grid)))
message(paste("Significant associations (FDR < 0.05):", sum(results$fdr.adjusted < 0.05, na.rm = TRUE)))
message(paste("Significant associations (Bonferroni < 0.05):", sum(results$bonferroni.adjusted < 0.05, na.rm = TRUE)))

# Top significant hits
top_hits <- results[order(p_value)][1:10]
print(top_hits[, .(compound, phenotype, analysis, Estimate, p_value, fdr.adjusted)])