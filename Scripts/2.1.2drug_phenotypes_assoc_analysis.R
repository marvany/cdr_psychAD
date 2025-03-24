# Environment
{
      library(data.table)
      library(ordinal)
      library(MASS)
      library(parallel)

      source("/sc/arion/projects/va-biobank/PROJECTS/ma_neurocog_grex/Scripts/GReX_to_pheno/acc_lsf_bash_RFunctions.R")
      source("/sc/arion/projects/va-biobank/PROJECTS/ma_neurocog_grex/Scripts/GReX_to_pheno/acc_resubmit_general.R")
      source("/sc/arion/projects/va-biobank/PROJECTS/ma_neurocog_grex/Scripts/GReX_to_pheno/accessoryScripts.R")
      p.stat_calcVarPart = "/sc/arion/projects/va-biobank/PROJECTS/ma_neurocog_grex/Scripts/GReX_to_pheno/stat_calcVarPart.R"
      p.stat_PseudoR2 = "/sc/arion/projects/va-biobank/PROJECTS/ma_neurocog_grex/Scripts/GReX_to_pheno/stat_PseudoR2.R"
      source(p.stat_calcVarPart)
      source(p.stat_PseudoR2)
      setwd("/sc/arion/projects/va-biobank/PROJECTS/ma_cdr_psychAD")


      create_coeffiecients_df = function(fit, fit.lm, phen, cmp, anl){
        
            coefs <- coef(summary(fit.lm))[-1,]

            coefs = as.data.frame(coef(summary(fit.lm))[-1, ])
            temp2 = transpose(coefs)
            colnames(temp2) = rownames(coefs)
            coefficients_df = temp2

              coefficients_df$phenotype = phen
              coefficients_df$compound = cmp
              
              #coefficients_df$model_ID= basename(outdir) # edit
              if(anl == 'logistic'){
                coefficients_df = coefficients_df[, c("compound","phenotype", "Estimate", "Std. Error", "z value", "Pr(>|z|)")]
              } else coefficients_df = coefficients_df[, c("compound","phenotype", "Estimate", "Std. Error", "t value", "Pr(>|t|)")]
              
              # Extract coefficients and standard errors
              coeff_summary <- coef(summary(fit.lm))
              estimates <- coeff_summary[-1, "Estimate"]
              std_errors <- coeff_summary[-1, "Std. Error"]

              # Calculate the critical t-value for 95% confidence intervals
              # Degrees of freedom = residual degrees of freedom from the model
              t_value <- qt(0.975, df.residual(fit.lm))

              # Lower and upper bounds of the confidence intervals
              coefficients_df$conf_int_low <- estimates - t_value * std_errors
              coefficients_df$conf_int_high <- estimates + t_value * std_errors

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
              

              # Extract multiple R-squared
              coefficients_df$multiple_r_squared <- summary(fit.lm)$r.squared

              # Extract adjusted R-squared
              coefficients_df$adjusted_r_squared <- summary(fit.lm)$adj.r.squared

              # Extract Residual Deviance
              rse = summary(fit.lm)$sigma
              freedom_degrees = df.residual(fit.lm) # degrees of freedom
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
              coefficients_df$third_q_residual_se =quantile(deviance_residuals, 0.75)
              coefficients_df$max_residual_se = max(deviance_residuals)

              # Calculate variance partitioning
              if(anl == 'linear') var_part = calcVarPart(fit) else var_part = calcVarPart(fit.lm)
              coefficients_df$varPart <- var_part[1]
              coefficients_df$residuals <- var_part[2]

              return(coefficients_df)
      }

      make_names <- function(x) {
        x <- gsub("\\+", ".plus.", x)
        x <- gsub("-", ".minus.", x)
        x <- make.names(x)
        return(x)
      }

      WriteProperTSV <- function(object, file.path) {
        write.table(object, file = file.path,
                    quote = FALSE, sep = '\t', row.names = FALSE) }

      ## parse argument (parse_NA_T_F_list)
      parse_NA_T_F_list <- function(x) {
        x <- strsplit(x,",")[[1]]
        if (length(x) == 1) {
          if (x %in% c("NA", "TRUE", "T", "FALSE", "F")) {
            x <- eval(parse(text = x))
          }
        }
        return(x)
      }
}



##############
# DEFINE PATHS
#p.drugs.phenos = 'Resources/result_analysis/test_merged_drugs_phenotypes.csv'
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


# These are the phenotype columns
psychAD_cols = names(psychAD)


# YOU CAN SKIP THIS PART
# these columns were excluded; we leave them here for quality checking purposes
if(FALSE){
    excluded_cols = c('SubID', 'Brain_bank', 'Source_Location', 'Death_Time', 'AMPAD_msbb_individualID','CMC_individual_ID',
  'primary_genotype', 'SNParray_HBBC', 'SNParray_CommonMind', 'WGS_CommonMind', 'WGS_RUSH', 'SNParray_Microglia',
  'SNParray_PsychAD', 'ADSP_SampleId', 'Imaging_XENum', 'Dx', 'snRNAseq_ID', 'Death_Season')
  inds = which(!(names(example) %in% excluded_cols))
  example = psychAD[1:2,names(example)[inds], with = FALSE]
}


###################
# MUNGE CDR RESULTS (attention you are munging the aggregated results from the jupyter notebooks part)

# LOAD MERGED DATAFRAME (contains combined psychAD and CDR results; created in Jupyter Notebook)
df.merged = fread(p.drugs.phenos) # everything will be sourced from here

# Reshape columns that need reshaping (see above)
df.merged[,'ApoE_gt'] = as.numeric(df.merged[['ApoE_gt']])
df.merged[,'Sex'] = ifelse(df.merged[['Sex']] == 'Male', TRUE, FALSE)
df.merged[,'Ethnicity'] = NULL
# Filter out 'character' columns
#'ind_ID', 'ind_ID_clin', 
df.merged = df.merged[,c(names(df.merged)[which(sapply(df.merged, class) != 'character')]), with = FALSE]


# Identify numeric and binary phenotypes (for linear vs. logistic)
log.ind = sapply(df.merged, is.logical)
log.phen = names(df.merged)[log.ind]
log.phen = log.phen[log.phen %in% psychAD_cols]

num.ind = sapply(df.merged, is.numeric)
num.phen = names(df.merged)[num.ind]
num.phen = num.phen[num.phen %in% psychAD_cols]

df.merged = data.table(sapply(df.merged, function(x){
    # If the column has TRUE/FALSE ; if it's NA return NA ; if it's TRUE/FALSE return 1/0
    if(is.logical(x)){ ifelse(is.na(x), x, ifelse(x, 1, 0))
    # If nothing has been returned, x is numeric, so scale it
    }else scale(x)
}))


###################
# CREATE COMBO GRID

# These are the columns that cdr results of drugs
id_and_drugs_cols = setdiff(names(df.merged), c('ind_ID', 'ind_ID_clin', psychAD_cols))


# Create cdr/phenotype combo grid
num.combo.grid = expand.grid(compound = id_and_drugs_cols, phenotype = num.phen)
log.combo.grid = expand.grid(compound = id_and_drugs_cols, phenotype = log.phen)


# Iterate over combo.grid; store results in dataframe
log.combo.grid$analysis = 'logistic'
num.combo.grid$analysis = 'linear'
combo.grid = rbind(num.combo.grid, log.combo.grid)


##########
# ANALYSIS

# Turn '-' to '_'
names(df.merged) = gsub('-', '_', names(df.merged))
combo.grid = data.table(sapply(combo.grid, function(x) gsub('-', '_', x)))

###################
# MAKE THIS ITERATE
message('combo.grid is of class ', class(combo.grid))
message('df.merged is of class ', class(df.merged))

# for debugging
if(FALSE){
  combo.grid$phenotype = "multiome_DLPFC"
  combo.grid$compound = "BRD_A65076780"
  combo.grid$logistic = "logistic"
  combo.grid = combo.grid[1,]
}

results = rbindlist(pbmcapply::pbmclapply(seq(nrow(combo.grid)), function(i){
  
    thisrow = combo.grid[i,]
    cmp = as.character(thisrow$compound)
    phen = as.character(thisrow$phenotype)
    anl = thisrow$analysis
    # debug logistic: 
    #phen = 'nps_LateInsomHxValue'
    #anl = "logistic"
    
    temp = na.omit(df.merged[, c(cmp, phen), with = FALSE])
    temp$phen <- as.numeric(as.character(temp[[phen]]))
    temp$cmp <- as.numeric(as.character(temp[[cmp]]))
    
    message('phenotype is ', phen, ' compound is ', cmp, ' analysis is ', anl)
    if(nrow(temp) > 1){
            if (anl == "linear") {
              # Linear Regression
              fit.lm <- lm(phen ~ cmp, data = temp) 
              fit <- glm(phen ~ cmp, data = temp, family=gaussian(link = "identity")) 
              dt = create_coeffiecients_df(fit = fit, fit.lm = fit.lm, phen = phen, cmp = cmp, anl = anl)
              dt$n_phen = nrow(temp)
              dt$analysis = anl 
              message('the table produced has row_col dimentions ', paste(dim(dt), collapse = ' '))
              return(dt)

          }else if(anl == "logistic"){

              # Ensure phen is a binary variable (as a factor)
              temp$phen <- as.factor(temp$phen)
              fit.lm <- glm(phen ~ cmp, data = temp, family = binomial) 
              dt = create_coeffiecients_df(fit = fit.lm, fit.lm = fit.lm, phen = phen, cmp = cmp, anl = anl)
              dt$n_phen = nrow(temp) 
              dt$analysis = anl
              message('the table produced has row_col dimentions ', paste(dim(dt), collapse = ' '))
              return(dt)
          }
    }else{
            thisrow
          }
}, mc.cores = parallel::detectCores() - 2), fill = TRUE)



# Export results 
fwrite(results, 'results/drug_pheno_assoc_anal_psychAD.csv')


