# Environment
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

ins <- function(dt, spec_value = NULL){
  to_return = data.table(
    variable = names(dt),
    class = sapply(dt, class),
    n_unique = sapply(dt, uniqueN),
    no_NA = sapply(dt, function(x) sum(is.na(x))),
    no_emptyChar = sapply(dt, function(x) sum(x %in% ''))
  )
  if(!is.null(spec_value)) cbind(to_return, data.table(no_spec_value = sapply(dt, function(x) sum(x %in% spec_value))))
  to_return
}
rm.var <- function(thisdt, these){
  thisdt[!(variable %in% these),]
}


print.a <- function(thisdt) print(thisdt, nrows = Inf)


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
              quote = FALSE, sep = '\t', row.names = FALSE)
              }

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

# Determines whether linear vs. logistic regression will be run for a given variable
determine_analysis <- function(dt) {
  if (!is.data.table(dt))
    stop("Input must be a data.table")

  # Count *distinct* non-missing values in the 3rd column
  k <- uniqueN(dt[[3]], na.rm = TRUE)

  if (k == 2L) "logistic" else "linear"
}



    # ===================================== #
    #       ASSOCIATION ANALYSIS PART       #
    # ===================================== #

## ------------------------------------------------------------------
##  Function that runs the full analysis for one phenotype -------
## ------------------------------------------------------------------
run_association_analysis <- function(pheno_name) {

  ph_dt <- to_analyze.list[[pheno_name]]          # id + Dx + phenotype
  anl   <- model_type_vec[[pheno_name]]           # “linear” or “logistic”

  ## ---- Merge with drug signatures on ind_ID_clin -----------------
  tmp   <- merge(
    ph_dt,
    dt.drugs,
    by       = "ind_ID_clin",
    all      = FALSE          # inner join → keep only matching IDs
  )

    ## ---- Identify drug columns -------------------------------------
    cmp_cols <- setdiff(names(tmp), c("ind_ID_clin", "Dx", pheno_name))
    
    ## ---- Run one regression per compound ---------------------------
    res_list <- pbmclapply(
      cmp_cols,
      mc.cores = parallel::detectCores() - 2,
      function(cmp) {
        tryCatch({
          
        temp <- tmp[ , .(
                  phen = get(pheno_name),
                  cmp  = get(cmp)
                )]

        temp <- temp[complete.cases(temp)]
            
        temp[, names(temp) := lapply(.SD, scale)]
        
        if (nrow(temp) < 3L) return(NULL)            # skip tiny samples

        if (anl == "linear") {
          fit.lm <- lm (phen ~ cmp, data = temp)
          fit    <- glm(phen ~ cmp, data = temp, family = gaussian(link = "identity"))

        } else {                                     # logistic
          temp[ , phen := factor(phen) ]
          fit    <- fit.lm <- glm(phen ~ cmp, data = temp, family = binomial)
        }
        out <- as.data.table(create_coeffiecients_df(
                fit  = fit,
                fit.lm = fit.lm,
                phen = pheno_name,
                cmp  = cmp,
                anl  = anl
              ))
        out[ , `:=`( n_phen  = nrow(temp),
                    analysis = anl ) ]
        out
       }, error = function(e){
          message(e)
      })
      })
      
      rbindlist(res_list, use.names = TRUE, fill = TRUE)
}






run_association_analysis <- function(pheno_name) {

  ## ---- 1.  Prepare data ---------------------------------------------
  ph_dt <- to_analyze.list[[pheno_name]]           # id + Dx + phenotype
  anl   <- model_type_vec[[pheno_name]]            # "linear" | "logistic"

  tmp <- merge(
    ph_dt,
    dt.drugs,
    by  = "ind_ID_clin",
    all = FALSE                                     # inner join
  )

  cmp_cols <- setdiff(names(tmp), c("ind_ID_clin", "Dx", pheno_name))

  ## ---- 2.  Melt to long format once ---------------------------------
  long_dt <- melt(
    tmp,
    id.vars       = c("ind_ID_clin", "Dx", pheno_name),
    measure.vars  = cmp_cols,
    variable.name = "cmp_name",       # identifier
    value.name    = "cmp"
  )


  ## keep rows where both phenotype & cmp are present
  long_dt <- long_dt[!is.na(get(pheno_name)) & !is.na(cmp)]

  ## ---- 3.  One regression per compound, but in a single DT pass -----
  res <- long_dt[ , {

      ## scale within-compound (replicates original behaviour)
      this <- copy(.SD)
      this[ , `:=`(
        phen = as.numeric(scale(get(pheno_name))),
        cmp  = as.numeric(scale(cmp))
      )]

      if (nrow(this) < 3L) return(NULL)            # skip tiny samples

      if (anl == "linear") {
        fit.lm <- lm (phen ~ cmp, data = this)
        fit    <- glm(phen ~ cmp, data = this, family = gaussian("identity"))
      } else {                                     # logistic
        this[ , phen := factor(phen) ]
        fit <- fit.lm <- glm(phen ~ cmp, data = this, family = binomial)
      }

      out <- as.data.table(create_coeffiecients_df(
        fit     = fit,
        fit.lm  = fit.lm,
        phen    = pheno_name,
        cmp     = .BY$cmp_name,
        anl     = anl
      ))

      out[ , `:=`(
        n_phen  = nrow(this),
        analysis = anl
      )]
    },
    by = cmp_name
  ]

  res[]                                            # return a tidy DT
}
