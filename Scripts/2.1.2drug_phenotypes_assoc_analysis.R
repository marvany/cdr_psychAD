source("/sc/arion/projects/va-biobank/PROJECTS/ma_cdr_psychAD/Scripts/2.1.3_assoc_an_helper.R")

    # =========================================== #
    #     PART 1: Process PsychAD Phenotypes      #
    # =========================================== #

##############
# DEFINE PATHS
#p.drugs.phenos = 'Resources/result_analysis/test_merged_drugs_phenotypes.csv'
p.psychAD = './Resources/psychAD/clinical_metadata.csv'
p.psychAD.IDmap = '/sc/arion/projects/roussp01a/deepika/merging_psychAD_SNParray_WGS/common_variants_psychAD/ancestry_pca_psychAD_1429_samples/psychAD_20PC_3_methods_ancestry.tsv'


##################
# MUNGE PHENOTYPES

# We will examine which columns that correspond to phenotypes need to be filtered out
# Confirm that all phenotypes are present
dt = fread(p.psychAD)
dt = dt[Ethnicity == 'EUR']
i.dt = ins(dt)
print(i.dt, nrows = Inf)
i.dt[, no_emptyChar := NULL]

print(i.dt, nrows = Inf)
# REMOVE PRS
crm = grep('prs_', i.dt[, variable], value = TRUE)
i.dt = rm.var(i.dt, crm)
print(i.dt, nrows = Inf)

# REMOVE ROWS WITH SEXUAL ANEUPLOIDY
rrm = which(dt[, Sex_chr_aneuploidy])j
dt <- dt[!rrm,]
print(i.dt, nrows = Inf)

# spec_keep are variables that may proveuseful in the future but shouldn't be analayzed
spec_keep = c('SubID', 'Dx')

# REMOVE HIGH-DIMENTIONAL CHARACTER COLS - KEEP Dx, SubID - 
crm = c(crm, i.dt[class == 'character' & n_unique > 10, variable])
crm = setdiff(crm, spec_keep)
i.dt = rm.var(i.dt, crm)
print(i.dt, nrows = Inf)

# REMOVE HIGH NAs
crm = c(crm, i.dt[no_NA > 400, variable])
i.dt = rm.var(i.dt, crm)
print(i.dt, nrows = Inf)

# MANUAL REMOVAL
crm = c(crm, c('Brain_bank', 'Source_Location', 'Death_Season')) # remove character columns that will increase dimentions for no good reason
i.dt = rm.var(i.dt, crm)
print.a(i.dt)

# FIRST FILTER VARIABLES 
dt.fil = dt[, i.dt[, variable], with = FALSE]

dt.fil[, Ethnicity := NULL] # HARD-CODED removal of Ethnicity; BEWARE!


# ENCODE FACTORS
char_cols = names(which(sapply(dt.fil, is.character)))
char_cols # pick characters to factor
to_factor = c('Sex', 'ApoE_gt')
#dt.fil[, (to_factor) := lapply(.SD, as.factor), .SDcols = to_factor]
dt.fil[ , (to_factor) := lapply(.SD, function(x) addNA(as.factor(x))), .SDcols = to_factor] # addNA protects dimentions

#dt.fil = cbind(dt.fil, as.data.table(model.matrix(~ ApoE_gt - 1, data = dt.fil)))
form <- reformulate(to_factor, intercept = FALSE)  # ~ Sex + ApoE_gt -1
dums <- as.data.table(model.matrix(form, data = dt.fil, na.action = na.pass))
dt.fil[, (names(dums)) := dums]

# FIND near-0 variance
print.a(ins(dt.fill))
num_cols <- names(which(sapply(dt.fil, is.numeric)))
i.dt2 <- data.table(
  variable = num_cols,
  variance = vapply(
    dt.fil[, ..num_cols],
    stats::var,                 # fully qualified
    numeric(1L),                # expected length-1 numeric result
    na.rm = TRUE
  )
)

i.dt2[, above_thresh := (variance > 0.01)]
print.a(i.dt2)
i.dt2 = i.dt2[above_thresh == TRUE, .(variable, variance)]
print.a(i.dt2)

# MANUAL REMOVE
i.dt2 <- rm.var(i.dt2, c('ApoE_gtNA', 'Ethnicity'))
print.a(i.dt2)

# SECOND FILTERING
to_keep = c(setdiff(char_cols, to_factor), i.dt2[, variable]) # factorized are now part of the matrix as numeric
dt.fil2 = dt.fil[, ..to_keep]
dim(dt.fil2)
length(to_keep)

fwrite(dt.fil2, 'Resources/psychAD/postProcess_clinical_metadata.csv')


    # =========================== #
    #     PART 2: Run Analysis    #
    # =========================== #             

# ------------------------ Create the table of drugs (dt.drugs) ; contains CDR-results + ID columns 

# We first load the merged.dt (contains combined psychAD and CDR results; created in Jupyter Notebook)
# We will filter out the columns that were added in the ipynb file
library(data.table)
library(pbmcapply)
p.drugs.phenos = 'Resources/result_analysis/cdr_psychAD_IC_microglia_v4.csv'
p.phenos = './Resources/psychAD/clinical_metadata.csv' # this is the original ipynb file used

dt.merged = fread(p.drugs.phenos) # everything will be sourced from here
dt.phenos = fread(p.phenos) # original raw phenotypes file

all.pheno.cols = names(dt.phenos)
all.merged.cols = names(dt.merged)

to_remove = intersect(all.merged.cols, all.pheno.cols) # we want to keep the ID columns
setdiff(all.pheno.cols, to_remove) # only SubID won't be removed, because it doesn't exist, but that's ok; you already have ID cols in merged

dt.drugs = copy(dt.merged)
dt.drugs[, (to_remove) := NULL]
dt.drugs[, 'ind_ID' := NULL]

length(names(dt.merged))
length(names(dt.phenos)) + length(names(dt.drugs)) # makes sense one extra variable for column SubID

# ------------------------ Prepare revised.phenos with matching ID
# Revised phenotypes
dt.revised.phenos = fread('Resources/psychAD/postProcess_clinical_metadata.csv') # Is the data created above

setnames(dt.revised.phenos, 'SubID', 'ind_ID_clin')
f.dt = dt.drugs[dt.revised.phenos, on = .(ind_ID_clin), nomatch = 0]
i.dt = ins(f.dt)
unique(i.dt$n_unique)

# We will iterate over the phenotypes and create tables
revised.pheno.names = setdiff(names(dt.revised.phenos), c('ind_ID_clin', 'Dx'))

to_analyze.list = pbmclapply(
  revised.pheno.names,
  mc.cores = parallel::detectCores() - 2,
  function(x){
    na.omit(f.dt[, c('ind_ID_clin', 'Dx', x), with = FALSE])
})

names(to_analyze.list) = revised.pheno.names
str(to_analyze.list[1])

# Determine what analysis will be run for each variable
model_type_vec <- pbmclapply(
  to_analyze.list,
  determine_analysis,
  mc.cores = parallel::detectCores() - 2
)
model_type_vec <- unlist(model_type_vec, use.names = TRUE)

# ------------------------ Prepare revised.phenos with matching ID
# Run the full analysis for each given phenotype

results_by_pheno <- pbmclapply(
  names(to_analyze.list),
  mc.cores = parallel::detectCores() - 2,
  run_association_analysis
)

final_dt <- rbindlist(results_by_pheno, use.names = TRUE, fill = TRUE)
fwrite(final_dt, "results/drug_vs_phenotype_coeffs.tsv", sep = "\t")

long_dt[sample(.N, 10)]
# Iterate over to_analyze.list with lapply
# merge each data.table of said list with the dt.drugs table  based on matching values in the ind_ID_clin
# extract for given data.table the responding type of analysis from the model_type_vec variable
# execute a separate linear/logistic regression for all of the table columns with the outcome variable of the 'given data.table' (the one that determined the type of anlaysis)
# here is the code for linear/logistic regression you should use, in conjunction with some additional code (as you can see anl is the varialbe that determines what will be executed)
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
# (create_coefficients_df is a custom made function, that you don't need to know anything about, just call it as instructed)



dt.drugs
head(dt.drugs[,1:5])

# TODO:
# Iterate over both lists ; and run linear or logistic regression (revise previously made functions). Ordinal can be skipped for now
# function create_coeffiecients_df() is the one needed to analyze


# CODE TO USE FOR ANALYSIS:
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
    }



###############################################################3
#   OLD

insp_psychAD <- ins(psychAD)
print(insp_psychAD, nrows = Inf)
keep_those = c()


insp_psychAD[class == 'character']
psychAD[, unique(Ethnicity)]
to_remove = insp_psychAD[n_unique < 30, variable]

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
ins(df.merged[, setdiff(names(df.merged), names(dt)), with = FALSE])














################### OLD VERSION BELOW
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


