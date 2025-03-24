######################
#### DEPENDENCIES ####

#working.dir = '/sc/arion/projects/va-biobank/PROJECTS/ma_neurocog_grex'
#setwd(working.dir)
#library(data.table)
#source('./Scripts/stat_calcVarPart.R')
#source('./Scripts/stat_PseudoR2.R')
library(ordinal)
library(MASS)
library(parallel)

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

####################################################
# OUTLINE
# PREPARE FOR JOB HANDLER (preparation before job handler in run_GReX_to_pheno_linear.R)
# JOB HANDLER (executed in run_GReX_to_pheno_linear.R)
# HELPER SCRIPT (called by )
# PRS_to_pheno FUNCTION
# LINEAR_PSEUDOR2 FUNCTION


################################################################################################
### DEVELOP GREX_TO_PHENO_JOB HANDLER (PENDING ITEMS; CREATE MASTER SCRIPT, CREATE BSUB COMMAND)
GReX_to_pheno_job_handler_v3 = function(
    path_to_spec_grex, # list that contains sublist 'df' where PRS, PCs and phenotypes are included
    phenotypes = pheno, # phenotypes to iterate over
    outDir = outDir, # output directory, when run for GReX over multiple tissues this should be unique for each tissue
    # additional arguments that arae passed into logit_pseudoR, see below:
    exclude.age = F, # exlude age?
    only.PRS = F, # "only.PRS" is deprecated as term, this will exclude other covariates 
    extra = NA, # extra covariates
    custom.covar = NA, # completely custom covariates e.g. c("age", "sex"),
    helper_script,
    working.dir,
    type.of.analysis,
    walltime,
    memory,
    n.cores,
    path_to_sex,
    path_to_age,
    path_to_pcs){

  # In case of large arrays this functions helps us in splitting them in chunks of 10000
  compute_values <- function(input_integer) {
    a <- input_integer %/% 10000 # divides the left operand by the right operand and returns only the integer part, discarding any fraction
    b <- input_integer %% 10000 # similary to above but return the fraction.
    list(a = a, b = b)
  }
  
  ############################################
  # PREPARE VARIABLES / PARSE NAMES
  #where all the job logs will go
  jobsdir <- paste0(outDir, "/small_script_jobs")
  jobsdir <- gsub('^\\.', getwd(), jobsdir) # we like full paths
  
  if (!dir.exists(jobsdir)){dir.create(jobsdir, recursive = TRUE)}
  
  # custom.covar may be of format c('age', 'sex')
  custom.covar = unlist(strsplit(custom.covar, split = ','))
  
  
  
  ##########################################
  # CREATE SMALL SCRIPTS
  # Creates all the small .sh files that each calls the helper script and passes it all the arguments (including ENSG12345_phenotype_model.txt)
  create_commandV2(
    these.gfs = these.gfs,
    jobsdir = jobsdir, helper_script=helper_script,
    all.rds.file = all.rds.file, pheno = pheno, outDir = outDir, exclude.age = exclude.age,
    only.PRS = only.PRS, extra = extra, custom.covar = custom.covar, working.dir = working.dir,
    type.of.analysis = type.of.analysis, n.cores = n.cores, path_to_sex = path_to_sex, path_to_age = path_to_age,
    path_to_pcs = path_to_pcs)
  
  ##########################################
  # CREATE MASTER SCRIPT
  # 'run.scripts.sh' changes directory to the 'path_to_core.scripts' (see above line) and submits the containing scripts according to the $LSB_JOBINDEX
  outDir.full.name = gsub("\\.", getwd(), outDir)
  master.array.script.sh <- paste0('#!/bin/bash\n',
                                   'cd ', paste0(outDir.full.name, '/small_script_files'), '\n', # inside there you will fines script_1.sh , script_2.sh etc..
                                   'scripts=($(ls))\n',
                                   'thisscript=${scripts[$LSB_JOBINDEX-1]}\n',
                                   'sh $thisscript')
  
  
  
  path.to.master.array.script = paste0(outDir.full.name, "/master_array_script.sh")
  # export the command
  writeLines(master.array.script.sh, path.to.master.array.script)
  
  # make the command executable
  system(paste0('chmod +x ', path.to.master.array.script))
  
  ##########################################
  # CREATE BSUB COMMAND
  # create bsub  # CURRENT 012025 03:29 AM IS PRELIMINARY
  # convert to './path/to/whatever' to '/path/to/whatever' (remove the dot)
  path.to.master.array.script = gsub('^\\.', getwd(), path.to.master.array.script)
  
  
  small.scripts.dir = paste0(outDir.full.name, "/small_script_files")
  no.small.scripts = length(list.files(small.scripts.dir)) # will define the size of the array
  
  integer.fraction = compute_values(no.small.scripts)
  
  total.iterations = integer.fraction$a + 1
  cells.remaining = no.small.scripts
  
  for(i in 1:total.iterations){
        
        start.ind =  1 + (i-1)*10000
        end.ind = min(i*10000, no.small.scripts)
          
        b.sub <- paste0('bsub -J ', paste0(basename(outDir), "[", start.ind,"-", end.ind, "]"),
                      ' -P acc_va-biobank',
                      ' -q premium',
                      ' -n ', n.cores,
                      ' -W ', walltime,
                      ' -R "span[hosts=1]"',
                      ' -R "affinity[core(1)]"',
                      ' -R \"rusage[mem=', round(memory), ']\"',
                      ' -oo ', paste0(jobsdir, '/script_%I.out'),
                      ' -eo ', paste0(jobsdir, '/script_%I.err'),
                      ' -L /bin/bash ')

        to_execute = paste0(b.sub, path.to.master.array.script)
        
        system(paste0("chmod +x ", path.to.master.array.script))
        system(to_execute)
  }
}





####################################################
# JOB HANDLER
### CALL JOB HANDLER (COPY ME PLEASE)
GReX_to_pheno_linear_job_handler_v2 = function(
  path_to_spec_grex, # list that contains sublist 'df' where PRS, PCs and phenotypes are included
  phenotypes = pheno, # phenotypes to iterate over
  outDir = outDir, # output directory, when run for GReX over multiple tissues this should be unique for each tissue
  # additional arguments that arae passed into logit_pseudoR, see below:
  exclude.age = F, # exlude age?
  only.PRS = F, # "only.PRS" is deprecated as term, this will exclude other covariates 
  extra = NA, # extra covariates
  custom.covar = NA, # completely custom covariates e.g. c("age", "sex"),
  helper_script,
  working.dir,
  type.of.analysis,
  walltime,
  memory){

  # where all the job logs will go
  jobsdir <- paste0(outDir, "/JOBS")
  jobsdir <- gsub('^\\.', getwd(), jobsdir)

  if (!dir.exists(jobsdir)){dir.create(jobsdir, recursive = TRUE)}

  # custom.covar may be of format c('age', 'sex')
  custom.covar = unlist(strsplit(custom.covar, split = ','))
  
  # I will split the jobs for each phenotype (one core per phenotype)
  #for (pheno in phenotypes){
    expected.file <- paste0(outDir, "/", pheno, ".prediction.tsv")
    if(TRUE){ #if (!file.exists(expected.file)){

      # make this the same as the antagonist # create_command() writes a command to be executed by bash
      # to "all.rds.file" to xrhsimopoioume kataxristika
     path_to_command = create_command(
              jobsdir = jobsdir, 
              helper_script=helper_script,
              all.rds.file = path_to_spec_grex, 
              pheno = pheno, 
              outDir = outDir, 
              exclude.age = exclude.age,
              only.PRS = only.PRS, extra = extra, custom.covar = custom.covar,
              working.dir = working.dir, type.of.analysis = type.of.analysis)

      # convert to './path/to/whatever' to '/path/to/whatever' (remove the dot)
      path_to_command = gsub('\\.', getwd(), path_to_command)
      path_to_command = paste0(path_to_command, '.sh')
      bsub = paste0('cd ', jobsdir, ' && bsub ', '-J "', pheno, '" -q premium -P acc_va-biobank -n 1 -W ',walltime,' -R "span[hosts=1]" -R "rusage[mem=', memory,']"',
              ' -oo ', pheno, '.out -eo ', pheno, '.err',' -L /bin/bash ')
      to_execute = paste0(bsub, path_to_command)
      system(paste0("chmod +x ", path_to_command))
      system(to_execute)
    } # generate the results file
  #} # all phenotypes processed
}


####################################################
####################################################
####################################################
####################################################
# JOB HANDLER
### CALL JOB HANDLER (COPY ME PLEASE)
# DEPRECATED!!!
GReX_to_pheno_linear_job_handler = function(
  all = prepared_all, # list that contains sublist 'df' where PRS, PCs and phenotypes are included
  phenotypes = pheno, # phenotypes to iterate over
  outDir = outDir, # output directory, when run for GReX over multiple tissues this should be unique for each tissue
  # additional arguments that arae passed into logit_pseudoR, see below:
  exclude.age = F, # exlude age?
  only.PRS = F, # "only.PRS" is deprecated as term, this will exclude other covariates 
  extra = NA, # extra covariates
  custom.covar = NA, # completely custom covariates e.g. c("age", "sex"),
  helper_script,
  working.dir,
  type.of.analysis){
  # If it is a character then read the file (preferred method of submission from the PRS_to_pheno_job_handler function
  if (is.character(all)) {all <- readRDS(all)}
  
  # now it has to be saved independently so it can be loaded by the jobs
  all.rds.file <- paste0(outDir, "/all.RDS") ### added phenotypes to all.rds.file file name, for when running multiple phenotypes separately 

  # creates outDir
  if (!dir.exists(dirname(all.rds.file))) dir.create(dirname(all.rds.file), recursive = TRUE)

  # saves all file
  saveRDS(all, all.rds.file)

  # where all the job logs will go
  jobsdir <- paste0(outDir, "/JOBS")
  jobsdir <- gsub('^\\.', getwd(), jobsdir)

  if (!dir.exists(jobsdir)){dir.create(jobsdir, recursive = TRUE)}

  # custom.covar may be of format c('age', 'sex')
  custom.covar = unlist(strsplit(custom.covar, split = ','))
  
  # I will split the jobs for each phenotype (one core per phenotype)
  #for (pheno in phenotypes){
    expected.file <- paste0(outDir, "/", pheno, ".prediction.tsv")
    if(TRUE){ #if (!file.exists(expected.file)){

      # make this the same as the antagonist # create_command() writes a command to be executed by bash
     path_to_command = create_command(jobsdir = jobsdir, helper_script=helper_script,
        all.rds.file = all.rds.file, pheno = pheno, outDir = outDir, exclude.age = exclude.age,
        only.PRS = only.PRS, extra = extra, custom.covar = custom.covar, working.dir = working.dir, type.of.analysis = type.of.analysis)

      # convert to './path/to/whatever' to '/path/to/whatever' (remove the dot)
      path_to_command = gsub('\\.', getwd(), path_to_command)
      path_to_command = paste0(path_to_command, '.sh')
      bsub = paste0('cd ', jobsdir, ' && bsub ', '-J "', pheno, '" -q premium -P acc_va-biobank -n 1 -W 01:00 -R "span[hosts=1]" -R "rusage[mem=5000]"',
              ' -oo ', pheno, '.out -eo ', pheno, '.err',' -L /bin/bash ')
      to_execute = paste0(bsub, path_to_command)
      system(paste0("chmod +x ", path_to_command))
      system(to_execute)
    } # generate the results file
  #} # all phenotypes processed
}
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
###################################
# CREATES COMMAND

create_commandV2 = function(
    these.gfs = these.gfs,
    jobsdir = jobsdir, helper_script=helper_script,
    all.rds.file = all.rds.file, pheno = pheno, outDir = outDir, exclude.age = exclude.age,
    only.PRS = only.PRS, extra = extra, custom.covar = custom.covar, working.dir = working.dir,
    type.of.analysis = type.of.analysis, ind = ind, n.cores = n.cores,
    path_to_sex = path_to_sex, path_to_age = path_to_age,
    path_to_pcs = path_to_pcs){
  
  pbmcapply::pbmclapply(seq_along(these.gfs), function(ind){
  #for(ind in seq_along(these.gfs)){
    path_to_script.sh = paste0(outDir, '/small_script_files/script_', ind,'.sh')
    if(!dir.exists(dirname(path_to_script.sh))) dir.create(dirname(path_to_script.sh), recursive = TRUE)
    
    if(!file.exists(path_to_script.sh)){
              p.to.genes.file = these.gfs[ind]
              
              ################
              # Create command
              thiscommand = paste0(
                "cd ", working.dir, "\n",
                "ml R", "\n",
                "Rscript --verbose ", helper_script, " ",
                #### GENISIS/Georgios local config
                ### 'Rscript /home/home1/vhabrxvoloug/Desktop/Protocols/PRScs/1.BED_pipeline/PRS_to_pheno_job_handler_helper_script.R ',
                p.to.genes.file, ' ', # path to csv file ; ENSG...phenotype...model.txt.gz
                pheno, ' ', # phenotype to be tested
                outDir, ' ', # output directory to save results
                exclude.age, ' ', # whether to exclude age from the analysis
                only.PRS, ' ', # whether to only perform PRS analysis without other covariates
                extra, ' ', # additional covariates to be included in the model
                paste0(custom.covar, collapse=","), ' ', # force custom covariates
                working.dir, ' ',
                type.of.analysis, ' ',
                path_to_sex, ' ',
                path_to_age, ' ',
                path_to_pcs
              )   
              ################
              # We write the command to /Results/biobank/pheno/model/ENSG00000123456_phenotype_w_pheno_GReX_model/executed_command
              
              # 'thiscommand' is exported to its respective .sh file (to be executed by run.scripts.sh)
              writeLines(
                paste0(
                  "#!/bin/bash", "\n",
                  thiscommand),
                path_to_script.sh
              )
    }

    system(paste0("chmod +x ", path_to_script.sh))
  }, mc.cores = parallel::detectCores() - 2)
}




# deprecated!!!!
create_command = function(jobsdir = jobsdir, helper_script=helper_script,
                    all.rds.file = all.rds.file, pheno = pheno, outDir = outDir, exclude.age = exclude.age,
                    only.PRS = only.PRS, extra = extra, custom.covar = custom.covar, working.dir = working.dir, type.of.analysis = type.of.analysis){
                command_dir = paste0(outDir, '/executed_command')
                        
                if(!dir.exists(command_dir)) dir.create(command_dir, recursive = TRUE)
                # This is the command that executes the core.function
                thiscommand <- paste0(
                    "cd ", jobsdir, "\n",
                    "ml R", "\n",
                    "Rscript --verbose ", helper_script, " ",
                    #### GENISIS/Georgios local config
                    ### 'Rscript /home/home1/vhabrxvoloug/Desktop/Protocols/PRScs/1.BED_pipeline/PRS_to_pheno_job_handler_helper_script.R ',
                    all.rds.file, ' ', # path to the RDS file containing the data to be processed by the jobs
                    pheno, ' ', # phenotype to be tested
                    outDir, ' ', # output directory to save results
                    exclude.age, ' ', # whether to exclude age from the analysis
                    only.PRS, ' ', # whether to only perform PRS analysis without other covariates
                    extra, ' ', # additional covariates to be included in the model
                    paste0(custom.covar, collapse=","), ' ', # force custom covariates
                    working.dir, ' ',
                    type.of.analysis
                    )
                    
                    # 'thiscommand' is exported to its respective .sh file (to be executed by run.scripts.sh)
                    writeLines(
                    paste0(
                        "#!/bin/bash", "\n",
                        thiscommand),
                    paste0(file.path(command_dir, paste0(pheno, '.sh')))
                    )
                    path_to_command = paste0(file.path(command_dir, pheno)) # we add '.sh' outside the function to avoid errors.
                    return(path_to_command) # return so that you can submit
}


########################################################################################################
#  SPLITS PER PHENOTYPE

# The whole PRScs wrapper
GReX_to_pheno_linear <- function(
  all = NA, # df where PRS, PCs and phenotypes are included
  phenotypes = NA, # if pheno is not NA then will do all phecodes (not recommended) or whatever is in all$phe
  outDir, # output directory
  type.of.analysis,
  ...){
    # additional arguments that are passed into logit_pseudoR, see below:
    #  exclude.age = F, # exlude age?
    #  only.PRS = F, # only do PRS do not see other covariates
    #  extra = NA, # extra covariates
    #  custom.covar = NA # completely custom covariates e.g. c("age", "sex")

    # If it is NA then load the phecodes (will only work as standalone)
    ### if (is.na(all)) {all <- readRDS("")} #### can add path to specific file here (won't be using this for Social Genomics project)
    # If it is a character then read the file (preferred method of submission from the PRS_to_pheno_job_handler function
    if (is.character(all)) {all <- readRDS(all)}
    # Do your magic
    #decide <- all$PRSinfo
    #if (is.na(traits)) {traits <- all[["PRS"]]} else {
    #  decide <- decide[trait %in% traits]
    #  traits <- decide$PRS }
    #if (!is.na(specify.phi)) {traits <- decide[phi == specify.phi]$PRS}
    if (is.na(phenotypes)) {phenotypes <- all$phe}
    message('There are ', length(phenotypes), ' phenotypes in total.')
    message('These are the phenotypes ',paste(phenotypes, sep = ' '))
    message('Type of analysis is ', type.of.analysis)
  
    # this list is to meet conventions (it replaces the original all list, that carries all[['df']])
    conventional_list = list(df = all)
    message('class of all is ', class(all))
    rm(all)
    message(paste0('class of conventional_list is: ',class(conventional_list[['df']])))
    # lapply(phenotypes, 
      #       FUN = function(x) { # x is the phenotype here
      x = phenotypes
  
  if(TRUE){#if (!file.exists(paste0(outDir, "/", x, ".prediction.tsv"))){
            #outputdir <- paste0(outDir, "/", x) # FOR THE ALL-GENES ANALYSIS THIS SHOULD BE RUN BUT FOR THE SINGLE GENE NO NEED TO SPLIT, EVERYTHING IS IN NAME
            outputdir = outDir
            message('x is ', x)
            message('printout is sent to outputdir which is: ', outputdir)
            if(!(type.of.analysis %in% c('linear', 'logistic'))) stop('For phenotype ', x, ' you have provided an ineligible type of analysis.')
            if(type.of.analysis == "linear") results <- GReX_linear_rsquared(ptt = 'GReX', x, conventional_list, outputdir, ...)
            if(type.of.analysis == "logistic") results <- GReX_logistic_pseudoR(ptt = 'GReX', x, conventional_list, outputdir, ...)
            
            #df <- linear_pseudoR_df_conv(results, phenotype = x)
            message('results = coefficients_df is exported to: ', outDir, ' with name ', x)
            #WriteProperTSV(results, paste0(outDir, "/", x, ".prediction.tsv")) 
            message('gg wp')
    } # only execute this if output does not exist
        #} # lapply function is done
  #) # lapply is done 
}


################

GReX_linear_rsquared <- function(
  ptt, # PRS to test e.g. SCZ.3 # e.g. ptt = "R3_imp_maf0.01_geno0.02_hwe.minus.10_rsq0.9...SCZ.2...phiauto"
  ttt, # trait to test phe_295_1 # e.g. for creativity project ttt = "wrktypart"
  all, # all object as generated 
  outdir, # outdir to save file
  exclude.age = F, # exlude age?
  only.PRS = F, # only do PRS do not see other covariates
  extra = NA,
  custom.covar = NA, # completely custom covariates e.g. c("age", "sex")
  gene.name = gene.name,
  path_to_sex = path_to_sex,
  path_to_age = path_to_age,
  path_to_pcs = path_to_pcs
  ){  # extra covariates
  # column names in data.table do not necessarily play well with glm
  # makesafe <- function(x) {return(gsub("-", "_", x))}
  # since the limitation is that bad names go as incompatible arguments in the lm
  # we will use my make_names(). I have also applied that in the PRS script before
  outdir = paste0(outdir, "/", gene.name)
  if (only.PRS) {outdir <- paste0(outdir,".onlyPRS")}
  if (!dir.exists(outdir)) dir.create(outdir)
  names(all[["df"]]) <- make_names(names(all[["df"]]))
  #all[["PRS"]] <- make_names(all[["PRS"]])
  ptt <- make_names(ptt)
  ttt <- make_names(ttt) # e.g. for creativity project ttt = "wrktypart"
  
  # It is important that the PRS is first to extract the measurements
  if(ptt == 'GReX'){
    ptt = colnames(all[["df"]])[grepl("ENSG", colnames(all[["df"]]))]
    }
  
  if (!any(is.na(custom.covar))) {
    frm <- paste(unlist(c(ptt, custom.covar)), collapse = "+")
  } else {
    if (is.na(extra)) {
      if (only.PRS) {
        frm <- paste(unlist(c(ptt)), 
                     collapse = "+")
      } else {
        frm <- paste(unlist(c(ptt, all[["PCs"]],all[["sex-age"]])), # MA: DO PCA
                     collapse = "+")
      }
    } else {
      if (only.PRS) {
        frm <- paste(unlist(c(ptt, extra)), 
                     collapse = "+")
      } else {
        frm <- paste(unlist(c(ptt, extra, all[["PCs"]],all[["sex-age"]])), # MA: MAKE COMPATIBLE FORMAT
                     collapse = "+")
      }
    }  
  } # end loop for if not custom covariates are given
  # remove age if specified
  if (exclude.age) {frm<-sub("\\+age\\+age\\.sq","",frm)}
  # glm does not like quotes 
  form <- reformulate(frm, response = ttt)
  #form  <- gsub("\"","\`",form)
  # rescale the PRS
  

  #######################
  # Add age sex PCs age^2
  # DEBUG
  # path_to_sex = './Resources/demographics/UKBB/ukbb_data_sex_values_updated.tsv'
  # path_to_sex = "./Resources/demographics/ABCD/reported_sex_ABCD.tsv"
  # path_to_age = './Resources/demographics/UKBB/ukbb_data_f21003_age_2020.tsv'
  # path_to_age = "./Resources/demographics/ABCD/age_ABCD.tsv"
  # path_to_pcs = "./Resources/demographics/UKBB/PC1_20.csv"
  # path_to_pcs = "./Resources/demographics/ABCD/PC1_20.csv"
  # df = fread("./Resources/GReX_to_pheno/ABCD/grex_pheno_comb/filtered_grex_pheno_comb/pea_ravlt_sd_trial_vi_tc/w_pheno_GReX_230721_MegaAnalysis_EUR_subclass_SMC_filtered/ENSG00000108654_pea_ravlt_sd_trial_vi_tc_w_pheno_GReX_230721_MegaAnalysis_EUR_subclass_SMC_filtered.txt.gz")
  # df = fread("./Resources/GReX_to_pheno/UKBB/grex_pheno_comb/filtered_grex_pheno_comb/trail_making_alphanumeric/w_pheno_GReX_230721_MegaAnalysis_EUR_subclass_SMC_filtered/ENSG00000006715_trail_making_alphanumeric_w_pheno_GReX_230721_MegaAnalysis_EUR_subclass_SMC_filtered.txt.gz")
  # all = list()  


  df = as.data.table(all[['df']])
  
  # Scale phenotype if unscaled
  this_phenotype = names(df)[3]
  message('Now I will check if phenotype', this_phenotype, ' is scaled.')
  if(sd(unlist(df[,3])) != 1){
      message('Phenotype ', names(df)[3], ' was not scaled, I will proceed with scaling (standard deviation is not 1)')
      df[[this_phenotype]] = as.numeric(scale(as.numeric(df[[this_phenotype]])))
      thissd = sd(unlist(df[,3]))
      message('After scaling phenotype ', this_phenotype, ' the standard deviation is ', sd(unlist(df[,3])))
  }  


  sex = fread(path_to_sex)
  age = fread(path_to_age)
  pcs = fread(path_to_pcs)

  df = merge(df, sex, by = "sample_id", all.x = TRUE)
  df = merge(df, age, by = "sample_id", all.x = TRUE)
  df = merge(df, pcs, by = "sample_id", all.x = TRUE)
  df$age_sq = df$age^2
  df$sex = ifelse(df$sex == "Male", 1, 0)
  df = na.omit(df)
  all[["df"]] = df

  message("Formula to be used: ", form)
  message("Names of columns are ", names(all[["df"]]))


  thesecols = colnames(all[["df"]])
  gene.col = names(all[['df']])[grepl('ENSG', names(all[['df']]))]
  if(all(as.numeric(df[[gene.col]]) == 0)){
    message('For this gene all entries are 0')
    fwrite(data.frame(gene.col = df[[gene.col]]), paste0(outdir, '/empty_gene_', ttt, "_safe.tsv"))
    return(NULL)
  }
  

  message('gene name is ', gene.col)
  for(thispredictor in thesecols){
    if(thispredictor %in% c(custom.covar, gene.col)){
      if(thispredictor != "sex"){
        message("Now scaling ", thispredictor)
        all[["df"]][[thispredictor]] <- as.numeric(scale(as.numeric(all[["df"]][[thispredictor]]))) 
        message("Scaled ", thispredictor)
        message("In total ", thispredictor, " has ", sum(is.na(all[["df"]][[thispredictor]])), " NA values.")
      }
    } else message("Column ", thispredictor, " will not be scaled.")
  }
  
  # Get a report of NA values
  na_report = list()
  for(thiscol in thesecols[!grepl("ENSG", thesecols)]){
   na_report[[thiscol]] =  sum(is.na(all[["df"]][[thiscol]]))
  }
  na_report = as.data.frame(na_report)
  fwrite(na_report, paste0(outdir, "/na_report_", ttt, ".prediction.tsv"))
  
  # get number of observations against number of predictors
  gene_cols = colnames(all[["df"]])[grepl("ENSG", colnames(all[["df"]]))]
  dimension_report = data.frame(
    n_genes = length(gene_cols),
    n_observations = nrow(all[["df"]]),
    n_non_NA_rows = nrow(na.omit(all[["df"]]))
  )
  fwrite(dimension_report, paste0(outdir, '/dimentions_report_', ttt, "_safe.tsv"))
  
  

  fit <- glm(formula = form, data = all[["df"]], family=gaussian(link = "identity")) 
  ### ^ main difference from logistic version above: https://www.statmethods.net/advstats/glm.html, also see help(family) and https://cran.r-project.org/web/packages/GlmSimulatoR/vignettes/exploring_links_for_the_gaussian_distribution.html
  fit.lm <- lm(formula = form, data = all[["df"]])
  
  # Extract the data that I need
  test <- list(
    #"model" = fit[names(fit) %in% c("coefficients", "effects", "R", "formula")], # file is too big otherwise (~450 MBs for EUR each)
    "model.summary" = summary(fit),
    "model.summary.lm" = summary(fit.lm),
    ### "PseudoR" = PseudoR2(fit), #### Not needed, because linear regression 
    "calcVarPart" = calcVarPart(fit)
  ) 
  
  file.target  <- paste0(outdir, "/model_summary_GReX_to_", ttt, ".txt")
  #model.target <- paste0(outdir, "/results/", ptt, "_to_", ttt, "_fitted_model.RDS")
  if (!dir.exists(dirname(file.target))) {
    dir.create(dirname(file.target), recursive=T)}
  #saveRDS(test, model.target)
  zz <- file(file.target, "w")
  sink(zz)
  print(form)
  print(test[["model.summary"]]) #### fixed from "model" only, pulled out added summary()
  ### print(summary(test[["model"]])) 
  ### print(test[["PseudoR"]]) #### Not needed, because linear regression
  sink()
  close(zz)



  ## estimate (Beta), z and p value from glm
  # MA: THIS IS WRONG FOR GREX, BECAUSE GREX IS NOT A SINGLE VALUE
  #x <- as.data.frame(coef(test$model.summary)[2,]) # 2 is for PRS
  #colnames(x) <- ptt

  ####################
  # MA: GREX ADDITIONS
  # Extract coefficients
  coefficients <- coef(summary(fit.lm))[-1,] # -1 removes intercept
  coefficients_df <- as.data.frame(coefficients, nrow = length(coefficients), byrow = TRUE)
  coefficients_df$phenotype = ttt
  coefficients_df$predictor = rownames(coefficients)
  coefficients_df$model_ID= basename(outdir)
  coefficients_df = coefficients_df[, c("phenotype", "model_ID", "predictor", "Estimate", "Std. Error", "t value", "Pr(>|t|)")]

  # Define linearly dependent predictors (so that you can exclude them from your analysis)
  linear_dependent_predictors = setdiff(names(fit.lm$coefficients)[-1], coefficients_df$predictor)

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
  coefficients_df$fdr.adjusted = p.adjust(coefficients_df$`Pr(>|t|)`, method = "fdr")
  coefficients_df$bonferroni.adjusted = p.adjust(coefficients_df$`Pr(>|t|)`, method = "bonferroni")
  coefficients_df$BH.adjusted = p.adjust(coefficients_df$`Pr(>|t|)`, method = "BH")

  # Extract multiple R-squared
  coefficients_df$multiple_r_squared <- summary(fit.lm)$r.squared

  # Extract adjusted R-squared
  coefficients_df$adjusted_r_squared <- summary(fit.lm)$adj.r.squared

  # Extract Residual Deviance
  rse = summary(fit.lm)$sigma
  freedom_degrees = df.residual(fit.lm) # degrees of freedom
  coefficients_df$deviance = rse^2 * freedom_degrees
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
  var_part = calcVarPart(fit)

  # Get residual variance (estimate of the noise, and it is the sum of the partioned variance)
  coefficients_df$residual_variance <- var_part[grepl('Residuals', names(var_part))]

  # Get partitioned variance for each predictor
  to_exclude = c('Residuals', linear_dependent_predictors)
  for(i in 1:length(to_exclude)){
  var_part = var_part[!grepl(to_exclude[i], names(var_part))]
  #message('i = ', i , ' to_exclude is ', to_exclude[i],' var_part has length ', length(var_part))
  }
  coefficients_df$varPart <- var_part
  
  fwrite(coefficients_df, paste0(outdir, "/linear_", gene.name, ".prediction.tsv"))

  return(coefficients_df)
}

GReX_logistic_pseudoR <- function(
  ptt, # PRS to test e.g. SCZ.3 # e.g. ptt = "R3_imp_maf0.01_geno0.02_hwe.minus.10_rsq0.9...SCZ.2...phiauto"
  ttt, # trait to test phe_295_1 # e.g. for creativity project ttt = "wrktypart"
  all, # all object as generated 
  outdir, # outdir to save file
  exclude.age = F, # exclude age?
  only.PRS = F, # only do PRS do not see other covariates
  extra = NA, # extra covariates
  custom.covar = NA # completely custom covariates e.g. c("age", "sex")
  ){
  outdir = paste0(outdir, "/", gene.name)
  # Ensure column names in data.table are compatible with glm
  if (only.PRS) {outdir <- paste0(outdir,".onlyPRS")}
  if (!dir.exists(outdir)) dir.create(outdir)
  names(all[["df"]]) <- make_names(names(all[["df"]]))
  all[["PRS"]] <- make_names(all[["PRS"]])
  ptt <- make_names(ptt) 
  ttt <- make_names(ttt) 
  
  # It is important that the PRS is first to extract the measurements
  if(ptt == 'GReX'){
    ptt = colnames(all[["df"]])[grepl("ENSG", colnames(all[["df"]]))]
    }
  
  custom.covar = unlist(strsplit(custom.covar, split = ','))
  
  # Construct the formula for the model
  if (!any(is.na(custom.covar))) {
    frm <- paste(unlist(c(ptt, custom.covar)), collapse = "+")
  } else {
    if (is.na(extra)) {
      if (only.PRS) {
        frm <- paste(unlist(c(ptt)), collapse = "+")
      } else {
        frm <- paste(unlist(c(ptt, all[["PCs"]],all[["sex-age"]])), collapse = "+")
      }
    } else {
      if (only.PRS) {
        frm <- paste(unlist(c(ptt, extra)), collapse = "+")
      } else {
        frm <- paste(unlist(c(ptt, extra, all[["PCs"]],all[["sex-age"]])), collapse = "+")
      }
    }  
  }
  
  # Remove age if specified # MA: ASK GV IF RESCALING IS NEEDED
  #if (exclude.age) {frm <- sub("\\+age\\+age\\.sq","",frm)}
  
  # Dissect the values of all[['df']][[ptt]] into 3 quantiles
  df = as.data.table(all[['df']])
  
  nrow_before = nrow(df)

  df = df[!is.na(get(ttt))]
  
  nrow_after = nrow(df)

  all[['df']] = df
  
  temp_df = data.frame(
    phenotype = ttt,
    before_filtering = nrow_before,
    after_filtering = nrow_after
  )
  
  fwrite(temp_df, paste0(outdir, "/number_of_rows_", ttt, ".tsv"))

  # Convert the outcome to factor (multiple levels)
  all[["df"]]$labeled_pheno <- as.factor(all[["df"]][[ttt]])
  # Convert labeled_pheno to an ordered factor
  all[["df"]][["labeled_pheno"]] <- factor(all[["df"]][["labeled_pheno"]], ordered = TRUE)

  # Create the formula object
  form <- reformulate(frm, response = "labeled_pheno")
  message(form)
  # Rescale the PRS
  #all[["df"]][[ptt]] <- scale(all[["df"]][[ptt]])

  if(length(levels(all[["df"]]$labeled_pheno)) == 2) {
    fit = glm(form, data = all[['df']], family=binomial(link="logit")) # THIS MAY NOT WORK
  }else{
    # requires mass package / library(mass)
    # BEWARE YOU NEED TO DOUBLE CHECK IF THE COMMANDS BELOW ACTUALLY FIT THE POLR MODEL
    tryCatch({
      fit = clm(form, data = all[["df"]], link = "logit", threshold = "flexible")
      }, error = function(e) {
        message('clm failed, trying polr')
      fit = polr(form, data = all[['df']], method = "logistic", Hess=TRUE)
      })
  }
  message('fit was done successfully')
  
  ## store table
  coefficients_df <- coef(summary(fit))
  coefficients_df = coefficients_df[!grepl('\\|', rownames(coefficients_df)),] # remove the 1|2 and 2|3 and etc..
  coefficients_df <- as.data.frame(coefficients_df, nrow = length(coefficients_df), byrow = FALSE)
  coefficients_df$predictors = rownames(coefficients_df)
  setnames(coefficients_df, 'Estimate', 'Value')
  
  message('add phenotype and model ID')
  # add phenotype and model ID
  coefficients_df$model_ID= basename(outDir)
  coefficients_df$phenotype = ttt

  message('calculate and store p values')
  if("t value" %in% names(coefficients_df)) coefficients_df$`Pr(>|t|)` <- pnorm(abs(coefficients_df[, "t value"]), lower.tail = FALSE) * 2
  coefficients_df = coefficients_df[,c(c('phenotype', 'model_ID', 'predictors'), setdiff(names(coefficients_df), c('model_ID','phenotype', 'predictors')))]
    
  crit_t <- qt(0.975, df.residual(fit))
 
  coefficients_df$deviance = deviance(fit)
  coefficients_df$n_genes = length(rownames(coefficients_df)[grepl('ENSG', rownames(coefficients_df))])
  
  # add n_pheno_levels
  coefficients_df$n_pheno_levels = length(unique(all[["df"]][[ttt]]))
  
  
  message('fit the null model.  Model with only an intercept (no predictors)')
  # fit the null model.  Model with only an intercept (no predictors)
  fit_null <- polr(labeled_pheno ~ 1,
                 data = all[["df"]],
                 method = "logistic")
  message('log-likelihood (also cost function of the model)')
  # log-likelihood (also cost function of the model) is the sum of the probability for each value
  # if you have 3 observations y (which are the correct ones) and for each the model assigns probability  0.1 , 0.2, 0.1
  # then the log-likelihood = log(0.1) + log(0.2) + log(0.1)
  ll_full <- logLik(fit)
  ll_null <- logLik(fit_null)
  
  no <- nobs(fit)  # Number of observations
  message('Calculate McFaddens R-squared')
  # Calculate McFadden's R-squared
  R2_McFadden <- 1 - (ll_full / ll_null)
  # Cox & Snell R²
  r2_cox_snell <- 1 - exp((2 / no) * (ll_null - ll_full))
  # Nagelkerke R² (Cragg & Uhler)
  r2_nagelkerke <- r2_cox_snell / (1 - exp((2 / no) * ll_null))

  coefficients_df$R2_McFadden <- R2_McFadden
  coefficients_df$R2_cox_snell <- r2_cox_snell
  coefficients_df$R2_nagelkerke <- r2_nagelkerke

  # Get AIC and BIC
  coefficients_df$AIC <- AIC(fit)
  coefficients_df$BIC <- BIC(fit)
  
  message('Extracted meaningful information')
  coef_df_file = paste0(outdir, "/flexible_logistic_coefficients_df_GReX_to_", gene.name, ".tsv")
  fwrite(coefficients_df, coef_df_file, sep = "\t", quote = F, row.names = F)
  return(coefficients_df)
  
  ###############################################################
  
  message('Define linearly dependent predictors (so that you can exclude them from your analysis)')
  # Define linearly dependent predictors (so that you can exclude them from your analysis)
  linear_dependent_predictors = setdiff(names(fit$coefficients), coefficients_df$predictors)

  ### Extract coefficients and standard errors
  estimates <- coefficients_df$Value
  std_errors <- coefficients_df$`Std. Error`

  # Calculate the critical t-value for 95% confidence intervals
  # Degrees of freedom = residual degrees of freedom from the model
  t_value <- qt(0.975, df.residual(fit))

  message("Calculate confidence intervals")
  # Lower and upper bounds of the confidence intervals
  coefficients_df$conf_int_low <- estimates - t_value * std_errors
  coefficients_df$conf_int_high <- estimates + t_value * std_errors
 
  message('odds ratios')
  # odds ratios
  coefficients_df$OR = exp(coefficients_df$Value)
  
  message('Get residual deviance')
  # Extract adjusted p-values
  if('z value' %in% names(coefficients_df)){
  coefficients_df$fdr.adjusted = p.adjust(coefficients_df$`Pr(>|z|)`, method = "fdr")
  coefficients_df$bonferroni.adjusted = p.adjust(coefficients_df$`Pr(>|z|)`, method = "bonferroni")
  coefficients_df$BH.adjusted = p.adjust(coefficients_df$`Pr(>|z|)`, method = "BH")
  
  # Get residual deviance
    ll <- logLik(fit)
    # Compute the deviance
    dev <- -2 * as.numeric(ll)
    coefficients_df$deviance = dev 
  }else if('t value' %in% names(coefficients_df)){
    
    coefficients_df$fdr.adjusted = p.adjust(coefficients_df$`Pr(>|t|)`, method = "fdr")
    coefficients_df$bonferroni.adjusted = p.adjust(coefficients_df$`Pr(>|t|)`, method = "bonferroni")
    coefficients_df$BH.adjusted = p.adjust(coefficients_df$`Pr(>|t|)`, method = "BH")
    
    # Compute the deviance
    coefficients_df$deviance = deviance(fit)
  }else{
     message('BE CAREFUL: NO t or z values in the coefficients_df')
  }


  message('add n_genes')
  # add n_genes
  coefficients_df$n_genes = length(rownames(coefficients_df)[grepl('ENSG', rownames(coefficients_df))])
  message('add n_pheno_levels')
  # add n_pheno_levels
  coefficients_df$n_pheno_levels = length(unique(all[["df"]][[ttt]]))
  message('fit the null model.  Model with only an intercept (no predictors)')
  # fit the null model.  Model with only an intercept (no predictors)
  fit_null <- polr(labeled_pheno ~ 1,
                 data = all[["df"]],
                 method = "logistic")
  message('log-likelihood (also cost function of the model)')
  # log-likelihood (also cost function of the model) is the sum of the probability for each value
  # if you have 3 observations y (which are the correct ones) and for each the model assigns probability  0.1 , 0.2, 0.1
  # then the log-likelihood = log(0.1) + log(0.2) + log(0.1)
  ll_full <- logLik(fit)
  ll_null <- logLik(fit_null)
  
  no <- nobs(fit)  # Number of observations
  message('Calculate McFaddens R-squared')
  # Calculate McFadden's R-squared
  R2_McFadden <- 1 - (ll_full / ll_null)
  # Cox & Snell R²
  r2_cox_snell <- 1 - exp((2 / no) * (ll_null - ll_full))
  # Nagelkerke R² (Cragg & Uhler)
  r2_nagelkerke <- r2_cox_snell / (1 - exp((2 / no) * ll_null))

  coefficients_df$R2_McFadden <- R2_McFadden
  coefficients_df$R2_cox_snell <- r2_cox_snell
  coefficients_df$R2_nagelkerke <- r2_nagelkerke

  # Get AIC and BIC
  coefficients_df$AIC <- AIC(fit)
  coefficients_df$BIC <- BIC(fit)

  # Get partitioned deviance
  mypredictors = coefficients_df$predictors
  #message('Now I will calculate partitioned Deviance (sapply(mypredictors,...))')
  # THIS IS RIDICULOUSLY EXPENSIVE LEAVE IT OUT UNLESS GV ASKS FOR IT
  #partitioned_deviance <- sapply(mypredictors, function(pred){
  
    # Fit the reduced model
  #  formula_reduced <- reformulate(paste(setdiff(mypredictors, pred), collapse = '+'), "labeled_pheno")
  #  fit_reduced <- polr(formula_reduced, data = all[["df"]], method = "logistic", Hess=TRUE)
    
    # Return the change in deviance
  #  deviance(fit_reduced) - coefficients_df$deviance[1] # this is the full deviance
 #})
  #message('I successfully calculated partitioned deviance.')
  #coefficients_df$partitioned_deviance <- partitioned_deviance

  coef_df_file = paste0(outdir, "/logistic_coefficients_df_GReX_to_", ttt, ".tsv")
  fwrite(coefficients_df, coef_df_file, sep = "\t", quote = F, row.names = F)

  return(coefficients_df)

  # DELETE MAYBE THE BELOW
  #############################################################
  #############################################################
  
  # Extract the necessary data
  #test <- list(
  #  "model.summary" = summary(fit),
  #  "PseudoR" = PseudoR2(fit),
  #  "calcVarPart" = calcVarPart(fit)
  #)
  
  # Save the results to a file
  #file.target  <- paste0(outdir, "/model_summary_GReX_to_", ttt, ".txt")
  #if (!dir.exists(dirname(file.target))) {
  #  dir.create(dirname(file.target), recursive=T)
  #}
  #zz <- file(file.target, "w")
  #sink(zz)
  #print(form)
  #print(test[["model.summary"]])
  #print(test[["PseudoR"]])
  #sink()
  #close(zz)

  ##################################################
  ##################################################
  # FOR LOGISTC REGRESSION
  # Extract coefficients and stats
  coefficients <- fit$coefficients # remove intercept
  coefficients_matrix <- matrix(coefficients, nrow = 1, ncol = length(coefficients), byrow = FALSE)
  coefficients_df <- as.data.frame(coefficients_matrix)
  colnames(coefficients_df) <- names(coefficients)
  coefficients_df$trait = ttt
  coefficients_df = coefficients_df[,c('trait', setdiff(names(coefficients_df), 'trait'))]
  
  # Get the residual variance
  coefficients_df$deviance = deviance(fit)

  # fit the null model.  Model with only an intercept (no predictors)
  fit_null <- polr(labeled_pheno ~ 1,
                 data = all[["df"]],
                 method = "logistic")


  ll_full <- logLik(fit)
  ll_null <- logLik(fit_null)
  
  no <- nobs(fit)  # Number of observations

  # Calculate McFadden's R-squared
  R2_McFadden <- 1 - (ll_full / ll_null)
  # Cox & Snell R²
  r2_cox_snell <- 1 - exp((2 / no) * (ll_null - ll_full))
  # Nagelkerke R² (Cragg & Uhler)
  r2_nagelkerke <- r2_cox_snell / (1 - exp((2 / no) * ll_null))

  coefficients_df$R2_McFadden <- R2_McFadden
  coefficients_df$R2_cox_snell <- r2_cox_snell
  coefficients_df$R2_nagelkerke <- r2_nagelkerke

  # Get AIC and BIC
  coefficients_df$AIC <- AIC(fit)
  coefficients_df$BIC <- BIC(fit)

  coef_df_file = paste0(outdir, "/logistic_coefficients_df_GReX_to_", ttt, ".predictions.tsv")
  fwrite(coefficients_df, coef_df_file, sep = "\t", quote = F, row.names = F)


  ################################################
  # Get partitioned deviance
  # Compute partitioned deviance

  mypredictors = colnames(coefficients_df)[grepl('ENSG|PC[1-9]|age|sex', colnames(coefficients_df))]
  mypredictors = setdiff(mypredictors, 'R2_nagelkerke')
  partitioned_deviance <- sapply(mypredictors, function(pred){
  
    # Fit the reduced model
    formula_reduced <- reformulate(paste(setdiff(mypredictors, pred), collapse = '+'), "labeled_pheno")
    
    fit_reduced <- polr(formula_reduced, data = all[["df"]], method = "logistic")
    
    # Return the change in deviance
    deviance(fit_reduced) - coefficients_df$deviance # this is the full deviance
  })

  partitioned_deviance_df <- data.frame(t(partitioned_deviance))
  partitioned_deviance_df$trait = ttt
  partitioned_deviance_df = partitioned_deviance_df[,c('trait', setdiff(names(partitioned_deviance_df), 'trait'))]

  fwrite(partitioned_deviance_df, paste0(outdir, "/logistic_partitioned_deviance_df_GReX_to_", ttt, ".tsv"), sep = "\t", quote = F, row.names = F)
  
  return(coefficients_df)
}
