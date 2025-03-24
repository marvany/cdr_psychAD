{            
            library(stringr)
            source("/sc/arion/projects/va-biobank/PROJECTS/ma_neurocog_grex/Scripts/GReX_to_pheno/acc_lsf_bash_RFunctions.R")
            source("/sc/arion/projects/va-biobank/PROJECTS/ma_neurocog_grex/Scripts/GReX_to_pheno/acc_resubmit_general.R")
            source("/sc/arion/projects/va-biobank/PROJECTS/ma_neurocog_grex/Scripts/GReX_to_pheno/accessoryScripts.R")
            create_lsf_command <- function(job_name, account, memory, walltime, cores, output_file, error_file, queue) {
                lsf_command <- paste0(
                    "bsub",
                    " -J ", job_name,
                    " -P ", account,
                    " -n ", cores,
                    " -R span[hosts=1]",
                    " -R rusage[mem=", memory, "]",
                    " -W ", walltime,
                    " -o ", output_file,
                    " -e ", error_file,
                    " -q ", queue,
                    " -L /bin/bash" 
                )
                
                return(lsf_command)
            }

            generate_sh_script <- function(script_path, output_sh_path, language, my_env = NA, ...) {
                if(language == 'python' && is.na(my_env)) stop("Please specify variable 'my_env'") 
                sh_content <- paste(
                    "#!/bin/bash",
                    paste0("cd ", getwd()),
                    ifelse(language == "R", "ml R", paste0("module load anaconda3\nsource activate ", my_env)),
                    ifelse(language == "R", paste("Rscript", script_path, ...), paste("python", script_path, ...)),
                    sep = "\n"
                )
                
                writeLines(sh_content, con = output_sh_path)
                system(paste("chmod +x", output_sh_path))
            }


            # calculate_memory <- function(file_size, alpha = 800, beta = 20, gamma = 0.97, k = 3000) {
                    # Formula: resources = α * log(file_size + 1) + β * file_size
            #        resources <- alpha * log(file_size + k) + beta * (file_size ^ gamma)

            #        return(resources)
            #}

            calculate_memory <- function(file_size, cutoff, alpha, beta){
            resources = alpha*file_size + beta
            if(resources > cutoff) resources = 500000
            return(resources)
            }
}

#####################
# Define paths
# INDEPENDENT
working.dir = "/sc/arion/projects/va-biobank/PROJECTS/ma_cdr_psychAD"
submission.dir = paste0(working.dir, '/drugs_pheno_assoc_anal') # this is where .err and .out will be stored
path_to_script = '/sc/arion/projects/va-biobank/PROJECTS/ma_cdr_psychAD/Scripts/2.1.2drug_phenotypes_assoc_analysis.R' # R file


#target_files = paste0(working.dir, '/Resources/GReX_to_pheno/UKBB/grex_pheno_comb/')
#target_files_pattern = 'filtered.txt.gz' # will be used to grep the files in target directory

# DEPENDENT (don't modify)
jobs_dir = paste0(submission.dir, '/jobs/')
logs_output_dir = paste0(submission.dir, '/logs')


#####################
# Create Directories
setwd(working.dir)
if(!dir.exists(logs_output_dir)) dir.create(logs_output_dir, recursive = TRUE)
if(!dir.exists(jobs_dir)) dir.create(jobs_dir, recursive = TRUE)

#####################
# Command parameters:
#job_name <- 'filter_NA' # filter_NA_grex_pheno
job_name <- 'drugs_pheno_asan'
account <- 'acc_va-biobank'
memory <- '3000' #calculated internally # 1 unit is 1 MB, 1000 is 1 GB, 10000 is 10 GB
cores <- '40'
walltime <- '72:00'

queue <- 'premium'

# MAKE SURE TO MODIFY output_file and error_file as well

####################################################################################
# !!!!!!!!!!!!!!!!!!! BEWARE !!!!!!!!!!!!!!!!!!! Define files that will be submitted 
#these_files = list.files(target_files, full.names = T)[grep(target_files_pattern, list.files(target_files))]
# this is recursive!!
#these_files = list.files(target_files, full.names = T, recursive = T)[grep(target_files_pattern, list.files(target_files))]




adjust_to_memory = F
# SUBMIT JOBS
#for(this_file in these_files){

if(adjust_to_memory){ # we request for 10x memory than the memory of the file
    full_file_path = these_files[grep(this_file, these_files)]
    memory = file.info(full_file_path)$size / (1024^2) # we want it in megabytes
    memory = calculate_memory(memory, alpha = 20, beta = 15000, cutoff = 500000)
}

output_file <- paste0(logs_output_dir, '/', 'assoc_anal', '_output.out')
error_file <- paste0(logs_output_dir, '/', 'assoc_anal', '_error.err')
if(!file.exists(dirname(output_file))) dir.create(dirname(output_file), recursive = TRUE)
if(!file.exists(dirname(error_file))) dir.create(dirname(error_file), recursive = TRUE)

# Create LSF command
bsub_command <- create_lsf_command(job_name,
                    account = account,
                    memory = memory,
                    cores = cores,
                    walltime = walltime,
                    output_file = output_file,
                    error_file = error_file,
                    queue = queue)

# Generate .sh script that will execute the Rscript
path_to_submission_file = paste0(jobs_dir, job_name, '.sh')
generate_sh_script(
    script_path = path_to_script, # /path/to/file.R which will be executed
    output_sh_path = path_to_submission_file, # Directory where .sh file will be stored
    language = 'R'
    )
to_submit <- paste0(bsub_command, ' ', path_to_submission_file)
job_id = system(to_submit, intern = TRUE)
job.id = gsub('queue.*', '', gsub('.*<', '', gsub('>.*', '', job_id)))
job.id
Sys.time()
#}
# 168940173 "2025-02-09 19:49:32 EST"
