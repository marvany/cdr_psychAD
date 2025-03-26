vaedir = "/sc/arion/projects/va-biobank/jamie/cdr/autoencoder/results/evaluation/test/bipolar-w23me-normrank/marios"
list.files(vaedir)
smc = fread("/sc/arion/projects/va-biobank/jamie/cdr/autoencoder/results/evaluation/test/bipolar-w23me-normrank/marios/EUR_MegaAnalysis_subclass_SMC_X.csv.gz")
head(smc)
a = fread('/sc/arion/projects/va-biobank/PROJECTS/marios_temp/parallel_antagonism/working.directories/bd_1_2/V1/results/GTP_CDR/intermediate.files/avgRank/results/meta.bip_sab.eur_w23me.unadjust/trt_cp_EUR_MegaAnalysis_subclass_SMC_X.signatures.AvgRank.csv.gz')
head(a)
