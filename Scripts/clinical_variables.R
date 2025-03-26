

# These are the variables we will use for this analysis
general <- c(
 "SubID",    # Unique identifier for each donor, prefixed by institution letter (M=MSSM, H=HBCC, R=RUSH)
 "Age",      # Age at death
 "Sex",      # Biological sex
 "Ethnicity", # Racial/ethnic background
 "Dx",       # Text field with diagnosis information (format varies by brain bank)
 "pH",       # Brain tissue pH value at collection, affects RNA quality and tissue preservation
 "PMI",      # Post-Mortem Interval: Time between death and brain tissue preservation (hours)
 "Source_Location"  # Anatomical location within the brain where the tissue sample was collected
)

alzheimer_vars <- c(
 "AD",                  # Alzheimer's disease
 "MCI",                 # Mild Cognitive Impairment
 "Dementia",            # Cognitive decline beyond normal aging
 "Plq_Mn",              # Mean plaque density score, overall amyloid plaque burden
 "Plq_Mn_MFG",          # Mean plaque density score in middle frontal gyrus
 "ArgyrophilicGrain",   # Argyrophilic Grain Disease
 "CerebralAtrophy",     # Brain tissue loss
 "Tauopathy",           # Tau protein accumulation disorders
 "ApoE_gt",             # ApoE genotype (Alzheimer's risk gene)
 "CERAD",               # Consortium to Establish a Registry for AD score (1-4)
 "BRAAK_AD",            # Braak staging for Alzheimer's pathology (0-6)
 "BRAAK_PD",            # Braak staging for Parkinson's pathology
 "CDRScore",            # Clinical Dementia Rating (overall score)
 "Others_Neurodegenerative"  # Other neurodegenerative conditions
)

cdr_vars <- c(
 "CDR_Memory",                    # Memory domain score
 "CDR_Orientation",               # Orientation domain score
 "CDR_Judgement",                 # Judgment/problem solving domain score
 "CDR_Community",                 # Community affairs domain score
 "CDR_HomeHobbies",               # Home and hobbies domain score
 "CDR_PersonalCare",              # Personal care domain score
 "CDR_SumBoxes",                  # Sum of all domain scores
 "Cognitive_Resilience",          # Measure of cognitive function relative to pathology
 "Cognitive_and_Tau_Resilience"   # Cognitive function relative to tau pathology
)

plaque_vars <- c(
 "HippoPlaquesValue",         # Amyloid plaque density in hippocampus
 "HippoPlaquesWCoresValue",   # Amyloid plaques with cores in hippocampus
 "EntorPlaquesValue",         # Amyloid plaque density in entorhinal cortex
 "EntorPlaquesWCoresValue",   # Amyloid plaques with cores in entorhinal cortex
 "AmygPlaquesValue",          # Amyloid plaque density in amygdala
 "AmygPlaquesWCoresValue",    # Amyloid plaques with cores in amygdala
 "MidPlaquesValue",           # Amyloid plaque density in middle frontal cortex
 "MidPlaquesWCoresValue",     # Amyloid plaques with cores in middle frontal cortex
 "SupPlaquesValue",           # Amyloid plaque density in superior temporal cortex
 "SupPlaquesWCoresValue",     # Amyloid plaques with cores in superior temporal cortex
 "InfPlaquesValue",           # Amyloid plaque density in inferior parietal cortex
 "InfPlaquesWCoresValue",     # Amyloid plaques with cores in inferior parietal cortex
 "OcciPlaquesValue",          # Amyloid plaque density in occipital cortex
 "OcciPlaquesWCoresValue"     # Amyloid plaques with cores in occipital cortex
)

# Neurofibrillary tangles
tangle_vars <- c(
 "HippoTanglesValue",    # Neurofibrillary tangles in hippocampus
 "EntorTanglesValue",    # Neurofibrillary tangles in entorhinal cortex
 "AmygTanglesValue",     # Neurofibrillary tangles in amygdala
 "MidTanglesValue",      # Neurofibrillary tangles in middle frontal cortex
 "MidGliosisValue"       # Gliosis in middle frontal cortex
)

gliosis_vars <- c(
 "SupGliosisValue",    # Gliosis in superior temporal cortex
 "InfGliosisValue",    # Gliosis in inferior parietal cortex
 "OcciGliosisValue"    # Gliosis in occipital cortex
)



neuro_vars <- c(
 "Sex_chr_aneuploidy",            # Boolean indicating sex chromosome abnormalities
 "NormPressHydrocephalus",        # Normal pressure hydrocephalus
 "MS",                            # Multiple sclerosis
 "PSP",                           # Progressive supranuclear palsy
 "Epilepsy",                      # Seizure disorder
 "Seizures",                      # Seizure events without epilepsy diagnosis
 "Migraine_headaches",            # Migraine condition
 "Vascular",                      # Vascular brain pathology
 "Others_Neurological",           # Other neurological conditions
 "Leucotomy",                     # Surgical white matter cutting
 "PD",                            # Parkinson's disease
 "PD_uncertain_plus_encephalitic", # Uncertain Parkinson's cases and encephalitic parkinsonism
 "DLBD",                          # Diffuse Lewy body disease
 "ALS",                           # Amyotrophic lateral sclerosis
 "Others_Neurodegenerative"       # Other neurodegenerative conditions
)

psych_vars <- c(
 "SCZ",                       # Schizophrenia
 "MDD",                       # Major depressive disorder
 "BD_unspecific",             # Bipolar disorder (unspecified type)
 "BD_I",                      # Bipolar disorder type I
 "BD_II",                     # Bipolar disorder type II
 "FTD",                       # Frontotemporal dementia
 "PTSD",                      # Post-traumatic stress disorder
 "ADHD",                      # Attention deficit hyperactivity disorder
 "OCD",                       # Obsessive-compulsive disorder
 "Schizoaffective_bipolar",   # Schizoaffective disorder bipolar type
 "Schizoaffective_depressive", # Schizoaffective disorder depressive type
 "Anorexia",                  # Anorexia nervosa
 "Bulimia",                   # Bulimia nervosa
 "Binge_Purge",               # Binge-purge behaviors
 "Eating_disorder",           # Other eating disorders
 "Anxiety",                   # Anxiety disorders
 "ASHCVD",                    # Atherosclerotic cardiovascular disease
 "TD_I",                      # Type 1 diabetes
 "TD_II",                     # Type 2 diabetes
 "Others_Neuropsychiatric"    # Other psychiatric conditions
)
