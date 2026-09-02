library(TCGAbiolinks)
library(SummarizedExperiment)
library(UCSCXenaTools)
library(dplyr)

tcga_rna <- GDCquery("TCGA-BRCA",
                     data.category = "Transcriptome Profiling",
                     access = "open",
                     experimental.strategy = "RNA-Seq",
                     workflow.type = "STAR - dplyr::counts"
)

GDCdownload(tcga_rna, method = "api", files.per.chunk = 5)


tcga_brca_data <- GDCprepare(tcga_rna, directory = "d:/tcga/GDCdata/")


brca_matrix <- assay(tcga_brca_data, "fpkm_unstrand")

# 1.- Metadata ------------------------------------------------------------



# 1.1 Generate and Query

data_query <- XenaGenerate(subset = XenaDatasets == "TCGA.BRCA.sampleMap/BRCA_clinicalMatrix") %>% 
  XenaQuery()

# 1.2 Download 

xe_download <- XenaDownload(data_query, destdir = "D:/tcga/GDCdata/Metadata")

# 1.3 Prepare (Load) the data

brca_clinical <- XenaPrepare(xe_download)


# TCIA metadata

metadata_tcia <- readxl::read_xlsx("./data/TCIA_TCGA-BRCA_09-16-2015-nbia-digest-1(1).xlsx")

# Obtain names of patients from tcia

names <- unique(metadata_tcia$`Patient ID`)

# 1.4 Create the Recurrence variables


refined_data <- 
  brca_clinical %>%
  mutate(
    SURVIVAL = ifelse(vital_status == "DECEASED", 1, 0),
    SURVIVAL_MON = ifelse(!is.na(days_to_death) & days_to_death >= 0,
                          days_to_death, # If they died they have this parameter
                          days_to_last_followup # Else they have this parameter
    ) / 30.4166667,
    RECURRENCE = ifelse(
      new_neoplasm_event_type %in% c("Locoregional Recurrence", "Distant Metastasis"),
      1,
      0
    ),
    RECURRENCE_MON = ifelse(
      RECURRENCE == 1,
      days_to_new_tumor_event_after_initial_treatment,
      coalesce(days_to_last_followup, days_to_death)
    ) / 30.4166667,
    Patient_ID = substr(start = 1, stop = 12, x = sampleID)
  ) 

metadata <- 
  refined_data %>% 
  filter(substr(start = 14, stop = 15, x = sampleID) == "01",
         Patient_ID %in% substr(start = 1, stop = 12, x = colnames(brca_matrix))) # Filter only the registry of the sample corresponding to primary


dim(metadata)


metadata %>% 
  dplyr::count(RECURRENCE) %>% 
  mutate(prop = n / sum(n))

# How many patients have normal adjacent tissue replicates

refined_data %>% 
  filter(Patient_ID %in% names,
         substr(start = 14, stop = 15, x = sampleID) == "11") %>% 
  dplyr::count(Patient_ID) %>% 
  arrange(desc(n))

# How many patients have metastasis samples

refined_data %>% 
  filter(Patient_ID %in% names,
         substr(start = 14, stop = 15, x = sampleID) == "06") %>% 
  dplyr::count(Patient_ID) %>% 
  arrange(desc(n))

# How many patients have normal adjacent tissue replicates and metastasis

refined_data %>% 
  filter(Patient_ID %in% names,
         substr(start = 14, stop = 15, x = sampleID) == "11" | substr(start = 14, stop = 15, x = sampleID) == "06") %>% 
  dplyr::count(Patient_ID) %>% 
  filter(n >= 2) %>% 
  arrange(desc(n))

# Keep only patients available in TCIA metadata

metadata_int <- metadata %>% 
  filter(Patient_ID %in% names)

# N of patients anbd general demographic

dim(metadata_int)

summary(metadata_int$age_at_initial_pathologic_diagnosis)


metadata_int %>% 
  filter(age_at_initial_pathologic_diagnosis < 40)

summary(metadata_int$age_at_initial_pathologic_diagnosis)

metadata_int %>% 
  dplyr::count(gender)

summary(metadata_int[metadata_int$SURVIVAL == 0, ]$SURVIVAL_MON)
summary(metadata_int[metadata_int$SURVIVAL == 1, ]$SURVIVAL_MON)

# N and prop of survival

metadata_int %>% 
  dplyr::count(SURVIVAL) %>% 
  mutate(prop = n / sum(n))

metadata_int %>% 
  filter(SURVIVAL == 0) %>% 
  summarise(
    mean_follow_up = mean(SURVIVAL_MON),
    min_follow_up = min(SURVIVAL_MON),
    max_follow_up = max(SURVIVAL_MON),
    sd_follow_up = sd(SURVIVAL_MON)
  )

metadata_int %>% 
  filter(SURVIVAL == 1) %>% 
  summarise(
    mean_follow_up = mean(SURVIVAL_MON),
    min_follow_up = min(SURVIVAL_MON),
    max_follow_up = max(SURVIVAL_MON),
    sd_follow_up = sd(SURVIVAL_MON)
  )

# N and prop of recurrence

metadata_int %>% 
  dplyr::count(RECURRENCE) %>% 
  mutate(prop = n / sum(n))

metadata_int %>% 
  filter(RECURRENCE == 0) %>% 
  summarise(
    mean_follow_up = mean(RECURRENCE_MON),
    min_follow_up = min(RECURRENCE_MON),
    max_follow_up = max(RECURRENCE_MON),
    sd_follow_up = sd(RECURRENCE_MON)
  )

metadata_int %>% 
  filter(RECURRENCE == 1) %>% 
  summarise(
    mean_follow_up = mean(RECURRENCE_MON),
    min_follow_up = min(RECURRENCE_MON),
    max_follow_up = max(RECURRENCE_MON),
    sd_follow_up = sd(RECURRENCE_MON)
  )


metadata_int %>% 
  dplyr::count(SURVIVAL,
        RECURRENCE)

# N and prop of Hist sub

metadata_int %>% 
  dplyr::count(histological_type) %>% 
  mutate(prop = n / sum(n))

metadata_int %>% 
  dplyr::count(histological_type_other) %>% 
  mutate(prop = n / sum(n))

metadata_int %>% 
  dplyr::group_by(histological_type) %>%
  dplyr::count(SURVIVAL) %>% 
  mutate(prop = n / sum(n))

metadata_int %>% 
  dplyr::group_by(histological_type) %>%
  dplyr::count(RECURRENCE) %>% 
  mutate(prop = n / sum(n))

# N and prop of PAM50

metadata_int %>% 
  dplyr::count(PAM50Call_RNAseq) %>% 
  mutate(prop = n / sum(n))

metadata_int %>% 
  dplyr::group_by(PAM50Call_RNAseq) %>%
  dplyr::count(SURVIVAL) %>% 
  mutate(prop = n / sum(n))

metadata_int %>% 
  dplyr::group_by(PAM50Call_RNAseq) %>%
  dplyr::count(RECURRENCE) %>% 
  mutate(prop = n / sum(n))


# N and prop of path stage

metadata_int %>% 
  dplyr::count(pathologic_stage) %>% 
  mutate(prop = n / sum(n))

metadata_int %>% 
  dplyr::group_by(pathologic_stage) %>%
  dplyr::count(SURVIVAL) %>% 
  mutate(prop = n / sum(n))

metadata_int %>% 
  dplyr::group_by(pathologic_stage) %>%
  dplyr::count(RECURRENCE) %>% 
  mutate(prop = n / sum(n))


# N and prop of ER

metadata_int %>% 
  dplyr::count(ER_Status_nature2012, HER2_Final_Status_nature2012) %>% 
  mutate(prop = n / sum(n))

metadata_int %>% 
  dplyr::group_by(ER_Status_nature2012, HER2_Final_Status_nature2012) %>%
  dplyr::count(SURVIVAL) %>% 
  mutate(prop = n / sum(n))

metadata_int %>% 
  dplyr::group_by(ER_Status_nature2012, HER2_Final_Status_nature2012) %>%
  dplyr::count(RECURRENCE) %>% 
  mutate(prop = n / sum(n))





patient_dce_names <- scan("output_data/dce_patient_names.csv", sep = ",", what = character())


metadata_int_dce <- 
  metadata_int %>% 
  filter(Patient_ID %in% patient_dce_names)

summary(metadata_int_dce$age_at_initial_pathologic_diagnosis)


metadata_int_dce %>% 
  filter(age_at_initial_pathologic_diagnosis < 40)

summary(metadata_int_dce$age_at_initial_pathologic_diagnosis)

metadata_int_dce %>% 
  dplyr::count(gender)

# N and prop of survival

metadata_int_dce %>% 
  dplyr::count(SURVIVAL) %>% 
  mutate(prop = n / sum(n))

metadata_int_dce %>% 
  filter(SURVIVAL == 0) %>% 
  summarise(
    mean_follow_up = mean(SURVIVAL_MON),
    min_follow_up = min(SURVIVAL_MON),
    max_follow_up = max(SURVIVAL_MON),
    sd_follow_up = sd(SURVIVAL_MON)
  )

metadata_int_dce %>% 
  filter(SURVIVAL == 1) %>% 
  summarise(
    mean_follow_up = mean(SURVIVAL_MON),
    min_follow_up = min(SURVIVAL_MON),
    max_follow_up = max(SURVIVAL_MON),
    sd_follow_up = sd(SURVIVAL_MON)
  )

# N and prop of recurrence

metadata_int_dce %>% 
  dplyr::count(RECURRENCE) %>% 
  mutate(prop = n / sum(n))

metadata_int_dce %>% 
  filter(RECURRENCE == 0) %>% 
  summarise(
    mean_follow_up = mean(RECURRENCE_MON),
    min_follow_up = min(RECURRENCE_MON),
    max_follow_up = max(RECURRENCE_MON),
    sd_follow_up = sd(RECURRENCE_MON)
  )

metadata_int_dce %>% 
  filter(RECURRENCE == 1) %>% 
  summarise(
    mean_follow_up = mean(RECURRENCE_MON),
    min_follow_up = min(RECURRENCE_MON),
    max_follow_up = max(RECURRENCE_MON),
    sd_follow_up = sd(RECURRENCE_MON)
  )


metadata_int_dce %>% 
  dplyr::count(SURVIVAL,
               RECURRENCE)

# N and prop of Hist sub

metadata_int_dce %>% 
  dplyr::count(histological_type) %>% 
  mutate(prop = n / sum(n))

metadata_int_dce %>% 
  dplyr::count(histological_type_other) %>% 
  mutate(prop = n / sum(n))

metadata_int_dce %>% 
  dplyr::group_by(histological_type) %>%
  dplyr::count(SURVIVAL) %>% 
  mutate(prop = n / sum(n))

metadata_int_dce %>% 
  dplyr::group_by(histological_type) %>%
  dplyr::count(RECURRENCE) %>% 
  mutate(prop = n / sum(n))

# N and prop of PAM50

metadata_int_dce %>% 
  dplyr::count(PAM50Call_RNAseq) %>% 
  mutate(prop = n / sum(n))

metadata_int_dce %>% 
  dplyr::group_by(PAM50Call_RNAseq) %>%
  dplyr::count(SURVIVAL) %>% 
  mutate(prop = n / sum(n))

metadata_int_dce %>% 
  dplyr::group_by(PAM50Call_RNAseq) %>%
  dplyr::count(RECURRENCE) %>% 
  mutate(prop = n / sum(n))


# N and prop of path stage

metadata_int_dce %>% 
  dplyr::count(pathologic_stage) %>% 
  mutate(prop = n / sum(n))

metadata_int_dce %>% 
  dplyr::group_by(pathologic_stage) %>%
  dplyr::count(SURVIVAL) %>% 
  mutate(prop = n / sum(n))

metadata_int_dce %>% 
  dplyr::group_by(pathologic_stage) %>%
  dplyr::count(RECURRENCE) %>% 
  mutate(prop = n / sum(n))

# N and prop of ER

metadata_int_dce %>% 
  dplyr::count(ER_Status_nature2012, HER2_Final_Status_nature2012) %>% 
  mutate(prop = n / sum(n))

metadata_int_dce %>% 
  dplyr::group_by(ER_Status_nature2012, HER2_Final_Status_nature2012) %>%
  dplyr::count(RECURRENCE) %>% 
  mutate(prop = n / sum(n))

metadata_int_dce %>% 
  dplyr::group_by(ER_Status_nature2012, HER2_Final_Status_nature2012) %>%
  dplyr::count(RECURRENCE) %>% 
  mutate(prop = n / sum(n))



