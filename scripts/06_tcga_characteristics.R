
library(dplyr)

tcia_pt <- readxl::read_xlsx("./data/TCIA_TCGA-BRCA_09-16-2015-nbia-digest-1(1).xlsx")

counts_raw <- readRDS("data/counts_raw_tcga.RDS")

refined_data <- readRDS("./data/metadata_tcga.RDS") 

patient_dce_names <- read_csv("./output_data/dce_patient_names.csv") %>% 
  as.data.frame() %>% 
  names()

names <- unique(tcia_pt$`Patient ID`[tcia_pt$`Patient ID` %in% substr(start = 1, stop = 12, colnames(counts_raw))])

metadata <- refined_data %>% 
  filter(substr(start = 14, stop = 15, x = sampleID) == "01")


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

metadata_int_dce %>% ls
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



