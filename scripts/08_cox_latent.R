library(dplyr)
library(survival)

train_latent <- read.csv("output_data/latent_train.csv", row.names = 1)

test_latent <- read.csv("output_data/latent_test.csv", row.names = 1)

refined_data <- readRDS("./data/metadata_tcga.RDS") 

metadata <- 
  refined_data %>% 
  filter(substr(start = 14, stop = 15, x = sampleID) == "01") # Filter only the registry of the sample corresponding to primary






metadata_train <- metadata[metadata$sampleID %in% substr(start = 1, stop = 15, x = rownames(train_latent)),]


train_surv <- 
  train_latent %>% 
  mutate(sampleID = substr(start = 1, stop = 15, x = rownames(train_latent))) %>% 
  left_join(metadata_train %>% 
              dplyr::select(RECURRENCE, RECURRENCE_MON, sampleID), by = "sampleID") %>% 
  column_to_rownames("sampleID") %>% 
  drop_na(RECURRENCE_MON)

surv_obj <- Surv(train_surv$RECURRENCE_MON, train_surv$RECURRENCE)

train_surv_for_cox <- train_surv %>% 
  dplyr::select(-c(RECURRENCE_MON, RECURRENCE))

cox_model <- coxph(surv_obj~., train_surv_for_cox)
summary(cox_model)






metadata_test <- metadata[metadata$sampleID %in% substr(start = 1, stop = 15, x = rownames(test_latent)),]

test_surv <- 
  test_latent %>% 
  mutate(sampleID = substr(start = 1, stop = 15, x = rownames(test_latent))) %>% 
  left_join(metadata_test %>% 
              dplyr::select(RECURRENCE, RECURRENCE_MON, sampleID), by = "sampleID") %>% 
  column_to_rownames("sampleID")


surv_obj_test <- Surv(test_surv$RECURRENCE_MON, test_surv$RECURRENCE)

test_surv_for_cox <- test_surv %>% 
  dplyr::select(-c(RECURRENCE_MON, RECURRENCE))

cox_test <- coxph(surv_obj_test~., test_surv_for_cox)
summary(cox_test)
cox.zph(cox_test)



test_surv$risk <- predict(cox_model, newdata = test_surv, type = "risk")

c_results <- concordance(surv_obj_test ~ risk, data = test_surv, reverse = TRUE)

