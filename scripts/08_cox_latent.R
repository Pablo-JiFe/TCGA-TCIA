library(dplyr)
library(tidyr)
library(tibble)
library(survival)
library(glmnet)
library(tidymodels)
library(censored)

train_latent <- read.csv("output_data/latent_train_vae.csv", row.names = 1)

test_latent <- read.csv("output_data/latent_test_vae.csv", row.names = 1)

refined_data <- readRDS("./data/metadata_tcga.RDS")

metadata <-
  refined_data %>%
  filter(substr(start = 14, stop = 15, x = sampleID) == "01") # Filter only the registry of the sample corresponding to primary






metadata_train <- metadata[metadata$sampleID %in% substr(start = 1,
                                                         stop = 15,
                                                         x = rownames(train_latent)), ]


train_surv <-
  train_latent %>%
  mutate(sampleID = substr(
    start = 1,
    stop = 15,
    x = rownames(train_latent)
  )) %>%
  left_join(metadata_train %>%
              dplyr::select(RECURRENCE, RECURRENCE_MON, sampleID),
            by = "sampleID") %>%
  column_to_rownames("sampleID") %>%
  drop_na(RECURRENCE_MON)

train_surv <- train_surv %>%
  mutate(RECURRENCE_MON = pmax(RECURRENCE_MON, 1e-3))

train_surv$surv_obj <- Surv(train_surv$RECURRENCE_MON, train_surv$RECURRENCE)




metadata_test <- metadata[metadata$sampleID %in% substr(start = 1,
                                                        stop = 15,
                                                        x = rownames(test_latent)), ]

test_surv <-
  test_latent %>%
  mutate(sampleID = substr(
    start = 1,
    stop = 15,
    x = rownames(test_latent)
  )) %>%
  left_join(metadata_test %>%
              dplyr::select(RECURRENCE, RECURRENCE_MON, sampleID),
            by = "sampleID") %>%
  column_to_rownames("sampleID")


test_surv$surv_obj <- Surv(test_surv$RECURRENCE_MON, test_surv$RECURRENCE)



# Penalized ---------------------------------------------------------------


lr_rec <-
  recipe(surv_obj ~ . ,data = train_surv) %>%
  update_role(RECURRENCE, RECURRENCE_MON, new_role = "Result") %>% 
  step_zv(all_predictors()) %>%
  step_nzv(all_predictors())


lr_mod <- proportional_hazards(
  penalty = tune(),
  mixture = tune()     
) %>%
  set_engine("glmnet", cox.ties = "breslow") 


lr_wf <- workflow() %>%
  add_model(lr_mod) %>%
  add_recipe(lr_rec)


set.seed(123)

# 3.1 Folds for evaluating with resamples

folds <- vfold_cv(
  train_surv,
  v = 10,
  strata = RECURRENCE
)

# 3.2 Grid for penalizing range

grid <- grid_regular(
  penalty(range = c( - 4, 1)),   
  mixture(range = c(0, 1)),
  levels = 10
)

# 3.3 Running the different penalization methods

res_ml <- tune_grid(
  lr_wf,
  resamples = folds,
  grid = grid,
  metrics = metric_set(concordance_survival), # Evaluates the different penalizing methods by c score
  control = control_grid(save_pred = TRUE)
)


collect_metrics(res_ml)

# 3.3.2 Object with best parameters for penalizing

best_params <- select_best(res_ml, metric = "concordance_survival")

# 4.- Actual training -----------------------------------------------------

# 4.1 Final workflow with the selecting the best parameter tested previously

final_wf <- finalize_workflow(lr_wf, best_params)

# 4.2 Final fit with training data

final_fit <- fit(final_wf, data = train_surv)

# 4.2.2 Observing genes that are maintained after penalziation

coef_tbl <- tidy(final_fit) %>%
  filter(estimate != 0) %>%
  arrange(desc(abs(estimate)))

cat(coef_tbl$term, sep = ", ")


st_pred <- predict(final_fit, new_data = test_surv, type = "linear_pred")

c_results <- concordance(surv_obj ~ st_pred$.pred_linear_pred, data = test_surv, reverse = FALSE)

c_results

# 5.- Proportional Hazards Diagnostics ------------------------------------

# 5.1 Extract the preprocessed training data from recipe

prepped_recipe <- extract_recipe(final_fit)
processed_train <- bake(prepped_recipe, new_data = train_surv)

# 5.2 names of the features 
active_genes <- coef_tbl %>% 
  filter(term != "penalty") %>% # Exclude penalty metadata if present
  pull(term)

# 5.3 formula matching predictors
diag_formula <- as.formula(
  paste("surv_obj ~", paste(active_genes, collapse = " + "))
)

# 5.4 Fit an unpenalized Cox model on the selected variables
diag_model <- coxph(diag_formula, data = processed_train)

# 5.5 Run the proportional hazards test and plot
ph_test <- cox.zph(diag_model)
print(ph_test)
