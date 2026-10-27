library(dplyr)
library(ggplot2)
library(bench)
library(arrow)
library(limma)

counts_raw <- readRDS("data/counts_raw_tcga.RDS")
counts_vst <- read_parquet("./data/counts_vst_un.parquet")
refined_data <- readRDS("./data/metadata_tcga.RDS")
genes_tcga <- readRDS("./data/genes_tcga.RDS")
fpkm_counts <- readRDS("./data/counts_fpkm.RDS") %>%
  as.matrix()
counts_raw_GTEx <- readRDS("./data/counts_df_gte.RDS") %>%
  tidyr::drop_na(Symbol) %>%
  tibble::rownames_to_column("prev_id") %>%
  dplyr::select(-prev_id) %>%
  tibble::column_to_rownames("Symbol")


counts_raw <-
  counts_raw %>%
  tibble::column_to_rownames("gene_name")

counts_raw <- counts_raw[rownames(counts_raw_GTEx), ]

merged <-
  cbind(counts_raw_GTEx, counts_raw)

merged_metadata <-
  data.frame(patients = colnames(merged)) %>%
  mutate(cohort = ifelse(substr(
    start = 1, stop = 4, patients
  ) == "TCGA", "TCGA", "GTEx"))


# VST transform -----------------------------------------------------------
counts_raw <- as.matrix(counts_raw)

res <- bench::mark(
  n_100_5k  = DESeq2::vst(counts_raw[1:5000, 1:100]),
  n_200_5k  = DESeq2::vst(counts_raw[1:5000, 1:200]),
  n_500_5k  = DESeq2::vst(counts_raw[1:5000, 1:500]),
  n_1000_5k  = DESeq2::vst(counts_raw[1:5000, 1:1000]),
  n_100_10k  = DESeq2::vst(counts_raw[1:10000, 1:100]),
  n_200_10k  = DESeq2::vst(counts_raw[1:10000, 1:200]),
  n_500_10k  = DESeq2::vst(counts_raw[1:10000, 1:500]),
  n_1000_10k  = DESeq2::vst(counts_raw[1:10000, 1:1000]),
  n_100_30k  = DESeq2::vst(counts_raw[1:30000, 1:100]),
  n_200_30k  = DESeq2::vst(counts_raw[1:30000, 1:200]),
  n_500_30k  = DESeq2::vst(counts_raw[1:30000, 1:500]),
  n_1000_30k  = DESeq2::vst(counts_raw[1:30000, 1:1000]),
  n_100_50k  = DESeq2::vst(counts_raw[1:50000, 1:100]),
  n_200_50k  = DESeq2::vst(counts_raw[1:50000, 1:200]),
  n_500_50k  = DESeq2::vst(counts_raw[1:50000, 1:500]),
  n_1000_50k  = DESeq2::vst(counts_raw[1:50000, 1:1000]),
  iterations = 1,
  check = FALSE
)

res_tab <- bind_cols(tibble(expresion = names(res$expression)), res[, 2:9])


counts_raw_GTEx <- as.matrix(counts_raw_GTEx)

res_vst_gtex <- bench::mark(
  n_100_5k  = DESeq2::vst(counts_raw_GTEx[1:5000, 1:100]),
  n_178_5k  = DESeq2::vst(counts_raw_GTEx[1:5000, 1:178]),
  n_100_10k  = DESeq2::vst(counts_raw_GTEx[1:10000, 1:100]),
  n_178_10k  = DESeq2::vst(counts_raw_GTEx[1:10000, 1:178]),
  n_100_30k  = DESeq2::vst(counts_raw_GTEx[1:30000, 1:100]),
  n_178_30k  = DESeq2::vst(counts_raw_GTEx[1:30000, 1:178]),
  n_100_50k  = DESeq2::vst(counts_raw_GTEx[1:50000, 1:100]),
  n_178_50k  = DESeq2::vst(counts_raw_GTEx[1:50000, 1:178]),
  iterations = 1,
  check = FALSE
)

res_vst_gtex_tab <- bind_cols(tibble(expresion = names(res_vst_gtex$expression)), res_vst_gtex[, 2:9])


# Scale -------------------------------------------------------------------

counts_vst <-
  counts_vst %>%
  tibble::column_to_rownames("gene_name") %>%
  as.matrix()

res_scale <- bench::mark(
  n_100_5k  = scale(counts_vst[1:5000, 1:100]),
  n_200_5k  = scale(counts_vst[1:5000, 1:200]),
  n_500_5k  = scale(counts_vst[1:5000, 1:500]),
  n_1000_5k = scale(counts_vst[1:5000, 1:1000]),
  n_100_10k   = scale(counts_vst[1:10000, 1:100]),
  n_200_10k   = scale(counts_vst[1:10000, 1:200]),
  n_500_10k   = scale(counts_vst[1:10000, 1:500]),
  n_1000_10k  = scale(counts_vst[1:10000, 1:1000]),
  n_100_30k   = scale(counts_vst[1:30000, 1:100]),
  n_200_30k   = scale(counts_vst[1:30000, 1:200]),
  n_500_30k   = scale(counts_vst[1:30000, 1:500]),
  n_1000_30k  = scale(counts_vst[1:30000, 1:1000]),
  n_100_50k   = scale(counts_vst[1:50000, 1:100]),
  n_200_50k   = scale(counts_vst[1:50000, 1:200]),
  n_500_50k   = scale(counts_vst[1:50000, 1:500]),
  n_1000_50k  = scale(counts_vst[1:50000, 1:1000]),
  iterations = 1,
  check = FALSE
)


res_scale_tab <- bind_cols(tibble(expresion = names(res_scale$expression)), res_scale[, 2:9])


# log2 transform -----------------------------------------------------------


fpkm_counts <-
  fpkm_counts %>%
  as.data.frame() %>%
  tibble::column_to_rownames("gene_name") %>%
  mutate(across(everything(), as.numeric))



res_log <- bench::mark(
  n_100_5k   = log2(fpkm_counts[1:5000 , 1:100] + 1),
  n_200_5k   = log2(fpkm_counts[1:5000 , 1:200] + 1),
  n_500_5k   = log2(fpkm_counts[1:5000 , 1:500] + 1),
  n_1000_5k  = log2(fpkm_counts[1:5000 , 1:1000] + 1),
  n_100_10k  = log2(fpkm_counts[1:10000, 1:100] + 1),
  n_200_10k  = log2(fpkm_counts[1:10000, 1:200] + 1),
  n_500_10k  = log2(fpkm_counts[1:10000, 1:500] + 1),
  n_1000_10k = log2(fpkm_counts[1:10000, 1:1000] + 1) ,
  n_100_30k  = log2(fpkm_counts[1:30000, 1:100] + 1),
  n_200_30k  = log2(fpkm_counts[1:30000, 1:200] + 1),
  n_500_30k  = log2(fpkm_counts[1:30000, 1:500] + 1),
  n_1000_30k = log2(fpkm_counts[1:30000, 1:1000] + 1),
  n_100_50k  = log2(fpkm_counts[1:50000, 1:100] + 1),
  n_200_50k  = log2(fpkm_counts[1:50000, 1:200] + 1),
  n_500_50k  = log2(fpkm_counts[1:50000, 1:500] + 1),
  n_1000_50k = log2(fpkm_counts[1:50000, 1:1000] + 1),
  iterations = 1,
  check = FALSE
  
)


res_log_tab <- bind_cols(tibble(expresion = names(res_log$expression)), res_log[, 2:9])


# Batch correct TCGA GTex -------------------------------------------------

patient_select <-
  function(x, size, gene_size) {
    a <- sample(x$patients, replace = TRUE, size = size)
    output <- list()
    output[[1]] <- merged_metadata[merged_metadata$patients %in% a, ]
    output[[2]] <- merged[1:gene_size , colnames(merged) %in% a]
    return(output)
  }

batch_corr <- function(x, size, gene_size) {
  obj <- patient_select(x, size, gene_size)
  limma::removeBatchEffect(obj[[2]], batch = obj[[1]]$cohort)
}


batch_corr(merged_metadata, 100, 500)


res_batch <- bench::mark(
  n_100_5k    = batch_corr(merged_metadata, 100, 5000),
  n_200_5k    = batch_corr(merged_metadata, 200, 5000),
  n_500_5k    = batch_corr(merged_metadata, 500, 5000),
  n_1000_5k   = batch_corr(merged_metadata, 1000, 5000),
  n_100_10k    = batch_corr(merged_metadata, 100, 10000),
  n_200_10k    = batch_corr(merged_metadata, 200, 10000),
  n_500_10k    = batch_corr(merged_metadata, 500, 10000),
  n_1000_10k   = batch_corr(merged_metadata, 1000, 10000),
  n_100_30k    = batch_corr(merged_metadata, 100, 30000),
  n_200_30k    = batch_corr(merged_metadata, 200, 30000),
  n_500_30k    = batch_corr(merged_metadata, 500, 30000),
  n_1000_30k   = batch_corr(merged_metadata, 1000, 30000),
  n_100_50k    = batch_corr(merged_metadata, 100, 50000),
  n_200_50k    = batch_corr(merged_metadata, 200, 50000),
  n_500_50k    = batch_corr(merged_metadata, 500, 50000),
  n_1000_50k   = batch_corr(merged_metadata, 1000, 50000),
  iterations = 1,
  check = FALSE
  
)


res_batch_tab <- bind_cols(tibble(expresion = names(res_batch$expression)), res_batch[, 2:9])



table <- bind_rows(
  res_tab,
  res_vst_gtex_tab,
  res_scale_tab,
  res_log_tab,
  res_batch_tab
) 

table$analysis <- 
  c(
  rep(c("VST_TCGA"), each = 16),
rep(c("VST_GTEx"), each = 8),
rep(c("Scale VST", "Log2 transform FPKM", "Batch TCGA-GTEx crudos"), each = 16)
)

table <- table %>% 
  dplyr::relocate(analysis) %>%
  dplyr::select(-c(n_itr, n_gc))



table <- 
  table %>% 
  mutate(
    n_samples = gsub("^n_(.*)_(.*)k$", "\\1", expresion),
    n_genes = (as.numeric(gsub("^n_(.*)_(.*)k$", "\\2", expresion)) * 1000),
    Input = case_when(
      analysis == "VST_TCGA" ~ "Conteos crudos",
      analysis == "VST_GTEx" ~ "Conteos crudos",
      analysis == "Scale VST" ~ "Conteos transformados con VST",
      analysis == "Log2 transform FPKM" ~ "Conteos FPKM",
      analysis == "Batch TCGA-GTEx crudos" ~ "Conteos crudos, metadata con cohorte de pertenencia",
    ),
    Output = case_when(
      analysis == "VST_TCGA" ~ "Conteos transformados con VST",
      analysis == "VST_GTEx" ~ "Conteos transformados con VST",
      analysis == "Scale VST" ~ "Conteos transformados con VST escalado y centrado",
      analysis == "Log2 transform FPKM" ~ "Conteos FPKM log 2",
      analysis == "Batch TCGA-GTEx crudos" ~ "Conteos crudos con corrección de lote",
    )
  ) %>% 
  dplyr::relocate(c(analysis, Input, Output, n_samples, n_genes)) %>% 
  select(-expresion, -min, -median) %>% 
  rename(
    "Analisis" = "analysis",
    "Número de muestras" = "n_samples",
    "Número de genes" = "n_genes",
    "Tiempo total" = "total_time",
    "Iteraciones segundo"  = "itr.sec",
    "Memoria RAM" = "mem_alloc",
    "Tiempo de garbage collection" = "gc.sec")


table %>% 
  flextable::flextable() %>%
  flextable::autofit() %>% 
  theme_zebra() %>% 
  flextable::merge_at(i=c(1:16), j = 1) %>% 
  flextable::merge_at(i=c(17:24), j = 1) %>% 
  flextable::merge_at(i=c(25:40), j = 1) %>% 
  flextable::merge_at(i=c(41:56), j = 1) %>% 
  flextable::merge_at(i=c(57:72), j = 1) %>% 
  flextable::merge_at(i=c(1:16), j = 2) %>% 
  flextable::merge_at(i=c(17:24), j = 2) %>% 
  flextable::merge_at(i=c(25:40), j = 2) %>% 
  flextable::merge_at(i=c(41:56), j = 2) %>% 
  flextable::merge_at(i=c(57:72), j = 2) %>% 
  flextable::merge_at(i=c(1:16), j = 3) %>% 
  flextable::merge_at(i=c(17:24), j = 3) %>% 
  flextable::merge_at(i=c(25:40), j = 3) %>% 
  flextable::merge_at(i=c(41:56), j = 3) %>% 
  flextable::merge_at(i=c(57:72), j = 3) %>% 
  flextable::hline(j = c(4:9)) %>% 
  flextable::hline(i = 16, j = c(1:9) ,border = officer::fp_border(width = 2)) %>% 
  flextable::hline(i = 24, j = c(1:9) ,border = officer::fp_border(width = 2)) %>% 
  flextable::hline(i = 40, j = c(1:9) ,border = officer::fp_border(width = 2)) %>% 
  flextable::hline(i = 56, j = c(1:9) ,border = officer::fp_border(width = 2)) %>% 
  flextable::hline(i = 72, j = c(1:9) ,border = officer::fp_border(width = 2)) %>% 
  
  flextable::padding(padding.top = 1.5, padding.bottom = 1.5, part = "all")


write.csv(table, "./results/time_table.csv")



table <- read.csv("results/time_table.csv", row.names = 1)
