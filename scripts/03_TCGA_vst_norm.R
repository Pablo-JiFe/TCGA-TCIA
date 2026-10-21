library(dplyr)
library(ggplot2)
library(bench)
library(arrow)

counts_raw <- readRDS("data/counts_raw_tcga.RDS")
refined_data <- readRDS("./data/metadata_tcga.RDS")
genes_tcga <- readRDS("./data/genes_tcga.RDS")

plate_center <- data.frame(
  row.names = colnames(counts_raw),
  plate = gsub(
    pattern = "^TCGA-(.*)-(.*)-(.*)-(.*)-(.*)-(.*)$",
    replacement = "\\5",
    x = colnames(counts_raw)
  ),
  center = gsub(
    pattern = "^TCGA-(.*)-(.*)-(.*)-(.*)-(.*)-(.*)$",
    replacement = "\\6",
    x = colnames(counts_raw)
  )
)

# VST transform -----------------------------------------------------------


res <- bench::mark(
  n_100_5k  = DESeq2::vst(counts_raw[1:5000, 1:100]),
  n_500_5k  = DESeq2::vst(counts_raw[1:5000, 1:200]),
  n_100_5k  = DESeq2::vst(counts_raw[1:5000, 1:500]),
  n_500_5k  = DESeq2::vst(counts_raw[1:5000, 1:1000]),
  n_100_10k  = DESeq2::vst(counts_raw[1:10000, 1:100]),
  n_500_10k  = DESeq2::vst(counts_raw[1:10000, 1:200]),
  n_100_10k  = DESeq2::vst(counts_raw[1:10000, 1:500]),
  n_500_10k  = DESeq2::vst(counts_raw[1:10000, 1:1000]),
  n_100_30k  = DESeq2::vst(counts_raw[1:30000, 1:100]),
  n_500_30k  = DESeq2::vst(counts_raw[1:30000, 1:200]),
  n_100_30k  = DESeq2::vst(counts_raw[1:30000, 1:500]),
  n_500_30k  = DESeq2::vst(counts_raw[1:30000, 1:1000]),
  n_100_60k  = DESeq2::vst(counts_raw[1:60000, 1:100]),
  n_500_60k  = DESeq2::vst(counts_raw[1:60000, 1:200]),
  n_100_60k  = DESeq2::vst(counts_raw[1:60000, 1:500]),
  n_500_60k  = DESeq2::vst(counts_raw[1:60000, 1:1000]),
  iterations = 1,
  check = FALSE
)


cbind(data.frame(expresion = names(res$expression)), res[,2:9]) %>% 
  flextable::flextable() %>% 
  flextable::autofit()

counts_vst <- DESeq2::vst(object = counts_raw)



# PCA before correction ---------------------------------------------------

pca <- prcomp(t(counts_vst))

pca$x %>% 
  ggplot(aes(x = PC1, y = PC2, colour = plate_center$plate, shape = plate_center$center)) + 
  geom_point()


# Batch correction --------------------------------------------------------

limma::removeBatchEffect()



# Scale ---------------------------------------------------------------


res_scale <- bench::mark(
  n_100_5k  = scale(counts_vst[1:5000, 1:100]),
  n_500_5k  = scale(counts_vst[1:5000, 1:200]),
  n_100_5k  = scale(counts_vst[1:5000, 1:500]),
  n_500_5k  = scale(counts_vst[1:5000, 1:1000]),
  n_100_10k  = scale(counts_vst[1:10000, 1:100]),
  n_500_10k  = scale(counts_vst[1:10000, 1:200]),
  n_100_10k  = scale(counts_vst[1:10000, 1:500]),
  n_500_10k  = scale(counts_vst[1:10000, 1:1000]),
  n_100_30k  = scale(counts_vst[1:30000, 1:100]),
  n_500_30k  = scale(counts_vst[1:30000, 1:200]),
  n_100_30k  = scale(counts_vst[1:30000, 1:500]),
  n_500_30k  = scale(counts_vst[1:30000, 1:1000]),
  n_100_60k  = scale(counts_vst[1:60000, 1:100]),
  n_500_60k  = scale(counts_vst[1:60000, 1:200]),
  n_100_60k  = scale(counts_vst[1:60000, 1:500]),
  n_500_60k  = scale(counts_vst[1:60000, 1:1000]),
  iterations = 1,
  check = FALSE
)


res_scale

vst_scaled <- scale(counts_vst)



# Name genes and delete duplicate -----------------------------------------

counts_vst_var <- counts_vst %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("gene_id") %>% 
  left_join(genes_tcga, by = "gene_id") %>% 
  dplyr::select(-c(gene_id, gene_type))

counts_vst_var$var <- apply(as.matrix(counts_vst_var[, !colnames(counts_vst_var) %in% "gene_name"]), 1, var)

counts_vst_un <- 
  counts_vst_var %>% 
  group_by(gene_name) %>% 
  slice_max(order_by = var, n = 1, with_ties = FALSE) %>% 
  select(-var) %>% 
  ungroup()


genes_tok <- counts_vst_var %>% 
  select(gene_name, var) %>% 
  group_by(gene_name) %>% 
  slice_max(order_by = var, n = 1, with_ties = FALSE) %>% 
  arrange(desc(var)) 

top_5k_genes <- genes_tok$gene_name[1:5000]

# Keep the replicates identified as normal tissue as seen in https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes

counts_control <- counts_vst_un[, substr(colnames(counts_vst_un), 14, 15) == "11" | colnames(counts_vst_un) %in% c("gene_name", "gene_id")] # 11 for normal tissue

counts_disease <- counts_vst_un[, substr(colnames(counts_vst_un), 14, 15) == "01" | colnames(counts_vst_un) %in% c("gene_name", "gene_id")] # 01 for primary disease

counts_mets <- counts_vst_un[, substr(colnames(counts_vst_un), 14, 15) == "06" | colnames(counts_vst_un) %in% c("gene_name", "gene_id")] # 06 for metastasis

# Keep the biological replicate for each sample with highest variance

counts_disease_t <- 
  counts_disease[, !colnames(counts_disease) %in% "gene_name"] %>% 
  t() 

variance <- apply(counts_disease_t, 1, var)

deduplicate_names <- 
  variance %>% 
  as.data.frame() %>% 
  dplyr::rename("var" = ".") %>% 
  tibble::rownames_to_column("id") %>% 
  mutate(id_short = substr(start = 1, stop = 12, x = id)) %>% 
  group_by(id_short) %>% 
  slice_max(order_by = var, n = 1, with_ties = FALSE) %>% 
  tibble::column_to_rownames("id_short") 


counts_disease_final <- 
  counts_disease[, c(deduplicate_names$id, "gene_name")] 


colnames(counts_control) <- substr(start = 1, stop = 12, x = colnames(counts_control))

colnames(counts_mets) <- substr(start = 1, stop = 12, x = colnames(counts_mets))



write_parquet(counts_disease_final[counts_disease_final$gene_name %in% top_5k_genes, ], "data/counts_disease5k_tcga.parquet")
write_parquet(counts_control, "data/counts_control_tcga.parquet")
write_parquet(counts_disease_final, "data/counts_disease_tcga.parquet")
write_parquet(counts_mets, "data/counts_mets_tcga.parquet")
write_parquet(counts_vst_un, "./data/counts_vst_un.parquet")


