library(dplyr)

# GTEx loading ------------------------------------------------------------

gtex_ds <- open_dataset("D:/GTEx/GTEx_Analysis_2026-05-19_v11_RNASeQCv2.4.3_gene_reads.parquet")

metadata_gte <- readxl::read_xlsx("D:/GTEx/GTEx_Analysis_v11_Annotations_SampleAttributesDD.xlsx")

metadata_gte2 <- readxl::read_xlsx("D:/GTEx/GTEx_Analysis_v11_Annotations_SubjectPhenotypesDD.xlsx")

metadata <- read.delim("D:/GTEx/metadata_gte.txt")

metadata2 <- read.delim("D:/GTEx/metadata_2.txt")

genes_tcga <- readRDS("data/genes_tcga.RDS")



metadata_women <- metadata2 %>% 
  filter(SEX == 2)

metadata_breast <- 
  metadata %>% 
  filter(SMTSD == "Breast - Mammary Tissue") %>% 
  mutate(sample = gsub("^([^-]+-[^-]+).*", "\\1", SAMPID))


metadata_breast_woman <- metadata_breast[metadata_breast$sample %in% metadata_women$SUBJID, ] %>% 
  dplyr::rename("SUBJID" = "sample") %>% 
  left_join(metadata_women, by = "SUBJID")


selected_cols <- c(intersect(colnames(gtex_ds), metadata_breast_woman$SAMPID), "Description")


counts_df <- read_parquet(
  "D:/GTEx/GTEx_Analysis_2026-05-19_v11_RNASeQCv2.4.3_gene_reads.parquet",
  col_select = selected_cols
)




# Converting names to symbol ----------------------------------------------

genes_tcga <- genes_tcga %>%
  dplyr::select(gene_id,
                gene_name,
                gene_type ) %>%
  filter(!(gene_type == ""))


genes_tcga$gene_id <- gsub("\\..*", "" , genes_tcga$gene_id)

common_ensembl <- data.frame(gene_id = counts_df$Description[counts_df$Description %in% genes_tcga$gene_id])

common_ensembl_sym <- 
  common_ensembl %>% 
  left_join(genes_tcga,
            by = "gene_id")

counts_df$new_symbol <- common_ensembl_sym$gene_name[base::match(counts_df$Description, common_ensembl_sym$gene_id)]

counts_df <- 
  counts_df %>% 
  mutate(old_symbol = 
           ifelse(
             Description %in% genes_tcga$gene_name,
             Description,
             NA)
  )

counts_df <- 
  counts_df %>% 
  mutate(Symbol = coalesce(new_symbol, old_symbol)) %>% 
  dplyr::select(-c(
    new_symbol, 
    old_symbol
  ))



counts_df$variance <- apply(counts_df[, !colnames(counts_df) %in% c("Description", "Symbol")], MARGIN = 1, FUN = function(x){
  var(x, na.rm = TRUE)
})


counts_df <- 
  counts_df %>% 
  group_by(Description) %>% 
  slice_max(order_by = variance, n = 1, with_ties = FALSE) %>% 
  tibble::column_to_rownames("Description") %>% 
  ungroup() %>% 
  group_by(Symbol) %>% 
  slice_max(order_by = variance, n = 1, with_ties = FALSE) %>% 
  dplyr::select(-variance)



saveRDS(counts_df, "data/counts_df_gte.RDS")
saveRDS(metadata_breast_woman, "data/metadata_breast_gte.RDS")


rm(list = ls())
gc()
