library(dplyr)

# GTEx loading ------------------------------------------------------------

gtex_ds <- open_dataset("D:/GTEx/GTEx_Analysis_2026-05-19_v11_RNASeQCv2.4.3_gene_reads.parquet")

metadata_gte <- readxl::read_xlsx("D:/GTEx/GTEx_Analysis_v11_Annotations_SampleAttributesDD.xlsx")

metadata_gte2 <- readxl::read_xlsx("D:/GTEx/GTEx_Analysis_v11_Annotations_SubjectPhenotypesDD.xlsx")

metadata <- read.delim("D:/GTEx/metadata_gte.txt")

metadata2 <- read.delim("D:/GTEx/metadata_2.txt")


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


saveRDS(counts_df, "data/counts_df_gte.RDS")
saveRDS(metadata_breast_woman, "data/metadata_breast_gte.RDS")


rm(list = ls())
gc()
