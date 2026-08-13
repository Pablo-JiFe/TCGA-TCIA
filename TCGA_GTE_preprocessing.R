library(TCGAbiolinks)
library(dplyr)
library(arrow)
library(ggplot2)


# 1.- Loading data --------------------------------------------------------



# 1.1 Query  obtaining the desiered data 

tcga_rna <- GDCquery("TCGA-BRCA",
                     data.category = "Transcriptome Profiling",
                     access = "open",
                     experimental.strategy = "RNA-Seq",
                     workflow.type = "STAR - Counts"
)

#   GDCdownload(tcga_rna, method = "api", files.per.chunk = 5)

# 1.2 Prepare data for usage

tcga_brca_data <- GDCprepare(tcga_rna, directory = "D:/tcga/GDCdata",
                             summarizedExperiment = FALSE) # non SE to keep gene_name, gene_id and gene_type


counts_col <- grep(
  pattern = "^fpkm_unstranded_", # Normalized
  x = names(tcga_brca_data),
  value = TRUE)

# Keep all counts and gene ids

counts_fpkm <- tcga_brca_data %>% 
  dplyr::select(all_of(counts_col),
                gene_name,
                gene_id,
                gene_type) 
# Object with only gene ids

genes_tcga <- tcga_brca_data %>% 
  dplyr::select(c(gene_id,
                  gene_name,
                  gene_type))

#for ram

rm("tcga_brca_data")

# Eliminate ensembl version

genes_tcga$gene_id <- gsub("\\..*", "" , genes_tcga$gene_id)

# Eliminate fpkm_unstranded text form titles

colnames(counts_fpkm) <- gsub("^fpkm_unstranded_", "", colnames(counts_fpkm))

# Data frame

counts_fpkm <- 
  counts_fpkm %>% 
  as.data.frame()

# Keep the replicates identified as normal tissue as seen in https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes

counts_control <- counts_fpkm[, substr(colnames(counts_fpkm), 14, 15) == "11" | colnames(counts_fpkm) %in% c("gene_name", "gene_id")] # 11 for normal tissue

counts_disease <- counts_fpkm[, substr(colnames(counts_fpkm), 14, 15) == "01" | colnames(counts_fpkm) %in% c("gene_name", "gene_id")] # 01 for primary disease

# GTEx loading ------------------------------------------------------------

# gte_counts <- read_parquet("D:/GTEx/GTEx_Analysis_2026-05-19_v11_RNASeQCv2.4.3_gene_reads.parquet")



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


metadata_breast_woman <- metadata_breast[metadata_breast$sample %in% metadata_women$SUBJID, ]


selected_cols <- c(intersect(colnames(gtex_ds), metadata_breast_woman$SAMPID), "Description")

counts_df <- gtex_ds %>%
  dplyr::select(all_of(selected_cols)) %>%
  collect()

counts_df <- read_parquet(
  "D:/GTEx/GTEx_Analysis_2026-05-19_v11_RNASeQCv2.4.3_gene_reads.parquet",
  col_select = selected_cols
)


saveRDS(counts_fpkm, "data/counts_fpkm.RDS")
saveRDS(counts_control, "data/counts_control.RDS")
saveRDS(counts_disease, "data/counts_disease.RDS")

saveRDS(counts_df, "data/counts_df.RDS")
saveRDS(metadata_breast_woman, "data/metadata_breast_woman.RDS")

saveRDS(genes_tcga, "data/genes_tcga.RDS")

rm(list = ls())
