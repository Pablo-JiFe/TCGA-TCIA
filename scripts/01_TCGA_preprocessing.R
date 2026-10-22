library(TCGAbiolinks)
library(dplyr)
library(arrow)
library(ggplot2)
library(SummarizedExperiment)
library(UCSCXenaTools)
library(DESeq2)


# 1.- Loading data --------------------------------------------------------

directory <- ""

# 1.1 Query  obtaining the desiered data 

tcga_rna <- GDCquery("TCGA-BRCA",
                     data.category = "Transcriptome Profiling",
                     access = "open",
                     experimental.strategy = "RNA-Seq",
                     workflow.type = "STAR - Counts"
)

GDCdownload(tcga_rna, method = "api", files.per.chunk = 5, directory = directory)

# 1.2 Prepare data for usage

tcga_brca_data <- GDCprepare(tcga_rna, directory = directory,
                             summarizedExperiment = TRUE) 

counts_raw <- assay(tcga_brca_data, "unstranded")

counts_fpkm <- assay(tcga_brca_data, "fpkm_unstrand")

genes_tcga <- rowData(tcga_brca_data) %>% 
  as.data.frame() %>% 
  dplyr::select(c("gene_id", "gene_name", "gene_type"))

#for ram

rm("tcga_brca_data")


# Data frame

counts_raw <- 
  counts_raw %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("gene_id") %>% 
  left_join(genes_tcga, by = "gene_id") %>% 
  dplyr::select(-c(gene_id, gene_type))

counts_raw$var <- apply(as.matrix(counts_raw[, !colnames(counts_raw) %in% "gene_name"]), 1, var)

counts_raw <- 
  counts_raw %>% 
  group_by(gene_name) %>% 
  slice_max(order_by = var, n = 1, with_ties = FALSE) %>% 
  select(-var) %>% 
  ungroup()

counts_fpkm <- 
  counts_fpkm %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("gene_id") %>% 
  left_join(genes_tcga, by = "gene_id") %>% 
  dplyr::select(-c(gene_id, gene_type))

counts_fpkm$var <- apply(as.matrix(counts_fpkm[, !colnames(counts_fpkm) %in% "gene_name"]), 1, var)

counts_fpkm <- 
  counts_fpkm %>% 
  group_by(gene_name) %>% 
  slice_max(order_by = var, n = 1, with_ties = FALSE) %>% 
  select(-var) %>% 
  ungroup()


# 1.- Metadata ------------------------------------------------------------

# 1.1 Generate and Query

data_query <- XenaGenerate(subset = XenaDatasets == "TCGA.BRCA.sampleMap/BRCA_clinicalMatrix") %>% 
  XenaQuery()

# 1.2 Download 

xe_download <- XenaDownload(data_query, destdir = "./data/")

# 1.3 Prepare (Load) the data

brca_clinical <- XenaPrepare(xe_download)

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



saveRDS(counts_raw, "./data/counts_raw_tcga.RDS")
saveRDS(refined_data, "./data/metadata_tcga.RDS")
saveRDS(genes_tcga, "./data/genes_tcga.RDS")
saveRDS(counts_fpkm, "./data/counts_fpkm.RDS")

rm(list = ls())
gc()
