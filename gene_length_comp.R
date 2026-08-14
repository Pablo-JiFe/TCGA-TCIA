library(dplyr)
library(ggplot2)

counts_fpkm <- readRDS("data/counts_fpkm.RDS")
counts_control <- readRDS("data/counts_control.RDS")
counts_disease <- readRDS("data/counts_disease.RDS")
counts_df <- readRDS("data/counts_df.RDS")
metadata_breast_woman <- readRDS("data/metadata_breast_woman.RDS")
genes_tcga <- readRDS("data/genes_tcga.RDS")


genes_tcga <- counts_fpkm %>%
  dplyr::select(gene_id,
                gene_name,
                gene_type ) %>%
  filter(!(gene_type == ""))


genes_tcga$gene_id <- gsub("\\..*", "" , genes_tcga$gene_id)



# GTE gene names

gene_names_gte <- counts_df$Description # It mixes hugo symbols fore some genes and ENSEMBL for other

# Thus the indexing includes searching for the symbol and ensembl counterparts

length(genes_tcga$gene_name[((genes_tcga$gene_name %in% gene_names_gte) | (genes_tcga$gene_id %in% gene_names_gte))]) # Number of genes in TCGA that can also be found on gte

length(genes_tcga$gene_name[!((genes_tcga$gene_name %in% gene_names_gte) | (genes_tcga$gene_id %in% gene_names_gte))]) # Number of genes in TCGA that can are missing on gte

# Keep the names of the missing genes

remaining_names <- genes_tcga$gene_name[!((genes_tcga$gene_name %in% gene_names_gte) | (genes_tcga$gene_id %in% gene_names_gte))]

# Keep ids of missing genes

remaining_id <- genes_tcga$gene_id[!((genes_tcga$gene_name %in% gene_names_gte) | (genes_tcga$gene_id %in% gene_names_gte))]

# Keep counts in tcga of genes missing in gte

pre_counts_remain <- counts_fpkm[counts_fpkm$gene_name %in% remaining_names, ] 

# Delete the columns of ensembl and gene type 

counts_remain <- pre_counts_remain %>%
  as.data.frame()  %>% 
  select(-c(gene_id, gene_type))

# Variancw

counts_remain$var <- apply(counts_remain[!colnames(counts_remain) == "gene_name"], 1, function(i){var(i, na.rm = TRUE)})

# Plot

ggplot(counts_remain, aes(x = var + 1)) +
  geom_density() +  
  scale_x_log10()

# Summary of variance

summary(counts_remain$var)



counts_remain <- counts_remain %>% 
  tidyr::drop_na(var)

# Gene names with highest variance

counts_remain$gene_name[counts_remain$var == max(counts_remain$var)]

# Gene names of highest variance genes

hig_var_gen <- counts_remain %>% 
  filter(var >= quantile(counts_remain$var, 0.9995)) %>% 
  arrange(desc(var)) %>% 
  select(gene_name)


# Genes in only control and only diseased patients with highest variance

counts_cont_filt <- counts_control[counts_control$gene_name %in% remaining_names,] %>% 
  as.data.frame()

counts_dis_filt <- counts_disease[counts_disease$gene_name %in% remaining_names,] %>% 
  as.data.frame()

# Once again variance

counts_cont_filt$var <- apply(counts_cont_filt[, setdiff(colnames(counts_cont_filt), c("gene_id", "gene_name", "gene_type"))], 1, function(i){var(i, na.rm = TRUE)}) 

# Summary

summary(counts_cont_filt$var)

# Plot

ggplot(counts_cont_filt, aes(x = var + 1)) +
  geom_density() +  
  scale_x_log10()

hig_var_gen_cont <- counts_cont_filt %>% 
  tidyr::drop_na(var) %>% 
  filter(var >= quantile(counts_cont_filt$var, 0.9995, na.rm = TRUE)) %>% 
  arrange(desc(var)) %>% 
  select(gene_name)


# Same

counts_dis_filt$var <- apply(counts_dis_filt[, setdiff(colnames(counts_dis_filt), c("gene_id", "gene_name", "gene_type"))], 1, function(i){var(i, na.rm = TRUE)}) 

summary(counts_dis_filt$var)

ggplot(counts_dis_filt, aes(x = var + 1)) +
  geom_density() +  
  scale_x_log10()

hig_var_gen_dis <- counts_dis_filt %>% 
  tidyr::drop_na(var) %>%  
  filter(var >= quantile(counts_dis_filt$var, 0.9995, na.rm = TRUE)) %>% 
  arrange(desc(var)) %>% 
  select(gene_name)

library(patchwork)


transposed_unique <- counts_dis_filt %>%
  tibble::rownames_to_column("del") %>% 
  select(- c(del,
             gene_id)) %>% 
  group_by(gene_name) %>% 
  slice_max(n = 1, order_by = var, with_ties = FALSE) %>% 
  ungroup() %>% 
  select(-var) %>% 
  tibble::column_to_rownames("gene_name") %>% 
  t() %>% 
  as.data.frame() 

transposed_unique_cont <- counts_cont_filt %>%
  tibble::rownames_to_column("del") %>% 
  select(- c(del,
             gene_id)) %>% 
  group_by(gene_name) %>% 
  slice_max(n = 1, order_by = var, with_ties = FALSE) %>% 
  ungroup() %>% 
  select(-var) %>% 
  tibble::column_to_rownames("gene_name") %>% 
  t() %>% 
  as.data.frame() 


genes_to_analyze <- unique(rbind(hig_var_gen, hig_var_gen_dis, hig_var_gen_cont))$gene_name


for (i in genes_to_analyze) {
  
  
  plot <- transposed_unique %>% ggplot(aes(x = .data[[i]] + 1)) +
    geom_boxplot() +  
    scale_x_log10() +
    labs(title = "Log10 scaled cancer samples")
  
  plot2 <- transposed_unique %>% ggplot(aes(x = .data[[i]])) +
    geom_boxplot() +
    labs(title = "Not scaled cancer samples")
  
  plot3 <- transposed_unique_cont %>% ggplot(aes(x = .data[[i]] + 1)) +
    geom_boxplot() +  
    scale_x_log10() +
    labs(title = "Log10 scaled control samples")
  
  plot4 <- transposed_unique_cont %>% ggplot(aes(x = .data[[i]])) +
    geom_boxplot() +
    labs(title = "Not scaled control samples")
  
  
  join_dis_plot <- plot/plot2
  
  join_cont_plot <- plot3/plot4
  
  cat("#############################\n", i, "\n------------------------\n")
  
  combined_plot <- join_dis_plot | join_cont_plot
  
  top_dis_counts <- head(sort(transposed_unique[[i]], decreasing = TRUE), n = 20)
  
  plot_cuantity_dis <- 
    top_dis_counts %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("index") %>% 
    dplyr::mutate(index = as.numeric(index)) %>%
    ggplot(aes(x = index, y = .)) + 
    geom_point() + 
    scale_x_continuous(breaks = 1:20) +
    labs(title = paste("Counts for the highest count disease samples for ", i))
  
  top_cont_counts <- head(sort(transposed_unique_cont[[i]], decreasing = TRUE), n = 20)
  
  plot_cuantity_cont <- 
    top_cont_counts %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("index") %>% 
    dplyr::mutate(index = as.numeric(index)) %>%
    ggplot(aes(x = index, y = .)) + 
    geom_point() + 
    scale_x_continuous(breaks = 1:20) +
    labs(title = "Counts for the highest count control samples")
  
  
  combined_plot_2 <- plot_cuantity_dis + plot_cuantity_cont
  
  print(combined_plot / combined_plot_2)  
  
  cat("Cancer samples\n")
  print(rownames(transposed_unique)[transposed_unique[[i]] == max(transposed_unique[[i]])])
  print(rownames(transposed_unique)[transposed_unique[[i]] == top_dis_counts[2]])
  cat("Non cancer samples\n")
  print(rownames(transposed_unique_cont)[transposed_unique_cont[[i]] == max(transposed_unique_cont[[i]])])
  print(rownames(transposed_unique_cont)[transposed_unique_cont[[i]] == top_cont_counts[2]])
  
  
  
}




