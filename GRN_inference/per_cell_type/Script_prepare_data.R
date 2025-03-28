# import packages
library(Seurat)
library(Signac)
library(tidyverse)
library(Pando)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(doParallel)
registerDoParallel(10)
library(dplyr)
library(ggplot2)

# set working folder
setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/GRN_on_more_filtered_data/GRN_per_cell_type")
# load data
combined <- readRDS("../Res_chip_infection_multiome.rds")
rna <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/more_filterings/Dat_salmonella_colon_epithelium_scRNA-seq_data_with_subtypes.rds")
combined$cluster_group <- rna@meta.data[colnames(combined), "cluster_group"]
combined$cluster_group_per_condition <- rna@meta.data[colnames(combined), "cluster_group_per_condition"]
saveRDS(combined, file="/home/yuq22/Bacteria_TRM/used_object/Res_chip_infection_multiome.rds")
## subset to individual cell types
selected_cell_types <- c("Colonocyte_4", "Colonocyte_2", "Colonocyte_1", "Colonocyte_3", "Goblet", "Colon_SC")
for(x in selected_cell_types){
  subset_data <- subset(combined, subset = cluster_group ==x)
  saveRDS(subset_data, file=paste0("Res_chip_infection_multiome_", x, ".rds"))
}
