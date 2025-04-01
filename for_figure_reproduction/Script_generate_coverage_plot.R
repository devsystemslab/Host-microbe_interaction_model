setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/C16/coverage_plot")

library(Seurat)
library(Signac)
library(dplyr)

# load multiome data
combined_rna <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/ATAC_on_more_filtered_data/motif_analysis/Res_salmonella_scMultiome_with_chromVAR_and_motif_res.rds")
all <- combined_rna
size <- table(combined_rna$RNA_subtype_group_per_condition)
combined_subset <- subset(combined_rna, subset = RNA_subtype_group_per_condition %in% names(size)[which(size>20)] & Coarse_cell_type != "TA")
combined_rna <- combined_subset
saveRDS(combined_rna, file="Res_salmonella_scMultiome_with_chromVAR_and_motif_res_subset.rds")
# get regions next to the feature genes
peak_anno <- readRDS('/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/peak_annotation/Data_frame_peak_gene_pairs.rds')

# load DARs
combined_dar_res <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/salmonella/multiome_analysis/ATAC_on_more_filtered_data/Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_combined_control_cells_DAR_all_res_on_RNA_cell_type_level.rds")
# focused on Salmonella infection up-regulated regions in colonocytes and stem cells
sal_up_region <- combined_dar_res[grepl("Salmonella", combined_dar_res$group) & grepl("SC|colonocyte", combined_dar_res$group),] %>% filter(logFC>0) %>% group_by(group) %>% top_n(200, wt=logFC)
union_dar <- unique(sal_up_region$feature)

## load DEGs induced upon infection
combined_deg_res <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/more_filterings/Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_matched_control_cells_and_cocluster_control_cells_DEG_sig_res.rds")
sal_deg <- combined_deg_res %>% group_by(group) %>% top_n(20, wt=logFC)
union_deg <- unique(sal_deg$feature)

# get regions next to feature genes, defined as DARs and next to DEGs
# generate coverage plot
selected_genes <- c("STAT3", "SOCS3")
selected_genes <- c("DMBT1", "PLA2G2A")
selected_genes <- c("TRIM40", "LURAP1L") # direct target of STAT3
setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/C16/coverage_plot/STAT3_direct_target")
selected_genes <- c("HSPA5",  "TXN",      "SLC7A11",  "CDC42BPA") # 2nd direct target of STAT3
setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/C16/coverage_plot/STAT3_2nd_direct_target")
idx <- which(peak_anno$external_gene_name %in% selected_genes)
peaks <- peak_anno$name[idx]
length(peaks)


# generate coverage plot
library(Seurat)
library(Signac)
library(IRanges)
library(GenomicRanges)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
# load evolutionary signature list
evo_list <- readRDS("/projects/site/pred/ihb-intestine-evo/evo_signature/summary_Feb_2024/Dat_SNC_GWAS_HAR_PS_HAQER.rds")

colors <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/used_object/cols.rds")
group <- sort(unique(combined_rna$RNA_subtype_group_per_condition))
mat <- do.call('rbind', strsplit(group, split=":"))
ct <- mat[,1]
names(ct) <- group
group_cols <- setNames(colors$cell_type[ct], group)
source("/home/yuq22/ihb-intestine-evo/common_script/plotting/Script_plotting_functions.R")
p_list <- coverage_plot_with_evo_sig(seu_obj=combined_rna, 
                                     peak_assay="peaks_merged", 
                                     rna_assay="RNA", 
                                     group="RNA_subtype_group_per_condition", 
                                     colors.use=group_cols, 
                                     peaks=peaks, 
                                     evo_list_path="/projects/site/pred/ihb-intestine-evo/evo_signature/summary_Feb_2024/Dat_SNC_GWAS_HAR_PS_HAQER.rds", 
                                     peak_anno_path="/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/peak_annotation/Data_frame_peak_gene_pairs.rds", 
                                     peak_anno_peak_col="name",
                                     peak_anno_gene_col="external_gene_name",
                                     extend_up=10000,
                                     extend_down=10000,
                                     window_size=500,
                                     color="#303030",
                                     gene_mode="ChIPseeker", 
                                     g_vec=NULL,
                                     do_plot_combined=FALSE,
                                     do_plot_individual=TRUE,
                                     combined_plot_name=NULL,
                                     plot_suffix="-3")
