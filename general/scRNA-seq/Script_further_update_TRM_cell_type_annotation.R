setwd("/home/yuq22/Bacteria_TRM/TRM/TRM_from_coculture/cell_type_annotation/coarse_cell_type_annotation/update_annotation_after_discussion")
library(ggplot2)
library(Seurat)
#trm <- readRDS("/home/yuq22/group_share/Bacteria_TRM/used_object/Dat_TRM_with_updated_cell_type_anno.rds")
#
## update cell type annotation
#ct1 <- c("8:Activated CD4 T", "3:Activated CD8 T", "4:gdT", "5:ILC3", "9:MAIT-like CD4", "1:Memory B", "2:NK cell", "7:Treg", "6:G2M CD4")
#ct2 <- setNames(c("8:CD4 TRM", "3:CD8 TRM", "4:gd T", "5:ILC", "9:CD4 TRM Th17", "1:Memory B", "2:NK cell", "7:Treg", "6:Cycling CD4"), ct1)
#trm$Updated_cell_type_index_v2 <- ct2[as.character(trm$Updated_cell_type_index)] 
#trm$Updated_cell_type_index <- NULL
#trm$Updated_cell_type <- NULL
#trm$Cell_type <- NULL
#trm$Updated_cell_type_index_v2[which(trm$Updated_cell_type_index_v2=="6:Cycling CD4")] <- "9:Cycling CD4"
#trm$Updated_cell_type_index_v2[which(trm$Updated_cell_type_index_v2=="9:CD4 TRM Th17")] <- "6:CD4 TRM Th17"
#trm$Updated_cell_type <- NULL
#saveRDS(trm, file="/home/yuq22/group_share/Bacteria_TRM/used_object/Dat_TRM_with_updated_cell_type_anno_after_discussion.rds")


marker_list <- list(
    "Memory B"=c("CD79B", "LINC01781", "TNFRSF13B"),
    "NK"=c("GZMB", "GZMK", "NKG7","EOMES"),
    "CD8 TRM"=c("XCL1", "CD8A", "CD8B"),
    "gdT"=c("CD160", "KLRC2", "KLRD1"),
    "ILC"=c("KIT", "ZFP36L1"),
    "MAIT/Th17"=c("SLC4A10", "KLRB1","RORC","IL23R","IL26","IL17A"), 
    "Treg"=c("LTB", "ARID5B", "CCR7", "CTLA4", "FOXP3", "TWIST1"),
    "CD4 TRM"=c("RGCC", "GPR183", "INPP4B"),
    "G2M"=c("CD4", "MKI67")
    
)

final_marker_list <- list(
    "Memory B"=c("CD79B", "TNFRSF13B"),
    "NK"=c("GZMB", "GZMK", "NKG7","EOMES"),
    "CD8 TRM"=c("CD8A", "CD8B"),
    "gdT"=c("CD160", "KLRC2", "KLRD1"),
    "ILC"=c("KIT"),
    "Th17"=c("RORC","IL23R","IL26","IL17A","IL22"), 
    "Treg"=c("CTLA4", "FOXP3", "TWIST1"),
    "CD4 TRM"=c("INPP4B"),
    "G2M"=c("CD4", "MKI67")
    
)

# load in vitro TRM data
trm <- readRDS("/home/yuq22/group_share/Bacteria_TRM/used_object/Dat_TRM_with_updated_cell_type_anno_after_discussion.rds")
p2 <- SCpubr::do_DotPlot(sample = trm, 
             features = final_marker_list,
             scale=T,
             group.by="Updated_cell_type_index_v2")
p2
ggsave(p2, filename="Plot_dotplot_adult_human_inVitro_immune_cell_type_markers.pdf", width=10, height=6)


# subset tissue immune data to those matched with in invitro data
tissue <- readRDS("/home/yuq22/Bacteria_TRM/intestine_multi_region_scRNA-seq_atlas/immune_cells/Res_gut_cell_atlas_adult_human_immune_no_batch_correction.rds")
trm$Source <- "inVitro"
ct1 <- c("Activated CD4 T", "Activated CD8 T", "gdT", "ILC3", "Th17", "Memory B", "NK cell", "Treg") 
tissue_subset <- subset(tissue, cells=colnames(tissue)[which(tissue$Cell_type%in%ct1)])
tissue_subset$Source <- "tissue"
ct2 <- setNames(c("8:CD4 TRM", "3:CD8 TRM", "4:gd T", "5:ILC", "6:CD4 TRM Th17", "1:Memory B", "2:NK cell", "7:Treg"), ct1)
tissue_subset$Updated_cell_type_index <- ct2[as.character(tissue_subset$Cell_type)]
saveRDS(tissue_subset, file="Dat_gut_cell_atlas_adult_human_immune_no_batch_correction_subset.rds")

p2 <- SCpubr::do_DotPlot(sample = tissue_subset, 
             features = final_marker_list,
             scale=T,
             group.by="Updated_cell_type_index")
p2
ggsave(p2, filename="Plot_dotplot_adult_human_tissue_immune_cell_type_markers.pdf", width=10, height=6)

trm$Updated_cell_type_index <- trm$Updated_cell_type_index_v2
combined <- merge(x=trm, y=tissue_subset)
combined$Cell_type_per_source <- paste(combined$Updated_cell_type_index, combined$Source, sep=":")
combined$Cell_type_per_source <- factor(combined$Cell_type_per_source, levels=sort(unique(combined$Cell_type_per_source)))
Idents(combined) <- combined$Cell_type_per_source
saveRDS(combined, file="Dat_combined_immune_inVitro_and_tissue_data_v2.rds")
saveRDS(combined, file="/home/yuq22/group_share/Bacteria_TRM/used_object/Dat_combined_immune_inVitro_and_tissue_data_v2.rds")

combined <- readRDS("/home/yuq22/group_share/Bacteria_TRM/used_object/Dat_combined_immune_inVitro_and_tissue_data_v2.rds")
p2 <- SCpubr::do_DotPlot(sample = combined, 
             features = final_marker_list,
             scale=T)
p2

ggsave(p2, filename="Plot_dotplot_adult_human_tissue_and_inVitro_immune_cell_type_markers_with_IL22.pdf", width=8, height=6)


