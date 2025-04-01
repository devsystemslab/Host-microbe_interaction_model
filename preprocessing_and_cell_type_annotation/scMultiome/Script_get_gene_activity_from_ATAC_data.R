setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/ATAC_on_more_filtered_data/get_gene_activity")
library(Seurat)
library(Signac)
library(ggplot2)
library(cowplot)

source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")

# load multiome data
combined_rna <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/ATAC_on_more_filtered_data/motif_analysis/Res_salmonella_scMultiome_with_chromVAR_and_motif_res.rds")
gene.activities <- GeneActivity(combined_rna)
saveRDS(gene.activities, file="Data_gene_activity_from_ATAC.rds")
# add the gene activity matrix to the Seurat object as a new assay and normalize it
combined_rna[['ATAC_gene_activity']] <- CreateAssayObject(counts = gene.activities)
combined_rna <- NormalizeData(
  object = combined_rna,
  assay = 'ATAC_gene_activity',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined_rna$nCount_RNA)
)
combined_rna$Coarse_cell_type_per_condition <- paste(combined_rna$Coarse_cell_type, combined_rna$condition, sep=":")
saveRDS(combined_rna, file="/home/yuq22/Bacteria_TRM/used_object/Res_salmonella_scMultiome_with_gene_activity.rds")

# get cell type per condition gene activity
ave_expr <- getAveExpr(seu.obj=combined_rna, feature.to.calc="Coarse_cell_type_per_condition", size.cutoff=20, colname.prefix=NULL, assay.type="ATAC_gene_activity")
saveRDS(ave_expr, file="Data_ATAC_coarse_cell_type_per_condition_gene_activity.rds")

# plot gene activity of selected genes
heatmap_input <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/NFKB_barrier_innate_sensing/combined_with_Ruben_gene_set/Res_Ruben_genes_heatmap_input.rds")
genes <- rownames(heatmap_input$heatmap_input)

col_order <- c("Colon_SC:Salmonella", "Colon_SC:Control",  "Colonocyte:Salmonella",  "Colonocyte:Control",  "Goblet:Salmonella" , "Goblet:Control")
expr_mat <- t(apply(ave_expr[genes, col_order], 1, function(vec){
    (vec-min(vec))/(max(vec)-min(vec))
}))

pdf("Plot_heatmap_ruben_gene_ATAC_based_gene_activity.pdf", height=10)
gplots::heatmap.2(expr_mat, col = beach.col.heatmap,  trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 0.3, cexCol=1, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
dev.off()

# subset to genes which show correlated gene activity and RNA expression
atac_data <- expr_mat
rna_data <- heatmap_input$heatmap_input

cor_vec <- sapply(seq(nrow(atac_data)), function(i){
    cor(atac_data[i,], rna_data[i,])
})
names(cor_vec) <- rownames(atac_data)

idx <- names(cor_vec)[which(cor_vec > 0.5)]

pdf("Plot_heatmap_ruben_gene_ATAC_based_gene_activity_subset.pdf", height=10)
gplots::heatmap.2(atac_data[idx,], col = beach.col.heatmap,  trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 0.3, cexCol=1, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
dev.off()

pdf("Plot_heatmap_ruben_gene_expr_subset.pdf", height=10)
gplots::heatmap.2(rna_data[idx,], col = beach.col.heatmap,  trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 0.3, cexCol=1, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
dev.off()

# get nearby or correlated regions of selected genes
# load region annotation
peak_anno <- readRDS('/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/peak_annotation/Data_frame_peak_gene_pairs.rds')
nearby_region_df <- peak_anno[which(peak_anno$external_gene_name %in% genes),]
# get coarse_cell_type_per_condition average accessibility data
ave_access <- getAveExpr(seu.obj=combined_rna, feature.to.calc="Coarse_cell_type_per_condition", size.cutoff=20, colname.prefix=NULL, assay.type="peaks_merged")
ave_access_sub <- t(apply(ave_access[, col_order], 1, function(vec){
    (vec-min(vec))/(max(vec)-min(vec))
}))

selected_genes <- intersect(genes, nearby_region_df$external_gene_name)
selected_peaks <- sapply(selected_genes, function(g){
    peak_idx <- which(nearby_region_df$external_gene_name == g)
    peaks <- nearby_region_df$name[peak_idx]
    if(length(peaks)==1){
        return(peaks)
    }else{
        cor_vec <- cor(t(ave_access_sub[peaks,]), rna_data[g,])[,1]
        return(names(cor_vec)[order(cor_vec, decreasing=T)[1]])
    }
    
})

cor_vec <- sapply(seq(nrow(nearby_region_df)), function(i){
    p <- nearby_region_df$name[i]
    g <- nearby_region_df$external_gene_name[i]
    cor(ave_access_sub[p,], rna_data[g,])
})
names(cor_vec) <- nearby_region_df$name

selected_cor_peaks <- selected_peaks[which(cor_vec[selected_peaks]>0.5)]
selected_cor_genes <- nearby_region_df$external_gene_name[which(nearby_region_df$name %in% selected_cor_peaks)]

rna_input <- heatmap_input$heatmap_input[which(rownames(heatmap_input$heatmap_input) %in% selected_cor_genes),]
df <- nearby_region_df[which(nearby_region_df$name %in% selected_cor_peaks),]
ordered_peaks <- df$name[match(rownames(rna_input), df$external_gene_name)]
atac_input <- ave_access_sub[ordered_peaks,]

pdf("Plot_heatmap_ruben_gene_ATAC_region_acc_v2.pdf", height=10)
gplots::heatmap.2(atac_input, col = pink.col.heatmap,  trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 0.3, cexCol=1, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
dev.off()


# also use sidebar to show the gene category
ruben_genes <- list(
    "Cytokines"=c("CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL8", "CXCL16", "CCL20", "CCL28", "IL1A", "IL1B", "IL32"),
    "Cytokine_receptors"=c("CXCR4", "IL1RAP", "IL2RG", "IL4R", "IL6R", "IL6ST", "IL10RB", "IL17RA", "IL22RA1"),
    "NFKB"=c("REL", "RELB", "NFKB1", "NFKB2", "BIRC3", "LTB"),
    "Interferon_family"=c("IRF1", "IRF3", "IFNGR1", "IFNGR2"),
    "TNF"=c("TNF", "TNFAIP3", "TNIP1", "TNFSF15"),
    "Toll_like_receptor_pathway"=c("IRAK2", "TRAF6", "TLR2", "TLR3", "TLR4","TLR5"),
    "Antimicrobials_ROS"=c("LCN2", "SAA1", "SAA2", "DUOXA2", "DUOX2", "GPX2", "PIGR"),
    "Mucins_TFFs_glycosylation"=c("TFF1", "TFF2", "MUC1", "MUC13", "ST3GAL1", "ST3GAL4"),
    "Junctions_barrier_membrane"=c("ACTB", "ACTN4", "CLDN1", "CLDN23", "TJP1", "CEACAM1", "CEACAM5", "CEACAM7"),
    "Hypoxia_angiogenesis"=c("VEGFA"),
    "Others"=c("PLAUR", "PI3"),
    "Antigen_presentation"=c("CD74", "HLA-A", "B2M", "HLA-E", "TAPBP")
)

gene_anno <- setNames(rep(names(ruben_genes), sapply(ruben_genes, length)), unlist(ruben_genes))
row_sidebar_values <- gene_anno[rownames(rna_input)]
gene_cat <- c(
    "Cytokines",                  
    "Cytokine_receptors",        
    "NFKB" ,                      
    "Interferon_family" ,        
    "TNF" ,                       
    "Toll_like_receptor_pathway",
    "Antimicrobials_ROS" ,        
    "Mucins_TFFs_glycosylation" ,
    "Junctions_barrier_membrane",
    "Hypoxia_angiogenesis",                        
    "Antigen_presentation" , 
    "Others"              
)
gene_cat_colors <- c(
    "#A4D37190",
    "#A4D371",
    "#CA6778",
    "#7BCAA4",
    "#AB99CF",
    "#FDD884",
    "#A84798",
    "#5E4FA2",
    "#4DA7B0",
    "#197636",
    "#C12886",
    "#969696"

)
names(gene_cat_colors) <- gene_cat
row_sidebar_colors <- gene_cat_colors[row_sidebar_values]


pdf("Plot_heatmap_ruben_gene_expr_v2.pdf", height=10)
gplots::heatmap.2(rna_input, col = beach.col.heatmap,  RowSideColors=row_sidebar_colors, trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 0.3, cexCol=1, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
dev.off()


