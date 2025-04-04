setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/GRN_on_more_filtered_data/regulation_analysis_for_fig3")

library(Seurat)
library(Signac)
library(Pando)
library(ggplot2)
library(ggrepel)


genes <- c("CASP1", "CASP4", "CASP5", "CASP8", "IL18", #caspases, inflammosome
           "STAT3", "SOCS3", "SBNO2", "TNIP3", "DOC2B", #STAT3, and anti-inflammatory 
           "OSMR","ANXA1", "LAMP3", "WARS", "ABCA1", "IGFBP1","TRIM40", "MUC4","LURAP1L", "PELI2", "SPINK1", "CUX1", "SLC3A2","AHRR", #other functions
           "HSPA5","PDIA4", #ER stress
           "XDH", "NOS2","TXN", "TXNRD1", "SLC7A11", "PRKG1", "GCLM", "MOCOS", "AKR1C2", "AKR1C1", "NIBAN1","GPX2", #oxidative stress
           "CDC42BPA","FMNL2","RHOQ", "DSG3", "CLDN1", "CLIP4", "ARHGAP29", #cytoskeleton
           "PLA2G2A","LCN2", "DMBT1", "REG1A" #antimicrobial
           )


peaks <- "chr11-112156998-112157728"  # IL18
peaks <- "chr2-227817856-227819213" # CCL20
peaks <- c(
    "chr15-45112365-45114682", # DUOXA2
    "chr15-45091064-45091884", # DUOX2
    "chr15-45089457-45090469", # DUOX2
    "chr11-112156998-112157728"  # IL18
)

peaks <- "chr4-138259791-138260135" # SLC7A11
combined_rna <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/ATAC_on_more_filtered_data/motif_analysis/Res_salmonella_scMultiome_with_chromVAR_and_motif_res.rds")

motifs_in_regions <- names(which(combined_rna@assays$peaks_merged@motifs@data[peaks,]))
data(motif2tf)
tfs_in_regions <- sort(unique(motif2tf$tf[which(motif2tf$motif%in% motifs_in_regions)]))


# plot expression logFC of TFs in C16
all_fc_data <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/NFKB_barrier_innate_sensing/combined_with_Ruben_gene_set/Dat_expr_logFC_infection_control.rds")
c16_fc_tf <- all_fc_data[intersect(unique(motif2tf$tf), rownames(all_fc_data)),"Colonocyte_2_infection_specific_C16_C16_Salmonella"]
c17_fc_tf <- all_fc_data[intersect(unique(motif2tf$tf), rownames(all_fc_data)),"Colonocyte_2_infection_specific_C17_C17_Salmonella"]
c14_fc_tf <- all_fc_data[intersect(unique(motif2tf$tf), rownames(all_fc_data)),"Colonocyte_2_infection_specific_C14_C14_Salmonella"]

par(mfrow=c(3,2))
par(mar=c(15,5,5,5))
for(g in c("ZNF292", "BHLHE40", "STAT3", "BACH1", "RBPJ")){
    barplot(all_fc_data[g,], las=2, main=g)
}

tf_highlight <- names(c16_fc_tf)[order(c16_fc_tf, decreasing=T)[1:5]]
plot(seq(length(c16_fc_tf)), sort(c16_fc_tf, decreasing=T), xlab="Rank", ylab="expr. logFC (infection - control)", pch=16)
df <- data.frame(
  'tf_name'=names(c16_fc_tf),
  'logFC'=c16_fc_tf,
  'order'=rank(c16_fc_tf, decreasing=T),
  stringsAsFactors = F)
df$highlight <- df$tf_name %in% tf_highlight
p1 <- ggplot(df, aes(x=order, y=logFC)) +
  geom_point(shape=16, color="#303030") +
  ggrepel::geom_text_repel(data=df[df$highlight,], aes(label=tf_name), size=6, xlim  = c(100, NA), arrow = arrow(length = unit(0.01, "npc"), type = "closed", ends = "first")) +
  theme_minimal() +
  labs(y="expr. logFC(infection - control)", x="Rank")
ggsave(p1, filename="Plot_c16_tf_condition_expr_logFC.pdf")# the file is in /home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/GRN_on_more_filtered_data/GRN_on_Sal_colonocytes


# show selected genes in colonocytes between conditions
module_combined <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/GRN_on_more_filtered_data/incoorperate_multiple_GRN/Res_combined_GRN_modules_no_all_epi_cell_type_module.rds")
# get gene regualted by STAT3 or BACH1 and involved in the manually selected genes
gs1 <- module_combined$target[module_combined$tf%in%c("STAT3","BACH1") & module_combined$target %in% genes]

# load M3 genes
g <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/C16/gene_expr/Res_Salmonella_infection_DEG_cluster.rds")
# get gene regualted by STAT3 or BACH1 and defined as M3 DEGs
gs2 <- unique(module_combined$target[module_combined$tf%in%c("BACH1","STAT3") & module_combined$target %in% names(g)[which(g==3)]])


# highlighted genes in text with specific functions
gs3 <- c("CASP4", "CASP5", "IL18", "SOCS3", "SBNO2", "TNIP3")


# load KEGG pathway genes
path <- "/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/more_filterings/more_DEG_characterization/"
gene_list <- readRDS(paste(path, "Res_KEGG_and_HALLMARK_gene_list.rds", sep="/"))

selected_terms <- c("Amino sugar and nucleotide sugar metabolism", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "Fatty acid biosynthesis", "Protein processing in endoplasmic reticulum", "Ferroptosis", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE")
term_genes <- unique(unlist(gene_list[selected_terms]))
term_ds2 <- intersect(term_genes, gs2)


# load scMultiome data
combined_rna <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/ATAC_on_more_filtered_data/motif_analysis/Res_salmonella_scMultiome_with_chromVAR_and_motif_res.rds")
# subset to colonocytes
colonocytes <- subset(combined_rna, subset= Coarse_cell_type=="Colonocyte")
# exclude the per cell type per cluster group with less than 20 cells
size <- table(colonocytes$RNA_subtype_group_per_condition)
colonocytes <- subset(colonocytes, subset = RNA_subtype_group_per_condition %in% names(size)[which(size>20)])
saveRDS(colonocytes, file="~/Bacteria_TRM/used_object/Res_salmonella_scMultiome_with_chromVAR_and_motif_res_colonocytes.rds")
# generate dotplot of selected genes
DefaultAssay(colonocytes) <- "RNA"
gs <- c("LURAP1L",   "IL18",  "SBNO2","SLC7A11","CASP4", "CASP5","SOCS3","TNIP3", "TRIM40") # SEC31A is involved in unfolder protein response, and is a predicted target of BACH1, LURAP1L is the predicted target of STAT3 and BACH1, TRIM40 is the predicted target of STAT3
group_order <- c("Non_infection_specific_colonocyte:Control","Non_infection_specific_colonocyte:Salmonella", "Infection_specific_colonocyte_subtype_1:Salmonella", "Infection_specific_colonocyte_subtype_2:Salmonella")
colonocytes$RNA_subtype_group_per_condition <- factor(colonocytes$RNA_subtype_group_per_condition, levels=group_order)
p1 <- SCpubr::do_DotPlot(colonocytes,  features=gs, group.by="RNA_subtype_group_per_condition", scale=T, flip=T)
ggsave(p1, filename="Plot_dotplot_selected_C16_enriched_genes_v2.pdf", height=6, width=2)

# check whether hypoxia genes are regualted by BHLHE40
gene_list <- list(
    "Antimicrobial_activity" = c("LCN2", "DMBT1", "SAA1", "SAA2", "PLA2G2A", "PIGR"),
    "Oxidative_stress" = c("DUOX2", "DUOXA2", "XDH", "NOS2", "GPX2"),
    "Cytokine" = c("CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL8", "CCL20", "CCL28", "IL32", "IL18", "IL1A", "IL1B"),
    "Cytokine_receptors" = c("CXCR4", "IL4R", "IL22RA1", "IL17RA", "IL6R"),
    "Cytoskeleton_Junctions" = c("ACTB", "RAC1", "CDC42", "CDC42BPA", "DST",  "CEACAM1", "CEACAM5", "CEACAM7", "KRT18"),
    "Hypoxia_Angiogenesis" = c("VEGFA", "HIF1A", "ADM", "NDRG1", "HILPDA", "ANKRD37"),
    "Interferon family" = c("IFNGR1", "IFNGR2", "IRF1", "IRF3"),
    "NFKB pathway" = c("RELB", "NFKB2", "BIRC3", "REL", "NFKB1", "TLR4", "TLR3", "MYD88", "TRAF6", "IRAK2"),
    "Mucus_TFFs" = c("MUC1", "MUC13", "MUC4", "MUC2", "TFF1", "TFF2"),
    "TNF" = c("TNF", "LTB", "TNFAIP3", "TNIP1"),
    "Others" = c("PI3", "PLAUR", "PLCG2", "BMP2", "ERBIN", "ANXA1", "CD74", "B2M", "HLA-A", "HLA-E", "EGFR", "SLC2A1")
)

module_combined[which(module_combined$tf=="BHLHE40" & module_combined$target %in% unlist(gene_list)),]
hypoxia_regulating_tf <- sort(unique(module_combined$tf[which(module_combined$target %in% gene_list[["Hypoxia_Angiogenesis"]])]))
intersect(hypoxia_regulating_tf, names(g)[which(g%in%c(1,6))])

# load scMultiome data
colonocytes <- readRDS("~/Bacteria_TRM/used_object/Res_salmonella_scMultiome_with_chromVAR_and_motif_res_colonocytes.rds")
combined_rna <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/ATAC_on_more_filtered_data/motif_analysis/Res_salmonella_scMultiome_with_chromVAR_and_motif_res.rds")
# exclude the per cell type per cluster group with less than 20 cells
size <- table(combined_rna$RNA_subtype_group_per_condition)
combined_rna_sub <- subset(combined_rna, subset = RNA_subtype_group_per_condition %in% names(size)[which(size>20)])
# exclude TA cells, BEST4+ cells and EEC cells
cells <- colnames(combined_rna_sub)[which(!combined_rna_sub$Coarse_cell_type %in% c("TA", "BEST4+_cell", "EEC"))]
combined_rna_sub <- subset(combined_rna_sub, cells=cells)
table(combined_rna_sub$RNA_subtype_group_per_condition)
saveRDS(combined_rna_sub, file="~/Bacteria_TRM/used_object/Res_salmonella_scMultiome_with_chromVAR_and_motif_res_excluding_small_size_groups_and_cell_type.rds")
seu_obj <- combined_rna_sub


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
group <- sort(unique(as.character(seu_obj$RNA_subtype_group_per_condition)))
mat <- do.call('rbind', strsplit(group, split=":"))
ct <- mat[,1]
names(ct) <- group
group_cols <- setNames(colors$cell_type[ct], group)
source("/home/yuq22/ihb-intestine-evo/common_script/plotting/Script_plotting_functions.R")
#peaks <- c(
#    "chr15-45112365-45114682", # DUOXA2
#    "chr15-45091064-45091884", # DUOX2
#    "chr15-45089457-45090469", # DUOX2
#    "chr11-112156998-112157728"  # IL18
#)
#selected_genes <- c("NQO1", "GCLC", "GCLM", "SLC7A11", "HO1")
selected_genes <- "VEGFA"
setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/GRN_on_more_filtered_data/regulation_analysis_for_fig3/")
peak_anno <- readRDS('/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/peak_annotation/Data_frame_peak_gene_pairs.rds')
idx <- which(peak_anno$external_gene_name %in% selected_genes)
peaks <- peak_anno$name[idx]
length(peaks)


p_list <- coverage_plot_with_evo_sig(seu_obj=seu_obj, 
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
                                     plot_suffix="-1")

# get the position of the GWAS regions
extend_up <- 0
extend_down <- 0
p_human <- "chr15-45112365-45114682"
p_human_vec <- strsplit(p_human, split="-")[[1]]
region <- GRanges(
  seqnames=p_human_vec[1],
  ranges=IRanges(
    start=as.numeric(p_human_vec[2])-extend_up,
    end=as.numeric(p_human_vec[3])+extend_down
  )
)
## plot the overlapped sites
pos_list <- lapply(names(evo_list), function(x){
  # subset to covered range
  peak.intersect <- subsetByOverlaps(x = evo_list[[x]], ranges = region)
  peak.df <- as.data.frame(x = peak.intersect)
})
names(pos_list) <- names(evo_list)


gwas_site <- list(
    "IL18"=c("11:112154583", "11:112153104", "11:112154583")
)

peaks <- "chr4-138259791-138260135" # SLC7A11
motifs_in_regions <- names(which(colonocytes@assays$peaks_merged@motifs@data[peaks,]))
data(motif2tf)
tfs_in_regions <- sort(unique(motif2tf$tf[which(motif2tf$motif%in% motifs_in_regions)]))

# update heatmap of general response genes
setwd("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/salmonella/multiome_analysis/GRN_on_more_filtered_data/regulation_analysis_for_fig3/")
ruben_genes <- list(
    "Antimicrobial_activity" = c("LCN2", "DMBT1", "SAA1", "SAA2", "PLA2G2A", "PIGR"),
    "Oxidative_stress" = c("DUOX2", "DUOXA2", "XDH", "NOS2", "GPX2"),
    "Cytokine" = c("CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL8", "CCL20", "CCL28", "IL32", "IL18", "IL1A", "IL1B"),
    "Cytokine_receptors" = c("CXCR4", "IL4R", "IL22RA1", "IL17RA", "IL6R"),
    "Cytoskeleton_Junctions" = c("ACTB", "RAC1", "CDC42", "CDC42BPA", "DST",  "CEACAM1", "CEACAM5", "CEACAM7", "KRT18"),
    "Hypoxia_Angiogenesis" = c("VEGFA", "HIF1A", "ADM", "NDRG1", "HILPDA", "ANKRD37"),
    "Interferon_family" = c("IFNGR1", "IFNGR2", "IRF1", "IRF3"),
    "NFKB_pathway" = c("RELB", "NFKB2", "BIRC3", "REL", "NFKB1", "TLR4", "TLR3", "MYD88", "TRAF6", "IRAK2"),
    "Mucus_TFFs" = c("MUC1", "MUC13", "MUC4", "MUC2", "TFF1", "TFF2"),
    "TNF" = c("TNF", "LTB", "TNFAIP3", "TNIP1"),
    "Others" = c("PI3", "PLAUR", "PLCG2", "BMP2", "ERBIN", "ANXA1", "CD74", "B2M", "HLA-A", "HLA-E", "EGFR", "SLC2A1")
)

# also use sidebar to show the gene category
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")
gene_cat <- c(
    "Cytokine",                  
    "Cytokine_receptors",        
    "NFKB_pathway" ,                      
    "Interferon_family" ,        
    "TNF" ,                       
    "Antimicrobial_activity" ,        
    "Mucus_TFFs" ,
    "Hypoxia_Angiogenesis",                        
    "Cytoskeleton_Junctions",
    "Oxidative_stress",   
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
    #"#C12886",
    #"#F07C28",
    "#969696"

)
names(gene_cat_colors) <- gene_cat

col_order <- c(
"Colon_SC_Salmonella", 
"Colon_SC_Control",                                  
"Colonocyte_2_infection_specific_C14_Salmonella",
"Colonocyte_2_infection_specific_C17_Salmonella",
"Colonocyte_2_infection_specific_C16_Salmonella",
"Non-infection-specific Colonocyte_Salmonella", 
"Non-infection-specific Colonocyte_Control",                        
"Goblet_Salmonella", 
"Goblet_Control"                                    
)  
ave_expr <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/NFKB_barrier_innate_sensing/combined_with_Ruben_gene_set/Dat_SC_GC_Colonocyte_cell_type_group_per_condition_ave_expr.rds")
genes <- unlist(ruben_genes)
expr_mat <- ave_expr[genes, col_order]
res <- prepareTreeAndHeatmapInput(expr.mat = expr_mat, hc.method = "ward.D2",  norm.method = "max", column.reorder = FALSE)
saveRDS(res, file="Res_Ruben_genes_heatmap_input.rds")
# also use sidebar to show the gene category
gene_anno <- setNames(rep(names(ruben_genes), sapply(ruben_genes, length)), unlist(ruben_genes))
row_sidebar_values <- gene_anno[rownames(res$heatmap_input)]
row_sidebar_colors <- gene_cat_colors[row_sidebar_values]
pdf("Plot_heatmap_ruben_gene_expr_v4_gene_group_mixed_subtypes_max_normalized.pdf", height=10)
gplots::heatmap.2(res$heatmap_input, col = beach.col.heatmap,  RowSideColors=row_sidebar_colors, trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 0.3, cexCol=1, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
dev.off()

res <- prepareTreeAndHeatmapInput(expr.mat = expr_mat, hc.method = "ward.D2",  norm.method = "quantile", column.reorder = FALSE)
saveRDS(res, file="Res_Ruben_genes_heatmap_input_min-max_normalized.rds")
# also use sidebar to show the gene category
pdf("Plot_heatmap_ruben_gene_expr_v4_gene_group_mixed_subtypes_min-max_normalized.pdf", height=10)
gplots::heatmap.2(res$heatmap_input, col = beach.col.heatmap,  RowSideColors=row_sidebar_colors, trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 0.3, cexCol=1, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
dev.off()

# generate dot plot of selected genes
combined_rna_sub <- readRDS("~/Bacteria_TRM/used_object/Res_salmonella_scMultiome_with_chromVAR_and_motif_res_excluding_small_size_groups_and_cell_type.rds")
genes <- c("CASP1", "CASP8", #caspases
           "STAT3", "BACH1", "CUX1", #TFs
           "ANXA1", "WARS", "ABCA1", "IGFBP1", #other functions
           "XDH", "TXNRD1", "GCLM", "MOCOS", "NIBAN1", #oxidative stress
           "CDC42BPA","RHOQ", "ARHGAP29", #cytoskeleton
           "PLA2G2A", "DMBT1", #antimicrobial
           "NQO1", "GCLC",# reported BACH1 targets
           "SEC31A", "RPS6KA5", "CTNNA1", "NT5C2" # "CTNNA1" "NT5C2"  "SEC31A" are predicted BACH1 targets, "RPS6KA5" is the predicted STAT3 target
           )

# get BACH1 and STAT3 predicted targets that also belong to M3
module_combined <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/GRN_on_more_filtered_data/incoorperate_multiple_GRN/Res_combined_GRN_modules_no_all_epi_cell_type_module.rds")
g <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/C16/gene_expr/Res_Salmonella_infection_DEG_cluster.rds")
m3_targets <- unique(module_combined$target[module_combined$tf%in%c("BACH1","STAT3") & module_combined$target %in% names(g)[which(g==3)]])
#selected_genes <- union(genes, m3_targets)
bach1_targets <- unique(module_combined$target[module_combined$tf%in%c("BACH1")])
stat3_targets <- unique(module_combined$target[module_combined$tf%in%c("STAT3")])
intersect(bach1_targets, genes) #  
intersect(stat3_targets, genes) # "RPS6KA5"

selected_genes <- genes
length(selected_genes)
# generate dotplot of selected genes
DefaultAssay(combined_rna_sub) <- "RNA"
p1 <- SCpubr::do_DotPlot(combined_rna_sub,  features=selected_genes, group.by="RNA_subtype_group_per_condition", scale=T, flip=T)
ggsave(p1, filename="Plot_dotplot_selected_C16_enriched_genes_v3.pdf", height=9, width=3)
