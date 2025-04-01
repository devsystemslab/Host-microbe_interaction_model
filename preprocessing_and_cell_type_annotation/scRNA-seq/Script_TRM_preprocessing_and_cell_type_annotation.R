setwd("/home/yuq22/Bacteria_TRM/TRM/TRM_from_coculture")
library(Seurat)
library(ggplot2)
library(dplyr)
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")

# Load the data
folder <- "/projects/site/pred/ihb-g-deco/USERS/lopezsar/Bacteria_TRM_Manuscript/3rd Figure - Epithelium vs Epithelium+TRM"
x <- "TRM"
count <- Read10X(data.dir=paste0(folder, "/", x))
colnames(count) <- paste(x, colnames(count), sep="_")
combined_rna <- CreateSeuratObject(count)
combined_rna <- NormalizeData(object = combined_rna, normalization.method = "LogNormalize", scale.factor = 1e4)
# calculate Mt%
combined_rna[["percent.mt"]] <- PercentageFeatureSet(combined_rna, pattern = "^MT-")
combined_rna$log_nCount <- log(combined_rna$nCount_RNA)
saveRDS(combined_rna, file="Dat_merged_TRM_from_coculture.rds")
# show Mt%, nFeature, nCount distribution per sample
# Visualize QC metrics as a violin plot
p1 <- SCpubr::do_ViolinPlot(combined_rna, features = c("nFeature_RNA", "nCount_RNA", "log_nCount","percent.mt"), ncol = 4)
pdf("Plot_ViolinPlot_QC_metrics_per_sample.pdf" , width = 20, height = 5)
p1
dev.off()

x <- exp(10)
p1 <- SCpubr::do_ViolinPlot(combined_rna, features = "nCount_RNA")
p1+geom_hline(yintercept = x, linetype = "dashed", color = "red")
p1 <- SCpubr::do_ViolinPlot(combined_rna, features = "percent.mt")
p1+geom_hline(yintercept = 15, linetype = "dashed", color = "red")
p1 <- SCpubr::do_ViolinPlot(combined_rna, features = "nFeature_RNA")
p1+geom_hline(yintercept = 500, linetype = "dashed", color = "red")
# filter cells
cells <- colnames(combined_rna)[which(combined_rna$nFeature_RNA > 500 & combined_rna$nCount_RNA > 1000 & combined_rna$nCount_RNA<exp(10) & combined_rna$percent.mt < 15)]
seu_obj <- subset(combined_rna, cells=cells)
p2 <- SCpubr::do_ViolinPlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA", "log_nCount","percent.mt"), ncol = 4)
pdf("Plot_ViolinPlot_QC_metrics_per_sample_after_filtering.pdf" , width = 20, height = 5)
p2
dev.off()
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 3000)
seu_obj <- ScaleData(object = seu_obj, verbose = T)
seu_obj <- RunPCA(object = seu_obj, features = VariableFeatures(seu_obj), verbose = F, npcs = 50)
usefulPCs <- 1:20
seu_obj <- FindNeighbors(object = seu_obj, dims = usefulPCs)
seu_obj <- FindClusters(object = seu_obj, resolution = 1)
seu_obj <- RunUMAP(object = seu_obj, dims = usefulPCs)
saveRDS(seu_obj, file="Dat_merged_filtered_TRM_from_coculture.rds")

seu_obj <- readRDS("Dat_merged_filtered_TRM_from_coculture.rds")
# identify cell type markers
de_res <- presto::wilcoxauc(seu_obj, group_by="RNA_snn_res.1")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res %>% filter(padj<0.05 & pct_diff>10 & logFC>0.1)
top_res <- sig_res %>% group_by(group) %>% top_n(3, wt=logFC)
deg_res <- list("sig_res"=sig_res, "top_res"=top_res)
saveRDS(deg_res, file="Res_TRM_cell_type_markers.rds")

deg_res <- readRDS("Res_TRM_cell_type_markers.rds")
sig_res <- deg_res$sig_res
top_res <- sig_res[!grepl("ENSG", sig_res$feature),] %>% group_by(group) %>% top_n(3, wt=logFC)
genes <- unique(top_res$feature)
p2 <- SCpubr::do_DotPlot(sample = seu_obj, 
             features = genes,
             scale=T)
p2
ggsave(p2, filename="Plot_dotplot_TRM_cell_type_markers.pdf", width=9, height=6)

p2 <- SCpubr::do_DimPlot(seu_obj, group.by = "RNA_snn_res.1", pt.size=6, label=T, label.size=20, border.size=1.5)+NoLegend()
png("Plot_UMAP_colored_by_cluster.png", width = 2000, height = 2000)
p2
dev.off()

p1 <- SCpubr::do_ViolinPlot(seu_obj, features = "nFeature_RNA", group.by="RNA_snn_res.1")+NoLegend()
p2 <- SCpubr::do_ViolinPlot(seu_obj, features = "log_nCount", group.by="RNA_snn_res.1")+NoLegend()
p3 <- SCpubr::do_ViolinPlot(seu_obj, features = "percent.mt", group.by="RNA_snn_res.1")+NoLegend()
pdf("Plot_ViolinPlot_QC_metrics_per_cluster_after_filtering.pdf" , width = 5, height = 15)
p1/p2/p3
dev.off()

genes <- paste0("CD", c(3,4,8,69,"49a",103,"45RA",27,"45RO",39))
p2 <- SCpubr::do_DotPlot(sample = seu_obj, 
             features = genes,
             scale=T)
p2
ggsave(p2, filename="Plot_dotplot_TRM_CD_genes.pdf", width=3, height=5)

# exclude epi cluster
immune_obj <- subset(seu_obj, RNA_snn_res.1!=15)
SCpubr::do_DimPlot(immune_obj, group.by = "RNA_snn_res.1", pt.size=6, label=T, label.size=20, border.size=1.5)+NoLegend()
saveRDS(immune_obj, file="Dat_merged_filtered_TRM_from_coculture_no_epi.rds")
marker_list <- list(
    "NK" =c("XCL1","XCL2","GZMK","TYROBP","GNLY","KLRD1","GZMB"),
    "Stressed T"=c("DNAJB1","HSPA1A"),
    "Signalling T"=c("THEMIS","RASGRF2","AKT3","PAM"),
    "Naive T"=c("CCR7","SELL","LEF1","TSHZ2"),
    "Activated T"=c("NFKBIA", "CREM", "IL23A"),
    "CD4+ Treg"=c("TIGIT","IL2RA","FOXP3"),
    "Tis T"=c("ISG15","IFIT3","MX1"),
    "gamma-theta T"=c("TRDC","TRGC1"),
    "MAIT"=c("MAF","SLC4A10","KLRB1"),
    "Trm"=c("RORA","IL7R","CCR9","ITGA1","CD8A","PTPRK"),
    "General"=c("CD4", "CD27")
)
p2 <- SCpubr::do_DotPlot(sample = immune_obj, 
             features = marker_list,
             scale=T)
p2
ggsave(p2, filename="Plot_dotplot_TRM_natureTRM_paper_fig2_marker.pdf", width=10, height=7)


# get the cell type average expression
setwd("/home/yuq22/Bacteria_TRM/intestine_multi_region_scRNA-seq_atlas/immune_cells")
ct_expr <- getAveExpr(seu.obj=combined, feature.to.calc="Cell_type", colname.prefix=NULL)
saveRDS(ct_expr, file="Dat_gut_cell_atlas_adult_human_immune_cell_type_average_expression.rds")

# compare the transcriptome similarity between the immune cells from the TRM dataset and the immune cells from the tissue
setwd("/home/yuq22/Bacteria_TRM/TRM/TRM_from_coculture/")
dir.create("cell_type_annotation")
setwd("/home/yuq22/Bacteria_TRM/TRM/TRM_from_coculture/cell_type_annotation")
# load the TRM data
chip_trm <- readRDS("/home/yuq22/Bacteria_TRM/TRM/TRM_from_coculture/Dat_merged_filtered_TRM_from_coculture_no_epi.rds")
SCpubr::do_DimPlot(chip_trm, reduction = "umap", group.by = "RNA_snn_res.1", pt.size=4, border.size = 2, label=T, label.size=10, font.size=30, legend.icon.size=10)
# get cluster average expression
chip_trm_ct_expr <- getAveExpr(seu.obj=chip_trm, feature.to.calc="RNA_snn_res.1")
saveRDS(chip_trm_ct_expr, file="Dat_TRM_from_coculture_no_epi_cell_type_average_expression.rds")

expressed_genes <- intersect(rownames(ct_expr)[rowSums(ct_expr)>0], rownames(chip_trm_ct_expr)[rowSums(chip_trm_ct_expr)>0])

# read author provided cell type markers
library("readxl")
gut_cell_atlas_markers <- read_excel("/home/yuq22/Bacteria_TRM/intestine_multi_region_scRNA-seq_atlas/author_supplemental_tables/Table_S7_cell_type_markers.xls")
# subset to immune cell type markers
immune_markers <- gut_cell_atlas_markers[which(gut_cell_atlas_markers$lineage%in%c("myeloid","T_NK_cells","B_plasma_redbloodcells")),]
features <- intersect(immune_markers$genes, expressed_genes)
ref_expr <- ct_expr[features,]
query_expr <- chip_trm_ct_expr[features,]
scc_mat <- cor(ref_expr, query_expr, method="spearman")
pdf("Plot_heatmap_SCC_adult_tissue_immune_vs_TRM.pdf")
gplots::heatmap.2(scc_mat, scale="column", trace="none", dendrogram="both", col=darkBlue2Red.heatmap, density.info="none", keysize=0.8, mar=c(5,10))
dev.off()

max_ct <- setNames(rownames(scc_mat)[apply(scc_mat, 2, which.max)], colnames(scc_mat))
df <- data.frame(unique(immune_markers[,c("cluster", "lineage")]), stringsAsFactors = F)
max_lineage <- df[match(max_ct, df$cluster), "lineage"]


# only C9 belongs to B cells, others are T cells
# focus on comparison on T cells
# identify T cell subtype markers from the tissue data
combined <- readRDS("/home/yuq22/Bacteria_TRM/intestine_multi_region_scRNA-seq_atlas/immune_cells/Res_gut_cell_atlas_adult_human_immune_no_batch_correction.rds")
tissue_t_cells <- subset(combined, category=="T cells")
t_cell_expr <- getAveExpr(seu.obj=tissue_t_cells, feature.to.calc="Cell_type", colname.prefix=NULL, size.cutoff=1)
saveRDS(t_cell_expr, file="Dat_adult_tissue_T_cells_cell_type_average_expression_cell_type.rds")
# identify cell type markers
de_res <- presto::wilcoxauc(tissue_t_cells, group_by="Cell_type")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res %>% filter(padj<0.05 & pct_diff>10 & logFC>0.1)
top_res <- sig_res %>% group_by(group) %>% top_n(20, wt=logFC)
deg_res <- list("sig_res"=sig_res, "top_res"=top_res)
saveRDS(deg_res, file="Res_adult_tissue_T_cell_subtype_markers.rds")

chip_t_cells <- subset(chip_trm, RNA_snn_res.1!=9)
chip_t_cell_expr <- getAveExpr(seu.obj=chip_t_cells, feature.to.calc="RNA_snn_res.1")
saveRDS(chip_t_cell_expr, file="Dat_T_cells_from_coculture_no_epi_B_cell_type_average_expression.rds")

expressed_genes <- intersect(rownames(t_cell_expr)[rowSums(t_cell_expr)>0], rownames(chip_t_cell_expr)[rowSums(chip_t_cell_expr)>0])
features <- intersect(immune_markers$genes, expressed_genes)
ref_expr <- t_cell_expr[features,]
query_expr <- chip_t_cell_expr[features,]
scc_mat <- cor(ref_expr, query_expr, method="spearman")
pdf("Plot_heatmap_SCC_adult_tissue_vs_chip_T_cells.pdf")
gplots::heatmap.2(scc_mat, scale="column", trace="none", dendrogram="both", col=darkBlue2Red.heatmap, density.info="none", keysize=0.8, mar=c(5,10))
dev.off()

# subset to specific population and identify DEGs
seu_obj <- subset(chip_t_cells, cells=colnames(chip_t_cells)[which(chip_t_cells$RNA_snn_res.1%in%c(4,5,8))])
# identify cell type markers
de_res <- presto::wilcoxauc(seu_obj, group_by="RNA_snn_res.1")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res %>% filter(padj<0.05 & pct_diff>10 & logFC>0.1)
top_res <- sig_res %>% group_by(group) %>% top_n(10, wt=logFC)
genes <- unique(top_res$feature[which(top_res$group==8)])
p2 <- SCpubr::do_DotPlot(sample = chip_t_cells, 
             features = c("IFNGR1",""),
             scale=T, dot.scale=12)
p2


p1 <- SCpubr::do_DotPlot(sample = tissue_t_cells, 
             features = "SLC4A10",
             scale=T, dot.scale=12, group.by="Cell_type")
p1

ggsave(p2, filename="Plot_dotplot_TRM_cell_type_markers.pdf", width=20, height=20)

top_res <- sig_res %>% group_by(group) %>% top_n(100, wt=logFC)
features <- intersect(unique(top_res$feature), rownames(tissue_t_cells))
ref_expr <- t_cell_expr[features,]
query_expr <- chip_t_cell_expr[features,paste0("C", sort(unique(seu_obj$RNA_snn_res.1)))]
cor_mat <- cor(ref_expr, query_expr, method="spearman")
setNames(rownames(cor_mat)[apply(cor_mat, 2, which.max)], colnames(cor_mat))


gut_cell_atlas_genes <- c(
    "CD3D", "CD69", "FOS", "CD4","CCR7","SELL","CD8B","CD8A","CXCR3","TBX21","RORA","IFNG","TNF","AHNAK","CX3CR1","TRAV1-2","SLC4A10",
    "IFNGR1", "TRDC", "HOPX", paste0("TRGV", c(2,4,1,5,7,9),"TRDV2","PDCD1","FOXP3","IL2RA","TIGIT","IL7R","RORC","KIT","CXCR5","CCR6",
    "LST1","LINC00299","IL4I1","SCN1B","PYGL","ZBTB16","GSN","NCR2","CSF2","IL22","NCKAP1","IL17A","CA10","PTGDR2","HPGDS","IL9R","PRF1","NKG7",
    "GZMA","KIR2DL3")
)
chip_t_cell_expr <- readRDS("Dat_T_cells_from_coculture_no_epi_B_cell_type_average_expression.rds")
tissue_t_cell_expr <- readRDS("Dat_adult_tissue_T_cells_cell_type_average_expression_cell_type.rds")
expressed_genes <- intersect(rownames(tissue_t_cell_expr)[rowSums(tissue_t_cell_expr)>0], rownames(chip_t_cell_expr)[rowSums(chip_t_cell_expr)>0])

# identify cell type markers from the reference
seu_obj_list <- list("tissue"=tissue_t_cells, "chip"=chip_t_cells)
feature_vec <- setNames(c("Cell_type", "RNA_snn_res.1"),c("tissue", "chip"))
deg_res_list <- list()
for(x in names(seu_obj_list)){
    de_res <- presto::wilcoxauc(seu_obj_list[[x]], group_by=feature_vec[x])
    de_res$pct_diff <- de_res$pct_in - de_res$pct_out
    sig_res <- de_res %>% filter(padj<0.05 & pct_diff>10 & logFC>0.1)
    top_res <- sig_res %>% group_by(group) %>% top_n(200, wt=logFC)
    deg_res <- list("sig_res"=sig_res, "top_res"=top_res)
    deg_res_list[[x]] <- deg_res
}
tissue_markers <- unique(deg_res_list$tissue$top_res$feature)
chip_markers <- unique(deg_res_list$chip$top_res$feature)
denovo_markers <- intersect(tissue_markers, chip_markers)
features <- intersect(union(gut_cell_atlas_genes, denovo_markers), expressed_genes)
length(features)
ref_expr <- tissue_t_cell_expr[features,]
query_expr <- chip_t_cell_expr[features,]
scc_mat <- cor(ref_expr, query_expr, method="spearman")
max_ct <- setNames(rownames(scc_mat)[apply(scc_mat, 2, which.max)], colnames(scc_mat))

pdf("Plot_heatmap_SCC_adult_tissue_immune_vs_TRM_on_NatureGutCellAtlas_plus_denovo_markers_intersect.pdf", height=10, width=10)
gplots::heatmap.2(scc_mat, scale="column", cellnote=round(scc_mat,2), notecol="#202020", trace="none", dendrogram="both", col=darkBlue2Red.heatmap, density.info="none", keysize=0.8, mar=c(5,10))
dev.off()

# map the cell types between the two datasets based on 
# 1. transcriptome similarity using the de novo markers 1) from the reference; 2) from the query; 3) from the intersection of the two; 4) from the union of the two ; 5) paper provided markers; 6) union of the paper provided markers and the intersect of de novo markers 
# 2. significant overlapping between markers
# take the majority vote of all the comparisons

feature_list <- list(
    "tissue"=tissue_markers,
    "chip"=chip_markers,
    "intersect"=denovo_markers,
    "union"=union(tissue_markers, chip_markers),
    "paper"=gut_cell_atlas_genes,
    "paper_intersect"=intersect(gut_cell_atlas_genes, denovo_markers)
)
mapping_res <- list()
# map cell type based on transcriptome similarity
for(x in names(feature_list)){
    features <- intersect(feature_list[[x]], expressed_genes)
    ref_expr <- tissue_t_cell_expr[features,]
    query_expr <- chip_t_cell_expr[features,]
    scc_mat <- cor(ref_expr, query_expr, method="spearman")
    max_ct_list <- lapply(colnames(scc_mat), function(que_ct){
        idx <- order(scc_mat[,que_ct], decreasing=T)
        cor_diff_1_2 <- scc_mat[idx[1],que_ct] - scc_mat[idx[2],que_ct]
        cor_diff_1_3 <- scc_mat[idx[1],que_ct] - scc_mat[idx[3],que_ct]
        if(cor_diff_1_2>0.1){
            return(rownames(scc_mat)[idx[1]])
        }else if(cor_diff_1_3>0.1){
            max_ct <- c(rownames(scc_mat)[idx[1]], rownames(scc_mat)[idx[2]])
            return(max_ct)
        }else{
            max_ct <- c(rownames(scc_mat)[idx[1]], rownames(scc_mat)[idx[2]], rownames(scc_mat)[idx[3]])
            return(max_ct)
        }
    })
    names(max_ct_list) <- colnames(scc_mat)
    mapping_res[[x]] <- max_ct_list
}
# map cell type based on significant overlapping between markers


ref_expr <- tissue_t_cell_expr[features,]
query_expr <- chip_t_cell_expr[features,]
scc_mat <- cor(ref_expr, query_expr, method="spearman")
max_ct <- setNames(rownames(scc_mat)[apply(scc_mat, 2, which.max)], colnames(scc_mat))







p2 <- SCpubr::do_DotPlot(sample = chip_t_cells, 
             features = features,
             scale=T, dot.scale=12)

p1 <- SCpubr::do_DotPlot(sample = tissue_t_cells, 
             features = features,
             scale=T, dot.scale=12, group.by="Cell_type")
p_out <- p1+p2

ggsave(p_out, filename="Plot_dotplot_TRM_cell_type_markers.pdf", width=18, height=10)
setwd("/home/yuq22/group_share/Bacteria_TRM/TRM/TRM_from_coculture/cell_type_annotation/coarse_cell_type_annotation")
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(dplyr)
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")

# load TRM data
trm <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/TRM/TRM_from_coculture/Dat_merged_filtered_TRM_from_coculture_no_epi.rds")
# load Epi from coculture data
epi <- readRDS("/home/yuq22/Bacteria_TRM/TRM/Dat_merged_filtered_ileum_epithelium_TRM_coculture_noImmune_noLowQuality_clustered.rds")

# coarse grain cell type annotation on TRM data
SCpubr::do_DimPlot(trm, group.by="RNA_snn_res.1", label=T)
# based on the similarity to gut cell atlas immune cell types, merged C4/5/8, C1/2/3/6/7, repectively
cell_type_anno <- list(
    "B" = 9,
    "ILC3-like CD4 T"=0,
    "SELL+ CD4 T"=10,
    "CD4 type 1 T"=c(1,2,3,6,7),
    "CD8 type 1 T"=c(4,5,8),
    "CD8 Tmem/NK"=11,
    "TRGV3/4/5 CD8 T"=12,
    "ILC3"=13,
    "G2M CD4"=14
)
trm$Cell_type <- NA
for(x in names(cell_type_anno)){
    trm$Cell_type[which(trm$RNA_snn_res.1 %in% cell_type_anno[[x]])] <- x
}
p1 <- SCpubr::do_DimPlot(trm, group.by="Cell_type", label=T)+NoLegend()
ggsave(p1, filename="Plot_UMAP_TRM_with_coarse_grain_cell_type_anno.png", width=6, height=6)
saveRDS(trm, file="Dat_TRM_with_coarse_grain_cell_type_anno.rds")

trm <- readRDS("Dat_TRM_with_coarse_grain_cell_type_anno.rds")
p1 <- SCpubr::do_DimPlot(trm, group.by="Cell_type", label=T)+NoLegend()
ggsave(p1, filename="Plot_UMAP_TRM_with_coarse_grain_cell_type_anno.png", width=6, height=6)
# load tissue immune cell data
tissue <- readRDS("/home/yuq22/Bacteria_TRM/intestine_multi_region_scRNA-seq_atlas/immune_cells/Res_gut_cell_atlas_adult_human_immune_no_batch_correction.rds")
# subset to T cells
tissue_t_cells <- subset(tissue, category=="T cells")
tissue_t_cells$Cell_type <- as.character(tissue_t_cells$Integrated_05)
size <- sort(table(tissue_t_cells$Cell_type))
saveRDS(tissue_t_cells, file="Dat_gut_cell_atlas_adult_human_immune_T_cells.rds")
trm_t_cells <- subset(trm, Cell_type!="B")
saveRDS(trm_t_cells, file="Dat_TRM_T_cells.rds")

tissue_t_cells <- readRDS("Dat_gut_cell_atlas_adult_human_immune_T_cells.rds")
trm_t_cells <- readRDS("Dat_TRM_T_cells.rds")
# identify cell type markers in tissue and in vitro data separately
de_res_tissue <- presto::wilcoxauc(tissue_t_cells, group_by="Cell_type")
de_res_tissue$pct_diff <- de_res_tissue$pct_in - de_res_tissue$pct_out  
sig_res_tissue <- de_res_tissue %>% filter(padj<0.05 & pct_diff>10 & logFC>0.1)
sig_res_tissue$Cell_type <- paste("tissue", as.character(sig_res_tissue$group), sep=":")

de_res_trm <- presto::wilcoxauc(trm_t_cells, group_by="Cell_type")
de_res_trm$pct_diff <- de_res_trm$pct_in - de_res_trm$pct_out
sig_res_trm <- de_res_trm %>% filter(padj<0.05 & pct_diff>10 & logFC>0.1)
sig_res_trm$Cell_type <- paste("inVitro", as.character(sig_res_trm$group), sep=":")

sig_res <- rbind(sig_res_tissue, sig_res_trm)
saveRDS(sig_res, file="Res_tissue_and_inVitro_T_cell_subtype_markers.rds")

# get the overlap of cell type markers between tissue and in vitro data
tissue_ct <- sort(unique(grep("tissue:", sig_res$Cell_type, value=T)))
vitro_ct <- sort(unique(grep("inVitro:", sig_res$Cell_type, value=T)))
n1 <- sapply(tissue_ct, function(ct1){
    sapply(vitro_ct, function(ct2){
        tissue_markers <- sig_res$feature[which(sig_res$Cell_type==ct1)]
        vitro_markers <- sig_res$feature[which(sig_res$Cell_type==ct2)]
        length(intersect(tissue_markers, vitro_markers))
    })
})
# test the overlap significance
union_markers <- sort(unique(sig_res$feature))
pval <- sapply(tissue_ct, function(ct1){
    sapply(vitro_ct, function(ct2){
        tissue_ct1_markers <- sig_res$feature[which(sig_res$Cell_type==ct1)]
        vitro_ct2_markers <- sig_res$feature[which(sig_res$Cell_type==ct2)]
        overlap_markers <- intersect(tissue_ct1_markers, vitro_ct2_markers)
        ct1_non_ct2_markers <- setdiff(tissue_ct1_markers, vitro_ct2_markers)
        ct2_non_ct1_markers <- setdiff(vitro_ct2_markers, tissue_ct1_markers)
        non_ct1_ct2_markers <- setdiff(union_markers, c(tissue_ct1_markers, vitro_ct2_markers))
        a <- length(overlap_markers)
        b <- length(ct1_non_ct2_markers)
        c <- length(ct2_non_ct1_markers)
        d <- length(non_ct1_ct2_markers)
        fisher.test(matrix(c(a, b, c, d), nrow=2))$p.value
    })
})
saveRDS(pval, file="Res_cell_type_marker_overlap_significance_pval.rds")
pval <- readRDS("Res_cell_type_marker_overlap_significance_pval.rds")

res <- -log10(pval)-(-log10(0.05))
rownames(res) <- sub("inVitro:", "", rownames(res))
colnames(res) <- sub("tissue:", "", colnames(res))
pdf("Plot_heatmap_T_cell_type_marker_overlap_between_chip_and_tissue.pdf", height=10, width=10)
gplots::heatmap.2(res, scale="row", cellnote=round(res,1), notecol="#202020", trace="none", dendrogram="both", col=darkBlue2Red.heatmap, density.info="none", keysize=0.8, mar=c(10,10), cexRow=1, cexCol=1, xlab="Tissue", ylab="inVitro")
dev.off()

# for each invitro cell type, get the top 3 tissue cell type with the strongest overlap, unless the difference between the top 1 and top 2/3 is larger than 10
# do not include in vitro G2M CD4 for this analysis
vitro_ct <- sub("inVitro:", "", setdiff(all_vitro_ct, "inVitro:G2M CD4"))
matched_tissue_cts <- lapply(vitro_ct, function(ct){
    vec <- res[ct, ]
    top3 <- order(vec, decreasing=T)[1:3]
    if(vec[top3[1]]-vec[top3[2]]>10){
        idx <- top3[1]
    }else if(vec[top3[1]]-vec[top3[3]]>10){
        idx <- top3[1:2]
    }else{
        idx <- top3
    }
    names(vec)[idx]

})
names(matched_tissue_cts) <- vitro_ct
# CD4 cells won't be matched to any CD8 cells, so take top 1,2,4
matched_tissue_cts[["CD4 type 1 T"]] <- c(matched_tissue_cts[["CD4 type 1 T"]][1:2],"Treg")


# calculate the similarity of T cell subtypes between tissue and in vitro
tissue_expr <- getAveExpr(seu.obj=tissue_t_cells, feature.to.calc="Cell_type", colname.prefix="tissue:")
vitro_expr <- getAveExpr(seu.obj=trm_t_cells, feature.to.calc="Cell_type", colname.prefix="inVitro:")
saveRDS(tissue_expr, file="Res_tissue_T_cell_subtype_expression.rds")
saveRDS(vitro_expr, file="Res_inVitro_T_cell_subtype_expression.rds")

# get markers for each cell type
sig_res <- readRDS("Res_tissue_and_inVitro_T_cell_subtype_markers.rds")
vitro_marker_genes <- sort(unique(sig_res$feature[grep("inVitro", sig_res$Cell_type)]))
tissue_marker_genes <- sort(unique(sig_res$feature[grep("tissue", sig_res$Cell_type)]))
both_marker_genes <- intersect(vitro_marker_genes, tissue_marker_genes)

scc_mat <- cor(vitro_expr[both_marker_genes, ], tissue_expr[both_marker_genes, ], method="spearman")
res <- scc_mat
rownames(res) <- sub("inVitro:", "", rownames(res))
colnames(res) <- sub("tissue:", "", colnames(res))
pdf("Plot_heatmap_T_cell_type_SCC_between_chip_and_tissue.pdf", height=10, width=10)
gplots::heatmap.2(res, scale="row", cellnote=round(res,2), notecol="#202020", trace="none", dendrogram="both", col=darkBlue2Red.heatmap, density.info="none", keysize=0.8, mar=c(10,10), cexRow=1, cexCol=1, xlab="Tissue", ylab="inVitro")
dev.off()

# for each invitro cell type, get the top 3 most similar tissue cell type, unless the difference between the top 1 and top 2/3 is larger than 0.1
# do not include in vitro G2M CD4 for this analysis
vitro_ct <- sub("inVitro:", "", setdiff(all_vitro_ct, "inVitro:G2M CD4"))
scc_matched_tissue_cts <- lapply(vitro_ct, function(ct){
    vec <- res[ct, ]
    top3 <- order(vec, decreasing=T)[1:3]
    if(vec[top3[1]]-vec[top3[2]]>0.1){
        idx <- top3[1]
    }else if(vec[top3[1]]-vec[top3[3]]>0.1){
        idx <- top3[1:2]
    }else{
        idx <- top3
    }
    names(vec)[idx]

})
names(scc_matched_tissue_cts) <- vitro_ct
# CD4 cells won't be matched to any CD8 cells, so take top 1,2,4
scc_matched_tissue_cts[["CD4 type 1 T"]] <- c(scc_matched_tissue_cts[["CD4 type 1 T"]][1:2],"Activated CD4 T")

# get the overlap of overlap-based and similarity-based matched tissue cell types
overlap_matched_tissue_cts <- lapply(vitro_ct, function(ct){
    intersect(matched_tissue_cts[[ct]], scc_matched_tissue_cts[[ct]])
})
names(overlap_matched_tissue_cts) <- vitro_ct
# inVitro CD8 Tmem/NK also show high similarity to tissue NK cells (SCC=0.61, rank 4; rank 1 SCC=0.66)
# so manually match inVitro CD8 Tmem/NK to tissue NK
overlap_matched_tissue_cts[["CD8 Tmem/NK"]] <- "NK cell"

# update cell type annotation
ct1 <- c("CD4 type 1 T", "CD8 Tmem/NK", "CD8 type 1 T", "ILC3-like CD4 T", "SELL+ CD4 T", "TRGV3/4/5 CD8 T")
ct2 <- c("Activated CD4 T", "NK cell", "Activated CD8 T", "MAIT-like cell", "Treg", "TRGV3/4/5 gdT")
trm$Updated_cell_type <- trm$Cell_type
for(x in ct1){
    trm$Updated_cell_type[which(trm$Cell_type==x)] <- ct2[which(ct1==x)]
}
df <- unique(trm@meta.data[,c("Cell_type", "Updated_cell_type")])
trm$Cell_type <- NULL
p1 <- SCpubr::do_DimPlot(trm, group.by="Updated_cell_type", label=T)+NoLegend()
ggsave(p1, filename="Plot_UMAP_TRM_with_updated_coarse_grain_cell_type_anno.png", width=6, height=6)
saveRDS(trm, file="Dat_TRM_with_updated_cell_type_anno.rds")

# update the heatmap of cell type marker overlap
res <- -log10(pval)-(-log10(0.05))
rownames(res) <- sub("inVitro:", "", rownames(res))
colnames(res) <- sub("tissue:", "", colnames(res))
rownames(res) <- df$Updated_cell_type[match(rownames(res), df$Cell_type)]
pdf("Plot_heatmap_T_cell_type_marker_overlap_between_chip_and_tissue_updated_ct.pdf", height=10, width=10)
gplots::heatmap.2(res, scale="row", cellnote=round(res,1), notecol="#202020", trace="none", dendrogram="both", col=darkBlue2Red.heatmap, density.info="none", keysize=0.8, mar=c(10,10), cexRow=1, cexCol=1, xlab="Tissue", ylab="inVitro")
dev.off()

res <- scc_mat
rownames(res) <- sub("inVitro:", "", rownames(res))
colnames(res) <- sub("tissue:", "", colnames(res))
rownames(res) <- df$Updated_cell_type[match(rownames(res), df$Cell_type)]
pdf("Plot_heatmap_T_cell_type_SCC_between_chip_and_tissue_updated.pdf", height=10, width=10)
gplots::heatmap.2(res, scale="row", cellnote=round(res,2), notecol="#202020", trace="none", dendrogram="both", col=darkBlue2Red.heatmap, density.info="none", keysize=0.8, mar=c(10,10), cexRow=1, cexCol=1, xlab="Tissue", ylab="inVitro")
dev.off()

# for each invitro cell type, get the top 3 tissue cell type with the strongest overlap, unless the difference between the top 1 and top 2/3 is larger than 10
# do not include in vitro G2M CD4 for this analysis
vitro_ct <- sub("inVitro:", "", setdiff(all_vitro_ct, "inVitro:G2M CD4"))
matched_tissue_cts <- lapply(vitro_ct, function(ct){
    vec <- res[ct, ]
    top3 <- order(vec, decreasing=T)[1:3]
    if(vec[top3[1]]-vec[top3[2]]>10){
        idx <- top3[1]
    }else if(vec[top3[1]]-vec[top3[3]]>10){
        idx <- top3[1:2]
    }else{
        idx <- top3
    }
    names(vec)[idx]

})
names(matched_tissue_cts) <- vitro_ct
# CD4 cells won't be matched to any CD8 cells, so take top 1,2,4
matched_tissue_cts[["CD4 type 1 T"]] <- c(matched_tissue_cts[["CD4 type 1 T"]][1:2],"Treg")

# update the cell type annotation in the sig_res object
sig_res_2 <- sig_res
for(ct in df$Cell_type){
    sig_res_2$Cell_type[which(sig_res$Cell_type==paste0("inVitro:", ct))] <- paste("inVitro", df$Updated_cell_type[which(df$Cell_type==ct)], sep=":")
}
saveRDS(sig_res_2, file="Res_tissue_and_inVitro_T_cell_subtype_markers_ct_updated.rds")
sig_res <- sig_res_2
top_res <- sig_res %>% group_by(Cell_type) %>% top_n(20, logFC)
all_tissue_ct <- sort(unique(grep("tissue:", sig_res$Cell_type, value=T)))
all_vitro_ct <- sort(unique(grep("inVitro:", sig_res$Cell_type, value=T)))


tissue_ct <- "tissue:Activated CD4 T"
vitro_ct <- "inVitro:Activated CD4 T"
tissue_markers <- top_res$feature[which(top_res$Cell_type%in%tissue_ct)]
vitro_markers <- top_res$feature[which(top_res$Cell_type==vitro_ct)]
genes <- intersect(tissue_markers, vitro_markers)

tissue <- readRDS("/home/yuq22/Bacteria_TRM/intestine_multi_region_scRNA-seq_atlas/immune_cells/Res_gut_cell_atlas_adult_human_immune_no_batch_correction.rds")
# subset to non G2M B cells
cells <- colnames(tissue)[which(tissue$Cell_type%in%c("Memory B", "Naive B"))]
tissue_b_cells <- subset(tissue, cells=cells)
saveRDS(tissue_b_cells, file="Dat_gut_cell_atlas_adult_human_immune_B_cells.rds")
# identify cell type markers in tissue and in vitro data separately
de_res_tissue <- presto::wilcoxauc(tissue_b_cells, group_by="Cell_type")
de_res_tissue$pct_diff <- de_res_tissue$pct_in - de_res_tissue$pct_out  
sig_res_tissue <- de_res_tissue %>% filter(padj<0.05 & pct_diff>10 & logFC>0.1)
sig_res_tissue$Cell_type <- paste("tissue", as.character(sig_res_tissue$group), sep=":")
top_res <- sig_res_tissue %>% group_by(Cell_type) %>% top_n(10, logFC)
genes <- top_res$feature

# update cell type annotation in invivo data from "TRGV3/4/5 gdT" to "gdT"
trm$Updated_cell_type[which(trm$Updated_cell_type=="TRGV3/4/5 gdT")] <- "gdT"
trm$Updated_cell_type[which(trm$Updated_cell_type=="B")] <- "Memory B"
trm$Updated_cell_type[which(trm$Updated_cell_type=="MAIT-like cell")] <- "MAIT-like CD4"
saveRDS(trm, file="Dat_TRM_with_updated_cell_type_anno.rds")

# update the heatmap of cell type marker overlap
df$Updated_cell_type[which(df$Updated_cell_type=="TRGV3/4/5 gdT")] <- "gdT"
res <- -log10(pval)-(-log10(0.05))
rownames(res) <- sub("inVitro:", "", rownames(res))
colnames(res) <- sub("tissue:", "", colnames(res))
rownames(res) <- df$Updated_cell_type[match(rownames(res), df$Cell_type)]
pdf("Plot_heatmap_T_cell_type_marker_overlap_between_chip_and_tissue_updated_ct.pdf", height=10, width=10)
gplots::heatmap.2(res, scale="row", cellnote=round(res,1), notecol="#202020", trace="none", dendrogram="both", col=darkBlue2Red.heatmap, density.info="none", keysize=0.8, mar=c(10,10), cexRow=1, cexCol=1, xlab="Tissue", ylab="inVitro")
dev.off()

res <- scc_mat
rownames(res) <- sub("inVitro:", "", rownames(res))
colnames(res) <- sub("tissue:", "", colnames(res))
rownames(res) <- df$Updated_cell_type[match(rownames(res), df$Cell_type)]
pdf("Plot_heatmap_T_cell_type_SCC_between_chip_and_tissue_updated.pdf", height=10, width=10)
gplots::heatmap.2(res, scale="row", cellnote=round(res,2), notecol="#202020", trace="none", dendrogram="both", col=darkBlue2Red.heatmap, density.info="none", keysize=0.8, mar=c(10,10), cexRow=1, cexCol=1, xlab="Tissue", ylab="inVitro")
dev.off()

tissue_subset <- subset(tissue, cells=colnames(tissue)[which(tissue$category%in%c("B cells", "T cells"))])
marker_list <- list(
    "G2M CD4"=c("CD4", "MKI67"),
    "Memory B"=c("MZB1", "CD79B", "LINC01781", "CD27", "TNFRSF13B"),
    "ILC3"=c("KIT", "ZFP36L1"),
    "NK"=c("FCER1G", "GZMK",   "CTSW", "CCL4",   "TYROBP", "NKG7"),
    "MAIT"=c("SLC4A10", "KLRB1"),
    "Treg"=c("LTB", "ARID5B", "CCR7"),
    "gdT"=c("CD160", "GNLY", "KLRC2", "KLRD1", "TRG-AS1"),
    "Activated CD8 T"=c("XCL1", "CD8A", "CD8B", "IFNG"),
    "Activated CD4 T"=c("RGCC", "GPR183")
)

p1 <- SCpubr::do_DotPlot(trm, features=marker_list, group.by="Updated_cell_type")
p2 <- SCpubr::do_DotPlot(tissue_subset, features=marker_list, group.by="Cell_type")
p <- p1/p2
ggsave(p, filename="Plot_DotPlot_TRM_vs_tissue_immune_cell_type_markers.pdf", width=10, height=15)

# identify cell type markers
de_res <- presto::wilcoxauc(trm, group_by="Updated_cell_type")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res %>% filter(padj<0.05 & pct_diff>10 & logFC>0.1)
top_res <- sig_res %>% group_by(group) %>% top_n(20, wt=logFC)
deg_res <- list("sig_res"=sig_res, "top_res"=top_res)
saveRDS(deg_res, file="Res_duo_TRM_immune_cell_type_markers.rds")

# get cluster average expression
cluster_expr <- getAveExpr(seu.obj=trm, feature.to.calc="Updated_cell_type", colname.prefix=NULL)
saveRDS(cluster_expr, file="Dat_duo_TRM_immune_cell_type_expr.rds")

genes <- unique(top_res$feature)
expr_mat <- cluster_expr[genes,]
res <- prepareTreeAndHeatmapInput(expr.mat = expr_mat, hc.method = "ward.D2", norm.method = "quantile")
saveRDS(res, file="Res_duo_TRM_immune_cell_type_marker_expr_heatmap_input.rds")
pdf("Plot_heatmap_duo_TRM_immune_cell_type_marker_expr.pdf", height=10)
gplots::heatmap.2(res$heatmap_input, col = beach.col.heatmap, trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 1, cexCol=1, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
dev.off()

pdf("Plot_hc_duo_TRM_immune_cell_type_marker_expr.pdf")
plot(res$hc_col, main="Hierarchical clustering of adult tissue subcluster marker expression", hang=-1, cex=0.5)
dev.off()

# add cell type index to the trm data
idx <- setNames(paste(seq(ncol(res$heatmap_input)), colnames(res$heatmap_input), sep=":"), colnames(res$heatmap_input))
trm$Updated_cell_type_index <- idx[trm$Updated_cell_type]
trm$Cell_type <- trm$Updated_cell_type
saveRDS(trm, file="Dat_TRM_with_updated_cell_type_anno.rds")
saveRDS(trm, file="/home/yuq22/group_share/Bacteria_TRM/used_object/Dat_TRM_with_updated_cell_type_anno.rds")
p1 <- SCpubr::do_DimPlot(trm, group.by="Updated_cell_type_index")+NoLegend()
ggsave(p1, filename="Plot_UMAP_TRM_with_updated_coarse_grain_cell_type_anno.png", width=6, height=6)


# subset tissue immune data to those matched with in invitro data
trm$Source <- "inVitro"
selected_cell_type <- c("gdT", "ILC3", "Memory B", "Activated CD4 T", "Activated CD8 T", "Treg", "NK cell", "MAIT cell")
tissue_subset <- subset(tissue, cells=colnames(tissue)[which(tissue$Cell_type%in%selected_cell_type)])
tissue_subset$Source <- "tissue"

ct1 <- c("Activated CD4 T", "Activated CD8 T", "gdT", "ILC3", "MAIT cell", "Memory B", "NK cell", "Treg")  
ct2 <- setNames(c("8:Activated CD4 T", "3:Activated CD8 T", "4:gdT", "5:ILC3", "9:MAIT cell", "1:Memory B", "2:NK cell", "7:Treg"), ct1)
tissue_subset$Updated_cell_type_index <- ct2[as.character(tissue_subset$Cell_type)]
combined <- merge(x=trm, y=tissue_subset)
combined$Cell_type_per_source <- paste(combined$Updated_cell_type_index, combined$Source, sep=":")
combined$Cell_type_per_source <- factor(combined$Cell_type_per_source, levels=sort(unique(combined$Cell_type_per_source)))
Idents(combined) <- combined$Cell_type_per_source
saveRDS(combined, file="Dat_combined_immune_inVitro_and_tissue_data.rds")
saveRDS(combined, file="/home/yuq22/group_share/Bacteria_TRM/used_object/Dat_combined_immune_inVitro_and_tissue_data.rds")
marker_list <- list(
    "Memory B"=c("CD79B", "LINC01781", "TNFRSF13B"),
    "NK"=c("GZMB", "GZMK", "NKG7"),
    "Act. CD8"=c("XCL1", "CD8A", "CD8B"),
    "gdT"=c("CD160", "KLRC2", "KLRD1"),
    "ILC3"=c("KIT", "ZFP36L1"),
    "G2M"=c("CD4", "MKI67"),
    "Treg"=c("LTB", "ARID5B", "CCR7"),
    "Act. CD4"=c("RGCC", "GPR183"),
    "MAIT"=c("SLC4A10", "KLRB1")
)


p2 <- SCpubr::do_DotPlot(sample = combined, 
             features = marker_list,
             scale=T)
p2

ggsave(p2, filename="Plot_dotplot_adult_human_tissue_and_inVitro_immune_cell_type_markers.pdf", width=8, height=6)



# identify cell type markers
trm <- readRDS("/home/yuq22/group_share/Bacteria_TRM/used_object/Dat_TRM_with_updated_cell_type_anno.rds")
# focus on expressed genes not starting with ENSG and not from confounding genes
gene <- rownames(trm)[!grepl("^ENSG", rownames(trm))]
#remove ribosomal genes, mitochondrial genes, and genes located on sex chromosomes before normalization
confound_genes <- readRDS("/home/yuq22/ihb-intestine-evo/Annotation/confound_genes/Dat_confound_genes.rds")
genes_to_exclude <- unique(unlist(confound_genes[c(4,5,6)]))
genes <- setdiff(gene, genes_to_exclude)
X <- trm@assays$RNA@data[genes, ]
y <- trm$Updated_cell_type_index
de_res <- presto::wilcoxauc(X=X, y=y)
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res %>% filter(padj<0.05 & pct_diff>10 & logFC>0.1)
top_res <- sig_res %>% group_by(group) %>% top_n(20, wt=logFC)
deg_res <- list("sig_res"=sig_res, "top_res"=top_res)
saveRDS(deg_res, file="Res_duo_TRM_immune_cell_type_markers.rds")
write.table(top_res, file="Res_duo_TRM_immune_cell_type_markers.txt", sep="\t", quote=F, row.names=F)



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


