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
