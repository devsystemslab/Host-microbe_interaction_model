setwd("/home/yuq22/Bacteria_TRM/salmonella/RNA_analysis/")
library(Seurat)
library(Signac)
library(dplyr)
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")

#data_folder <- "/home/yuq22/Bacteria_TRM/salmonella/data/count_res/filtered_feature_bc_matrix"
#atac_anno <- readRDS("~/ihb-intestine-evo/Annotation/EnsDB_for_ATAC/data.ens98_annot_for_atac.rds")
#
## create a list to store the Seurat object
## first do analysis on RNA component
#sample_id <- list.files(data_folder)
#seu_obj_list <- list()
#for(x in sample_id){
#  path <- file.path(data_folder, x, "matrix.h5")
#  count_list <- Read10X_h5(path)
#  rna_count <- count_list$`Gene Expression`
#  id <- strsplit(x, split="_")[[1]][3]
#  colnames(rna_count) <- paste(id, colnames(rna_count), sep="_")
#  # Create Seurat object
#  seu_obj <- CreateSeuratObject(counts = rna_count)
#  group <- substr(id,1,nchar(id)-1)
#  seu_obj$condition <- group
#  seu_obj_list[[id]] <- seu_obj
#  #atac_count <- count_list$Peaks
#} 

data_folder <- "/home/yuq22/Bacteria_TRM/salmonella/data_2"
sample_id <- list.files(data_folder)
seu_obj_list <- list()
for(x in sample_id){
  path <- file.path(data_folder, x)
  rna_count <- Read10X_h5(path)
  id <- strsplit(x, split="_")[[1]][1]
  colnames(rna_count) <- paste(id, colnames(rna_count), sep="_")
  # Create Seurat object
  seu_obj <- CreateSeuratObject(counts = rna_count)
  group <- substr(id,1,nchar(id)-1)
  seu_obj$condition <- group
  seu_obj_list[[id]] <- seu_obj
} 


# RNA analysis
combined_rna <- merge(x=seu_obj_list[[1]], y=seu_obj_list[-1])
combined_rna <- NormalizeData(object = combined_rna, normalization.method = "LogNormalize", scale.factor = 1e4)
# calculate Mt%
combined_rna[["percent.mt"]] <- PercentageFeatureSet(combined_rna, pattern = "^MT-")
saveRDS(combined_rna, file="Dat_merged_ileum_epithelium_infection_ds2.rds")

#combined_rna$log_nCount <- log(combined_rna$nCount_RNA)
# show Mt%, nFeature, nCount distribution per sample
# Visualize QC metrics as a violin plot
p1 <- SCpubr::do_ViolinPlot(combined_rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pdf("Plot_ViolinPlot_QC_metrics_per_sample_ds2.pdf" , width = 20, height = 5)
p1
dev.off()
# filter cells
cells <- colnames(combined_rna)[which(combined_rna$nFeature_RNA > 200 & combined_rna$nCount_RNA<30000)]
seu_obj <- subset(combined_rna, cells=cells)
#p2 <- SCpubr::do_ViolinPlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA", "log_nCount","percent.mt"), ncol = 4)
#pdf("Plot_ViolinPlot_QC_metrics_per_sample_after_filtering.pdf" , width = 20, height = 5)
#p2
#dev.off()
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 3000)
seu_obj <- ScaleData(object = seu_obj, verbose = T)
seu_obj <- RunPCA(object = seu_obj, features = VariableFeatures(seu_obj), verbose = F, npcs = 50)
usefulPCs <- 1:20
seu_obj <- FindNeighbors(object = seu_obj, dims = usefulPCs)
seu_obj <- FindClusters(object = seu_obj, resolution = 1)
seu_obj <- RunUMAP(object = seu_obj, dims = usefulPCs)
saveRDS(seu_obj, file="Dat_merged_filtered_ileum_epithelium_infection_ds2.rds")

# examine the regional identity the epithelial cells
p1 <- SCpubr::do_DimPlot(seu_obj, group.by = "condition", pt.size=6, font.size=100, border.size=1.5, legend.icon.size=40)
p2 <- SCpubr::do_DimPlot(seu_obj, group.by = "RNA_snn_res.1", pt.size=6, label=T, label.size=20, border.size=1.5)+NoLegend()
png("Plot_UMAP_colored_by_condition_and_cluster_ds2.png", width = 2000*2, height = 2000)
p1+p2
dev.off()

plotFeature(seu.obj=seu_obj, dr="umap", col.num=5, genes.to.plot=c("EPCAM", "PTPRC", "PDX1", "ONECUT2", "CDX2", "OSR2", "HOXB7", "HOXA10","SATB2","FABP6", "APOA4","BEST4","RBPJ","LGR5","","MUC2","CHGA","MKI67","OLFM4","APOC3","AQP10"))
plotFeature(seu.obj=seu_obj, dr="umap", col.num=5, genes.to.plot=c("EPCAM", "PTPRC", "PDX1", "ONECUT2", "CDX2", "OSR2", "HOXB7", "HOXA10","SATB2","FABP6", "APOA4","BEST4","RBPJ","LGR5","","MUC2","CHGA","MKI67","OLFM4","APOC3","AQP10"))
SCpubr::do_FeaturePlot(seu_obj, split.by="group", feature="CCL20", order=T)


