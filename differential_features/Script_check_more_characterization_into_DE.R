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
library(ggpubr)

setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/more_filterings/more_DEG_characterization")
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")

rna <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/more_filterings/Dat_salmonella_colon_epithelium_scRNA-seq_data_with_subtypes.rds")
SCpubr::do_DimPlot(rna, group.by="RNA_snn_res.1", label=T)
rna$cluster_group <- rna$Subtype
for(i in c(14, 16, 17)){
    rna$cluster_group[which(rna$RNA_snn_res.1==i)] <- paste0("Colonocyte_2_infection_specific_C",i)
}
rna$cluster_group_per_condition <- paste0(rna$cluster_group, "_", rna$condition)
rna$Subtype_v2 <- rna$Subtype
rna$Subtype_v2[grep("Colonocyte_2_infection_specific_C", rna$cluster_group)] <- "Colonocyte_2_infection_specific"
saveRDS(rna, file="/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/more_filterings/Dat_salmonella_colon_epithelium_scRNA-seq_data_with_subtypes.rds")

# generate stacked bar plot of cell type composition per condition
n1 <- table(rna$condition, rna$cluster_group)
prop1 <- prop.table(n1, margin=1)
df <- data.frame('cell_type'=rep(colnames(n1), each=nrow(n1)),
                 'condition'=rep(rownames(n1), ncol(n1)),
                 'prop'=as.vector(prop1),
                 stringsAsFactors = F)
df$cell_type <- factor(df$cell_type, levels=colnames(n1)[c(1:3,5:7,4,8:ncol(n1))])
cols <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/used_object/cols.rds")
#ct_cols <- cols$cell_type
##ct_cols <- c(ct_cols, setNames("#a93226", "Colonocyte_2_infection_specific"))
#ct_cols <- c(ct_cols, setNames(c("#ec7063", "#cb4335", "#a93226"), c("Colonocyte_2_infection_specific_C14", "Colonocyte_2_infection_specific_C16", "Colonocyte_2_infection_specific_C17")))
#cols$cell_type <- ct_cols
#saveRDS(cols, file="/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/used_object/cols.rds")

p1 <- ggplot(df, aes(x=condition, y=prop, fill=cell_type))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=ct_cols[unique(rna$cluster_group)])+
    theme_light()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    labs(x="Condition", y="Cell number", fill="Cell type")
p1
ggsave(p1, filename="Plot_cell_type_composition_per_condition.pdf", width=5, height=5)

# generate UMAP on updated cell type annotation
p1 <- SCpubr::do_DimPlot(rna, group.by="cluster_group", colors.use=ct_cols, border.size=1.5)+NoLegend()
p1
ggsave(p1, filename="Plot_UMAP_cell_type_annotation.png", width=5, height=5)

p2 <- SCpubr::do_DimPlot(rna, group.by="condition", colors.use=setNames(c("#969696", "#202020"), c("Control", "Salmonella")), border.size=1.5)+NoLegend()
p2
ggsave(p2, filename="Plot_UMAP_RNA_condition.png", width=5, height=5)



# get per cluster group per condition average gene expression
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")
ave_expr <- getAveExpr(seu.obj=rna, feature.to.calc="cluster_group_per_condition", size.cutoff=20, colname.prefix=NULL)
saveRDS(ave_expr, file="Dat_ave_expr_per_cluster_group_per_condition.rds")

data('motifs')
data('motif2tf')
tfs <- sort(unique(motif2tf$tf))

# load non-TA s2E DEGs induced upon infection
combined_sig_res <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/more_filterings/Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_matched_control_cells_and_cocluster_control_cells_DEG_sig_res.rds")
do.call('rbind', strsplit(combined_sig_res$group, split="_")) %>% as.data.frame() %>% setNames(c("cluster_group", "condition")) -> df
combined_sig_res <- cbind(combined_sig_res, df)
# group the DEGs per cluster into cell groups
df <- unique(rna@meta.data[,c("cluster_group", "RNA_snn_res.1")])
df[,2] <- paste0("C", df[,2])
vec <- setNames(df[,1], df[,2])
combined_sig_res$cell_group <- vec[combined_sig_res$cluster_group]
combined_sig_res$cell_group_per_condition <- paste(combined_sig_res$cell_group, combined_sig_res$condition, sep="_")
saveRDS(combined_sig_res, file="Res_Salmonella_infection_on_colon_epi_per_salmonella_cluster_vs_matched_control_cells_and_cocluster_control_cells_DEG_sig_res.rds")
de_tf_res <- combined_sig_res[which(combined_sig_res$feature%in%tfs),]

# load goblet cell DEGs
goblet_deg_res <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/more_filterings/goblet/Res_Salmonella_infection_on_colon_goblet_cluster_DEG.rds")
goblet_sig_res <- goblet_deg_res %>% filter(padj<0.05 & pct_in>10 & logFC>0.1)
goblet_sig_res$cell_group <- "Goblet"
goblet_sig_res$condition <- sapply(goblet_sig_res$group, function(x){strsplit(x, split="_")[[1]][1]})
goblet_sig_res$cell_group_per_condition <- paste(goblet_sig_res$cell_group, goblet_sig_res$condition, sep="_")
colnames(combined_sig_res)[15] <- "Cluster"
epi_sig_res <- data.frame(rbind(combined_sig_res[,colnames(goblet_sig_res)], goblet_sig_res))
saveRDS(epi_sig_res, file="Res_Salmonella_infection_on_colon_epi_goblet_and_s2E_DEG_sig_res.rds")
saveRDS(epi_sig_res, file="/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/used_object/Res_Salmonella_infection_on_colon_epi_goblet_and_s2E_DEG_sig_res.rds")

# plot Salmonella upregulated genes per cell group
epi_sig_res <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/used_object/Res_Salmonella_infection_on_colon_epi_goblet_and_s2E_DEG_sig_res.rds")
epi_sal_up_res <- epi_sig_res[epi_sig_res$condition=="Salmonella",]
pdf("Plot_barplot_salmonella_upregulated_genes_per_cell_group.pdf", height=5)
barplot(sqrt(table(epi_sal_up_res$cell_group_per_condition)), las=2)
dev.off()

# generat barplot to show DEG numbers, separate by cell type and direction
epi_sig_res <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/used_object/Res_Salmonella_infection_on_colon_epi_goblet_and_s2E_DEG_sig_res.rds")
num1 <- sapply(sort(unique(epi_sig_res$cell_group)), function(i){
    sapply(c("Salmonella", "Control"), function(j){
        sqrt(sum(epi_sig_res$condition==j & epi_sig_res$cell_group==i))
    })
})

pdf("Plot_barplot_salmonella_DEG_number_per_cell_group.pdf", height=5)
barplot(num1["Salmonella", ], las=2, ylim=c(-max(num1["Control", ]), max(num1["Salmonella", ])), col="#202020")
barplot(-num1["Control", ], las=2, names="", add=T, col="#979697")
legend("topright", legend=c("Salmonella", "Control"), fill=c("#202020", "#979697"))
dev.off()


SCpubr::do_ViolinPlot(rna, group.by="cluster_group", feature="nCount_RNA")

input <- ave_expr[sort(unique(de_tf_res$feature)),]
res <- prepareTreeAndHeatmapInput(expr.mat = input, hc.method = "complete", norm.method = "quantile")
saveRDS(res, file="Res_salmonella_DE_TF_heatmap_input.rds")
pdf("Plot_heatmap_salmonella_DE_TF_expression.pdf", height=10)
gplots::heatmap.2(res$heatmap_input, col = beach.col.heatmap, trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 1, cexCol=0.5, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
dev.off()

pdf("Plot_hc_salmonella_DE_TF_expression.pdf")
plot(res$hc_col, main="Hierarchical clustering of adult tissue subcluster marker expression", hang=-1, cex=0.5)
dev.off()

# calculate gene set score on the epithelium before and after infection
library(msigdbr)
library(ggpubr)

all_gene_sets <- msigdbr(species = "Homo sapiens")
h_genes <- msigdbr(species = "Homo sapiens", category = "H")
h_gene_list <- tapply(h_genes$gene_symbol, h_genes$gs_name, list)
hsp_genes <- read.table("/home/yuq22/group_share/Annotation/HGNC_gene_family/Table_heat_shock_protein_genes.txt", sep="\t", stringsAsFactors=F, head=T, quote="")
hsp_genes <- sort(unique(hsp_genes$Approved.symbol))
h_gene_list[["HGNC_hsp_genes"]] <- hsp_genes
kegg_gene_list <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Annotation/KEGG/Res_KEGG_gene_list.rds")
gene_list <- c(h_gene_list, kegg_gene_list)
saveRDS(gene_list, file="Res_KEGG_and_HALLMARK_gene_list.rds")
selected_terms <- names(gene_list)
rna <- AddModuleScore(rna, features=gene_list, name="gs")
colnames(rna@meta.data)[grep("gs", colnames(rna@meta.data))] <- selected_terms
saveRDS(rna, file="/home/yuq22/Bacteria_TRM/used_object/Dat_salmonella_colon_epithelium_scRNA-seq_data_with_subtypes_and_feature_scores.rds")

# identify marker genes of cluster groups
# exclude groups with too few cells
size <- table(rna$cluster_group_per_condition)
rna_subset <- subset(rna, subset = cluster_group_per_condition %in% names(size)[which(size>20)])
de_res <- presto::wilcoxauc(rna_subset, group_by="cluster_group_per_condition", seurat_assay="RNA")
de_res$pct_diff <- de_res$pct_in - de_res$pct_out
sig_res <- de_res %>% filter(padj<0.05 & pct_in>10 & logFC>0.1)
table(sig_res$group)
top_res <- sig_res %>% group_by(group) %>% top_n(200, wt=logFC)
deg_res <- list("top_res"=top_res, "sig_res"=sig_res)
saveRDS(deg_res, file="Res_Salmonella_infection_on_colon_epi_cluster_group_per_condition_marker_gene.rds")

# get the marker gene overlap with feature gene sets
# marker genes overlap with gene sets
gene_list <- readRDS("Res_KEGG_and_HALLMARK_gene_list.rds")
sig_res <- readRDS("Res_Salmonella_infection_on_colon_epi_goblet_and_s2E_DEG_sig_res.rds")
expressed_prop <- getExpressedProp(seu.obj=rna, feature.to.calc="Cluster_per_condition", size.cutoff=20)
all_expressed_genes <- rownames(expressed_prop)[apply(expressed_prop, 1, max)>0.1]
saveRDS(all_expressed_genes, file="Dat_s2E_and_goblet_all_expressed_genes.rds")
all_overlap_genes <- intersect(unlist(gene_list), all_expressed_genes)

# check enrichment on colonocyte infection upregulated genes
idx <- which(sig_res$condition=="Salmonella" & grepl("Colonocyte", sig_res$cell_group))
col_degs <- unique(sig_res$feature[idx])
idx <- which(sig_res$condition=="Salmonella" & sig_res$cell_group=="Goblet")
gob_degs <- unique(sig_res$feature[idx])
degs <- intersect(col_degs, gob_degs)

union_pval <- sapply(names(gene_list), function(gs){
        gene_set <- gene_list[[gs]]
        marker_genes <- degs
        overlap_genes <- intersect(gene_set, marker_genes)
        marker_non_gs_genes <- setdiff(marker_genes, gene_set)
        gs_non_marker_genes <- setdiff(gene_set, marker_genes)
        non_gs_non_marker_genes <- setdiff(all_overlap_genes, union(gene_set, marker_genes))
        a <- length(overlap_genes)
        b <- length(marker_non_gs_genes)
        c <- length(gs_non_marker_genes)
        d <- length(non_gs_non_marker_genes)
        fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="g")$p.value
    })

# check enrichment on infection induced DEGs on stem cell, colonocyte and goblet cell
pval <- sapply(sort(unique(sig_res$cell_group_per_condition)), function(group){
    sapply(names(gene_list), function(gs){
        gene_set <- gene_list[[gs]]
        marker_genes <- sig_res$feature[sig_res$cell_group_per_condition==group]
        overlap_genes <- intersect(gene_set, marker_genes)
        marker_non_gs_genes <- setdiff(marker_genes, gene_set)
        gs_non_marker_genes <- setdiff(gene_set, marker_genes)
        non_gs_non_marker_genes <- setdiff(all_overlap_genes, union(gene_set, marker_genes))
        a <- length(overlap_genes)
        b <- length(marker_non_gs_genes)
        c <- length(gs_non_marker_genes)
        d <- length(non_gs_non_marker_genes)
        fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="g")$p.value
    })
})
saveRDS(pval, file="Res_cell_group_DEG_and_gene_set_overlap_enrichment_pval.rds")

or <- sapply(sort(unique(sig_res$cell_group_per_condition)), function(group){
    sapply(names(gene_list), function(gs){
        gene_set <- gene_list[[gs]]
        marker_genes <- sig_res$feature[sig_res$cell_group_per_condition==group]
        overlap_genes <- intersect(gene_set, marker_genes)
        marker_non_gs_genes <- setdiff(marker_genes, gene_set)
        gs_non_marker_genes <- setdiff(gene_set, marker_genes)
        non_gs_non_marker_genes <- setdiff(all_overlap_genes, union(gene_set, marker_genes))
        a <- length(overlap_genes)
        b <- length(marker_non_gs_genes)
        c <- length(gs_non_marker_genes)
        d <- length(non_gs_non_marker_genes)
        fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="g")$estimate
    })
})
rownames(or) <- sub(".odds ratio", "", rownames(or))
saveRDS(or, file="Res_cell_group_DEG_and_gene_set_overlap_enrichment_odds_ratio.rds")

# get average module score per cell group per condition
rna_with_score <- readRDS("/home/yuq22/Bacteria_TRM/used_object/Dat_salmonella_colon_epithelium_scRNA-seq_data_with_subtypes_and_feature_scores.rds")
X <- t(rna_with_score@meta.data[,23:ncol(rna_with_score@meta.data)])
y <- rna_with_score$cluster_group_per_condition
score_data <- list(
    "X"=X,
    "y"=y
)
saveRDS(score_data, file="Dat_salmonella_data_KEGG_hallmark_gene_set_module_score_data.rds")
ave_score <- getAveExpr(input.type="matrix", X=X, y=y, size.cutoff=20, colname.prefix=NULL)
saveRDS(ave_score, file="Dat_salmonella_data_KEGG_hallmark_gene_set_module_score_per_cell_type_per_condition.rds")

rownames(ave_score) <- sub("hallmark ", "", gsub("_", " ", tolower(rownames(ave_score))))
saveRDS(ave_score, file="Dat_salmonella_data_KEGG_hallmark_gene_set_module_score_per_cell_type_per_condition_lowcase_term_name.rds")


"Antigen processing and presentation"
ave_score <- readRDS("Dat_salmonella_data_KEGG_hallmark_gene_set_module_score_per_cell_type_per_condition.rds")
pval <- readRDS("Res_cell_group_DEG_and_gene_set_overlap_enrichment_pval.rds")
#padj <- apply(pval, 2, function(x) p.adjust(x, method="BH"))
term_idx <- rownames(pval)[which(rowSums(pval[,grep("Salmonella",colnames(pval))]<0.05)>0)]
all_sig_term <- -log10(pval[term_idx,grep("Salmonella",colnames(pval))])
saveRDS(all_sig_term, file="Dat_all_sig_terms_for_Salmonella_upregulated_genes.rds")
# exclude disease related terms
#term_idx <- term_idx[!grepl("disease", term_idx)]
length(term_idx)
  

term_idx <- c(
"HALLMARK_TNFA_SIGNALING_VIA_NFKB"                                                                                                       
,"HALLMARK_HYPOXIA"
,"HALLMARK_APOPTOSIS" 
,"HIF-1 signaling pathway"  
,"Fructose and mannose metabolism" 
,"Galactose metabolism"                                             
,"HALLMARK_IL2_STAT5_SIGNALING"                                                   
,"HALLMARK_TGF_BETA_SIGNALING"  
"HALLMARK_UNFOLDED_PROTEIN_RESPONSE"                                                                                                                                                                                                               
,"Viral protein interaction with cytokine and cytokine receptor"                                            
,"Adherens junction"                                            
,"Tight junction"                                               
,"Leukocyte transendothelial migration"                                                                
,"Pathogenic Escherichia coli infection"                        
,"Shigellosis"                                                                                             
,"Bacterial invasion of epithelial cells"  
,"Amoebiasis"                                                   
,"Toxoplasmosis"                                                                                                                          
)  

input <- -log10(pval[term_idx,grep("Salmonella",colnames(pval))])
input_or <- or[term_idx,grep("Salmonella",colnames(pval))]
df <- data.frame('gene_set'=rep(rownames(input), ncol(input)),
                 'cell_group'=rep(colnames(input), each=nrow(input)),
                 'pval'=as.vector(input),
                 'odd_ratio'=as.vector(input_or),
                 stringsAsFactors = F)

library(ggplot2)
# order the condition
ubiquitous_terms <- c(
"HALLMARK_TNFA_SIGNALING_VIA_NFKB"                                                                                                       
,"HALLMARK_HYPOXIA"
,"HALLMARK_APOPTOSIS" 
,"HIF-1 signaling pathway"  
,"Fructose and mannose metabolism" 
,"Galactose metabolism"                                             
,"HALLMARK_IL2_STAT5_SIGNALING"                                                   
,"HALLMARK_TGF_BETA_SIGNALING"  
)
term_idx <- c(
"HALLMARK_UNFOLDED_PROTEIN_RESPONSE"                                                                                                                                                                                                               
,"Viral protein interaction with cytokine and cytokine receptor"                                            
,"Adherens junction"                                            
,"Tight junction"                                               
,"Leukocyte transendothelial migration"                                                                
,"Pathogenic Escherichia coli infection"                        
,"Shigellosis"                                                                                             
,"Bacterial invasion of epithelial cells"  
,"Amoebiasis"                                                   
,"Toxoplasmosis"                                                                                                                          
)  
max_ct <- colnames(input)[apply(input[term_idx, ], 1, which.min)]
term_order <- unlist(lapply(colnames(input), function(x){
    which(max_ct==x)
}))
term_order <- c(ubiquitous_terms, term_idx[term_order])
term_order_lower <- sub("hallmark ", "", gsub("_", " ", tolower(term_order)))
cell_type_order <- sub("_Salmonella", "", colnames(input))
df$gene_set <- sub("hallmark ", "", gsub("_", " ", tolower(df$gene_set)))
df$gene_set <- factor(df$gene_set, levels=rev(term_order_lower))
df$cell_group <- sub("_Salmonella", "", df$cell_group)
df$cell_group <- factor(df$cell_group, levels=cell_type_order)
saveRDS(df, file="Dat_cell_group_DEG_and_gene_set_overlap_enrichment.rds")

# generate barplot to show infection upregulated genes enriched terms
df <- readRDS("Dat_cell_group_DEG_and_gene_set_overlap_enrichment.rds")
enrichment <- sapply(sort(unique(df$gene_set)), function(x){
    max(df$pval[which(df$gene_set==x)])
})
names(enrichment) <- sort(unique(df$gene_set))

pdf("Plot_barplot_selected_term_enrichment.pdf", height=5)
barplot(sqrt(sort(enrichment)), horiz=T, las=2, xlab="-log10(p-value)")
abline(v=sqrt(-log10(0.05)), lty=2)
dev.off()


###
# get infection-downregulated gene enriched terms
setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/more_filterings/more_DEG_characterization/")
pval <- readRDS("Res_cell_group_DEG_and_gene_set_overlap_enrichment_pval.rds")
down_term_idx <- rownames(pval)[which(rowSums(pval[,grep("Control",colnames(pval))]<0.05)>0)]
up_term_idx <- rownames(pval)[which(rowSums(pval[,grep("Salmonella",colnames(pval))]<0.05)>0)]
# exclude disease related terms
term_idx <- down_term_idx[!grepl("disease", down_term_idx)]
length(term_idx)

term_idx <- c("HALLMARK_ADIPOGENESIS"                   
,"HALLMARK_CHOLESTEROL_HOMEOSTASIS"        
,"HALLMARK_FATTY_ACID_METABOLISM"          
,"HALLMARK_MTORC1_SIGNALING"               
,"HALLMARK_MYC_TARGETS_V1"                 
,"HALLMARK_OXIDATIVE_PHOSPHORYLATION"      
,"HALLMARK_PEROXISOME"                     
,"HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"                       
,"Propanoate metabolism"                   
,"Butanoate metabolism"                    
,"Oxidative phosphorylation"               
,"Sulfur metabolism"                       
,"Fatty acid biosynthesis"                 
,"Steroid biosynthesis"                    
,"Porphyrin metabolism"                    
,"Terpenoid backbone biosynthesis"         
,"Sulfur relay system"                     
,"Antigen processing and presentation"     
,"Fc gamma R-mediated phagocytosis"        
,"PPAR signaling pathway")     


# generate barplot to show infection upregulated genes enriched terms
enrichment <- apply(-log10(pval[term_idx, grep("Control",colnames(pval))]), 1, max)

pdf("Plot_barplot_selected_infection_down_regulated_gene_enriched_terms.pdf", height=10)
barplot(sqrt(sort(enrichment)), horiz=T, las=2, xlab="-log10(p-value)")
abline(v=sqrt(-log10(0.05)), lty=2)
dev.off()


###



# get module score difference induced upon salmonella infection
score_diff <- sapply(cell_type_order, function(x){
    if(grepl("infection_specific", x)){
        diff <- ave_score[,paste0(x,"_Salmonella")]
    }else{
        diff <- ave_score[,paste0(x,"_Salmonella")] - ave_score[,paste0(x,"_Control")]
    }    
    return(diff)
    
}) 
saveRDS(score_diff, file="Dat_salmonella_data_KEGG_hallmark_gene_set_module_score_diff_per_cell_type.rds")

diff_vec <- sapply(seq(nrow(df)), function(i){
    gs <- df$gene_set[i]
    ct <- df$cell_group[i]
    score_diff[gs, ct]
})
df$module_score_diff <- diff_vec
saveRDS(df, file="Dat_cell_group_DEG_and_gene_set_overlap_enrichment.rds")


p1 <- ggplot(df, aes(x=cell_group, y=gene_set, size=pval, fill=odd_ratio))+
    geom_point(shape=21, col="#202020")+
    scale_size(range=c(1,10))+
    scale_fill_gradientn(colors=beach.col)+
    theme_light()+
    guides(x =  guide_axis(angle = 45))+
    labs(x="Cell group", y="Gene set", fill="Odd ratio", size="-log10(p-value)")
p1
ggsave(p1, filename="Plot_dot_selected_term_enrichment_2.pdf", width=8, height=7)

p2 <- ggplot(df, aes(x=cell_group, y=gene_set, size=module_score_diff, fill=pval))+
    geom_point(shape=21, col="#202020")+
    scale_size(range=c(1,10))+
    scale_fill_gradientn(colors=beach.col)+
    theme_light()+
    guides(x =  guide_axis(angle = 45))+
    labs(x="Cell group", y="Gene set", size="Odd ratio", fill="module_score_diff")
p2
ggsave(p2, filename="Plot_dot_selected_term_enrichment_3.pdf", width=8, height=7)


# plot mit%, gene number and transcript count per RNA cluster in the salmonella infection condition
sal <- subset(rna, subset=condition=="Salmonella")
p1 <- SCpubr::do_ViolinPlot(sal, group.by="cluster_group", feature=c("nCount_RNA", "nFeature_RNA", "percent.mt"), colors.use=cols$cell_type)&NoLegend()
ggsave(p1, filename="Plot_violin_salmonella_infection_only_sample_RNA_cluster_gene_number_and_mit_percent.pdf", width=20, height=5)

# load seurat object
rna <- readRDS("/home/yuq22/Bacteria_TRM/used_object/Dat_salmonella_colon_epithelium_scRNA-seq_data_with_subtypes_and_feature_scores.rds")
term <- "Antigen processing and presentation"
term <- c("NF-kappa B signaling pathway", "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
term <- c("Adherens junction", "Tight junction")


sig_res <- readRDS("Res_Salmonella_infection_on_colon_epi_goblet_and_s2E_DEG_sig_res.rds")
union_degs <- unique(sig_res$feature)
features <- intersect(gene_list[[term]], union_degs)
top_res <- sig_res %>% group_by(cell_group_per_condition) %>% top_n(50, wt=logFC)
union_top_degs <- unique(top_res$feature)
features <- intersect(unlist(gene_list[term]), union_top_degs)

# subset to goblet cells and s2E cells
# exclude groups with too few cells
size <- table(rna$cluster_group_per_condition)
selected_pops <- grep("SC|Colonocyte|Goblet", names(size)[which(size>20)], value=T)
rna_subset <- subset(rna, subset = cluster_group_per_condition %in% selected_pops)
saveRDS(rna_subset, file="/home/yuq22/Bacteria_TRM/used_object/Dat_salmonella_colon_epithelium_scRNA-seq_data_GC_SC_colonocyte.rds")
rna_subset <- readRDS("/home/yuq22/Bacteria_TRM/used_object/Dat_salmonella_colon_epithelium_scRNA-seq_data_GC_SC_colonocyte.rds")
p1 <- SCpubr::do_DotPlot(rna_subset, features=features, scale=T, group.by="cluster_group_per_condition")
ggsave(p1, filename="Plot_dotplot_junction_gene_expression.pdf", width=10, height=10)

feature_list <- list(
    "Antigen processing and presentation"=c("TNF", "CD74", "HLA-A",  "B2M","HLA-E", "TAPBP"),
    "NFKB"=c("LTB", "TNFAIP2", "BIRC3", "NFKBIA", "NAMPT", "SQSTM1"),
    "Junction"=c("NECTIN3", "PRKAG2", "LMO7", "ACTN4","TCF7L2","ACTB")
)
p1 <- SCpubr::do_DotPlot(rna_subset, features=feature_list, scale=T, group.by="cluster_group_per_condition")
ggsave(p1, filename="Plot_dotplot_selected_gene_expression_combined.pdf", width=10, height=10)




# subset to colonocyte 4, goblet cells and SC cells
#rna_subset <- subset(rna, subset = cluster_group %in% c("Colonocyte_4", "Goblet", "Colon_SC"))
#saveRDS(rna_subset, file="/home/yuq22/Bacteria_TRM/used_object/Dat_salmonella_colon_epithelium_scRNA-seq_data_GC_SC_colonocyte4.rds")
#p1 <- SCpubr::do_DotPlot(rna_subset, features=c("TNF", "CD74", "HLA-A",  "B2M","HLA-E", "TAPBP"), scale=T, group.by="cluster_group_per_condition")
#ggsave(p1, filename="Plot_dotplot_antigen_presenting_gene_expression_selected.pdf", width=5, height=5)

