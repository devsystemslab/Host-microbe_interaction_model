setwd("/home/yuq22/Bacteria_TRM/TRM/update_epi_anno_and_cellChat")
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(ggrepel)
library(CellChat)
library(Seurat)
source("/projects/site/pred/ihb-intestine-evo/common_script/Script_functions.R")

# load epi seurat object
epi <- readRDS("/home/yuq22/Bacteria_TRM/TRM/Dat_merged_filtered_ileum_epithelium_TRM_coculture_noImmune_noLowQuality_clustered_updated_v2.rds")

# examine the feature number distribution per cluster
SCpubr::do_BoxPlot(epi, feature='nFeature_RNA', group.by="anno_v1")

# exclude TA-1, which has much fewer features than other clusters
epi <- subset(epi, anno_v1 != "TA-1")
SCpubr::do_DimPlot(epi, group.by="Coarse_cell_type")
saveRDS(epi, file="/home/yuq22/Bacteria_TRM/TRM/Dat_merged_filtered_ileum_epithelium_TRM_coculture_noImmune_noLowQuality_clustered_updated_v2.rds")

epi_co <- subset(epi, orig.ident=="withTRM")
epi_co$Cell_type <- epi_co$Coarse_cell_type
# load immune cell seurat object
trm <- readRDS("/home/yuq22/Bacteria_TRM/TRM/TRM_from_coculture/cell_type_annotation/coarse_cell_type_annotation/Dat_TRM_with_updated_cell_type_anno.rds")
trm$Cell_type <- trm$Updated_cell_type_index
combined <- merge(x=epi_co, y=trm)
saveRDS(combined, file="Dat_cocultured_epi_TRM_combined_seurat_object.rds")

# run cellchat
combined <- readRDS("Dat_cocultured_epi_TRM_combined_seurat_object.rds")
data.input <- combined[["RNA"]]@data # normalized data matrix
labels <- combined$Cell_type
meta <- data.frame(labels = labels, row.names = names(labels))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
groupSize <- as.numeric(table(cellchat@idents))
cellchat@DB  <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 10) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.05)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file = "Res_cellchat_cocultured_epi_TRM.rds")

cellchat <- readRDS("Res_cellchat_cocultured_epi_TRM.rds")
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

df.net <- subsetCommunication(cellchat)
"CDH1_ITGAE_ITGB7" # CDH1 - ITGAE [CD103]
"CD160_TNFRSF14" # CD160 - TNFRSF14 [HVEM]

df <- df.net[which(df.net$interaction_name=="CDH1_ITGAE_ITGB7" & df.net$target=="4:gdT"),]
combined <- readRDS("Dat_cocultured_epi_TRM_combined_seurat_object.rds")
# focus on interactions with enterocytes
immune_cell_types <- setdiff(unique(combined$Cell_type[which(combined$orig.ident=="TRM")]), "6:G2M CD4")
# does enterocytes tend to be sender or receiver? 
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
ggsave(gg1, file="Plot_overall_signaling_role_scatter.pdf", width=6, height=6)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2


# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat, sources.use = "Prox_SI_Ent", targets.use =immune_cell_types, lab.cex = 0.5,legend.pos.y = 30)

# get KEGG pathways
kegg_df <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Annotation/KEGG/Res_KEGG_gene_list_df.rds")
# get genes involved in antigen processing and presentation
# load manually curated interacting LR pairs
manual_interactions <- read.table("/home/yuq22/group_share/Annotation/TRM_epi_interaction_annotation/Table_manually_curated_interaction.txt", sep="\t", head=T)
lympho_genes <- sort(unique(unlist(strsplit(manual_interactions$Intraepithelial.Lymphocytes, split="/"))))
epi_genes <- sort(unique(unlist(strsplit(manual_interactions$Epithelium, split="/"))))
sort(unique(kegg_df$pathway[which(kegg_df$gene %in% lympho_genes)]))

kegg_pathway <- "Antigen processing and presentation"  
antigene_pathway_genes <- sort(unique(kegg_df$gene[which(kegg_df$pathway==kegg_pathway)])) # antigen processing and presentation - APP
app_interactions <- sort(unique(df.net$interaction_name[grep(paste(antigene_pathway_genes,collapse="|"), df.net$interaction_name)]))
res <- t(apply(manual_interactions, 1, function(vec) {
    g1 <- gsub("/", "|", vec[1])
    g2 <- gsub("/", "|", vec[2])
    idx <- which(grepl(g1, df.net$interaction_name) & grepl(g2, df.net$interaction_name))
    lr_pairs <- as.character(sort(unique(df.net$interaction_name[idx])))
    pathways <- as.character(sort(unique(df.net$pathway_name[idx])))
    if(length(lr_pairs)>0) {
        return(c(paste(lr_pairs, collapse=","), paste(pathways, collapse=",")))
    }else{
        return(rep(NA,2))
    }

    
}))
manual_interactions$LR_pairs <- res[,1]
manual_interactions$Pathways <- res[,2]
saveRDS(manual_interactions, file="Res_manually_curated_interaction.rds")

manual_interactions <- readRDS("Res_manually_curated_interaction.rds")
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
p1 <- netVisual_chord_gene(cellchat, sources.use = "Prox_SI_Ent", targets.use =immune_cell_types, pairLR.use = data.frame('interaction_name'=setdiff(unique(manual_interactions$LR_pairs), NA)) , legend.pos.x = 8)
pdf("Plot_manual_annotated_significant_interaction_chord_plot.pdf")
p1
dev.off()

p2 <- netVisual_bubble(cellchat, sources.use = "Prox_SI_Ent", targets.use = immune_cell_types, pairLR.use = data.frame('interaction_name'=setdiff(unique(manual_interactions$LR_pairs), NA)), remove.isolate = TRUE)
pdf("Plot_manual_annotated_significant_interaction_bubble_plot.pdf")
p2
dev.off()


library(ggplot2)
lr_pairs <- sort(unique(unlist(strsplit(setdiff(unique(manual_interactions$LR_pairs), NA), split=","))))
# enterocytes as sender
prob_data <- matrix(0, ncol=length(immune_cell_types), nrow=length(lr_pairs))
rownames(prob_data) <- lr_pairs
colnames(prob_data) <- immune_cell_types
for(j in colnames(prob_data)){
    idx <- which(df.net$source=="Prox_SI_Ent" & df.net$target==j)
    prob_vec <- setNames(df.net$prob[idx], df.net$interaction_name[idx])
    sig_annotated_pairs <- intersect(names(prob_vec), lr_pairs)
    prob_data[sig_annotated_pairs, j] <- prob_vec[sig_annotated_pairs]
}
sender_prob_data <- prob_data

# enterocytes as receiver
prob_data <- matrix(0, ncol=length(immune_cell_types), nrow=length(lr_pairs))
rownames(prob_data) <- lr_pairs
colnames(prob_data) <- immune_cell_types
for(j in colnames(prob_data)){
    idx <- which(df.net$target=="Prox_SI_Ent" & df.net$source==j)
    prob_vec <- setNames(df.net$prob[idx], df.net$interaction_name[idx])
    sig_annotated_pairs <- intersect(names(prob_vec), lr_pairs)
    prob_data[sig_annotated_pairs, j] <- prob_vec[sig_annotated_pairs]
}
receiver_prob_data <- prob_data
colnames(receiver_prob_data) <- paste("From", colnames(receiver_prob_data), sep=":") # enterocyte as receiver, immune cell as sender
colnames(sender_prob_data) <- paste("To", colnames(sender_prob_data), sep=":") # enterocyte as sender, immune cell as receiver
combined_prob_data <- cbind(sender_prob_data, receiver_prob_data)
idx <- apply(combined_prob_data, 1, function(vec){
    max(vec)>0.0025
})
input <- combined_prob_data[idx,]
df <- data.frame(
    'interaction'=rep(rownames(input), ncol(input)),
    'partner_cell_type'=rep(colnames(input), each=nrow(input)), 
    'prob'=as.vector(input),
    stringsAsFactors = F)

p1 <- ggplot(df, aes(x=partner_cell_type, y=interaction, size=prob, fill=prob))+
    geom_point(shape=21, col="#202020")+
    scale_size(range=c(1,10))+
    scale_fill_gradientn(colors=beach.col)+
    theme_light()+
    guides(x =  guide_axis(angle = 45))+
    labs(x="Partner cell type", y="LR pairs", size="Probability", fill="Probability")
p1
ggsave(p1, filename="Plot_dot_sig_manual_interaction_v2.pdf", width=8, height=5)

# overlap with epithelial DE info
cocluster_sig_res <- readRDS("/home/yuq22/group_share/Bacteria_TRM/TRM/epi_DEG/Res_withTRM_vs_noTRM_cells_DEG_sig_res.rds")

# de novo identify top interactions between enterocytes and each of immune cells 
# load cellchat result
cellchat <- readRDS("Res_cellchat_cocultured_epi_TRM.rds")
df_net <- subsetCommunication(cellchat)
# enterocytes as receiver or sender
df_ent <- df_net[which(df_net$source=="Prox_SI_Ent" | df_net$target=="Prox_SI_Ent"),]
ent_receiver <- df_net[which(df_net$target=="Prox_SI_Ent"),]
immune_ct <- setdiff(grep(":", unique(ent_receiver$source), fixed=T, value=T), "6:G2M CD4")
all_interactions <- unique(df_ent$interaction_name)
dat <- matrix(0, ncol=length(immune_ct)*2, nrow=length(all_interactions))
rownames(dat) <- all_interactions
colnames(dat) <- paste(rep(c("From", "To"), each=length(immune_ct)), rep(immune_ct, 2), sep=":")
for(direction in c("From", "To")){
    for(immune_cell in immune_ct){
        if(direction=="To"){
            idx <- which(df_ent$source=="Prox_SI_Ent" & df_ent$target==immune_cell)
        }else{
            idx <- which(df_ent$target=="Prox_SI_Ent" & df_ent$source==immune_cell)
        }
        dat[as.character(df_ent$interaction_name[idx]), paste(direction, immune_cell, sep=":")] <- df_ent$prob[idx]
    }
}
saveRDS(dat, file="Res_enterocyte_interaction_prob.rds")

raw_dat <- readRDS("Res_enterocyte_interaction_prob.rds")
# get the interaction enrichment in particular cell type by subtracting the mean probability across all groups
dat <- t(scale(t(raw_dat)))

terms <- c("CCL25_CCR9", "CD160_TNFRSF14", "CD6_ALCAM", "CDH1_ITGAE_ITGB7", "HLA-E_KLRC1", "HLA-E_KLRK1", "IFNG_IFNGR1_IFNGR2")
colnames(dat)[apply(dat[terms, ], 1, which.max)]
# get top interactions from each column
idx <- apply(dat, 2, function(vec){
    order(vec, decreasing=T)[1:20]
})
union_idx <- union(rownames(dat)[sort(unique(as.vector(idx)))], terms)
input <- raw_dat[union_idx,]

res <- prepareTreeAndHeatmapInput(expr.mat=input, column.reorder=F)
saveRDS(res, file="Res_enterocyte_interaction_heatmap_input.rds")
pdf("Plot_heatmap_enterocyte_interaction_scaled_probability.pdf", height=10, width=8)
gplots::heatmap.2(res$heatmap_input, col = beach.col.heatmap, trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 0.1, cexCol=1, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
dev.off()

res <- readRDS("Res_enterocyte_interaction_heatmap_input.rds")
res$heatmap_input -> input
write.table(input, file="Table_denovo_top20_enterocyte_with_TRM_interaction_heatmap_input.txt", sep="\t", quote=F, col.names=T, row.names=T)


# update the analysis with the new annotation
setwd("/home/yuq22/Bacteria_TRM/TRM/update_epi_anno_and_cellChat")
# load data
cellchat <- readRDS("Res_cellchat_cocultured_epi_TRM.rds")
df.net <- subsetCommunication(cellchat)

# get the cellchat DB human database
all_interaction  <- as.character(sort(unique(CellChatDB.human$interaction$interaction_name)))
all_interaction <- gsub(":", "_", all_interaction)
# load manually curated interacting LR pairs
manual_interactions <- read.table("/home/yuq22/group_share/Annotation/TRM_epi_interaction_annotation/Table_manually_curated_interaction.txt", head=T)

interaction_component <- lapply(seq(length(all_interaction)), function(i){
    strsplit(all_interaction[i], split="_")[[1]]
})
names(interaction_component) <- all_interaction
component_genes <- sort(unique(unlist(interaction_component)))
component_idx <- matrix(0, nrow=length(component_genes), ncol=length(all_interaction))
rownames(component_idx) <- component_genes
colnames(component_idx) <- all_interaction
for(x in names(interaction_component)){
    component_idx[interaction_component[[x]], x] <- 1
}
saveRDS(component_idx, file="Res_interaction_component_idx.rds")

component_idx <- readRDS("Res_interaction_component_idx.rds")
db_name_vec <- c()
for(i in seq(nrow(manual_interactions))){
    # for genes1
    # get the gene symbol from the annotation
    genes1 <- strsplit(split="/", manual_interactions$Intraepithelial_Lymphocytes[i])[[1]]
    # get the annotated gene that are among the detected interaction component
    genes1 <- intersect(genes1, rownames(component_idx))
    
    # for genes2
    genes2 <- strsplit(split="/", manual_interactions$Epithelium[i])[[1]] 
    genes2 <- intersect(genes2, rownames(component_idx))
    
    if(length(genes1)==0 | length(genes2)==0){
        db_name <- NA
    }else{
        if(length(genes1)>1){
            names1 <- colnames(component_idx)[colSums(component_idx[genes1, ])>0]
        }else if(length(genes1)==1){
            names1 <- colnames(component_idx)[which(component_idx[genes1, ]>0)]
        }

        if(length(genes2)>1){
            names2 <- colnames(component_idx)[colSums(component_idx[genes2, ])>0]
        }else if(length(genes2)==1){
            names2 <- colnames(component_idx)[which(component_idx[genes2, ]>0)]
        }

        shared_interaction_name <- intersect(names1, names2)
        if(length(shared_interaction_name)==0){
            db_name <- NA
        }else{
            db_name <- paste(shared_interaction_name, collapse=",")
        }
    }
    db_name_vec <- c(db_name_vec, db_name)
}
manual_interactions$db_name <- db_name_vec
saveRDS(manual_interactions, file="Res_manually_curated_interaction_with_cellChatDB_interaction_name.rds")
setwd("/home/yuq22/Bacteria_TRM/TRM/update_epi_anno_and_cellChat/")
# load manually curated interacting LR pairs
manual_interactions <- readRDS("Res_manually_curated_interaction_with_cellChatDB_interaction_name.rds")
# check whether the epithelial components are DEGs in the TRM vs noTRM comparison
cocluster_sig_res <- readRDS("/home/yuq22/group_share/Bacteria_TRM/TRM/epi_DEG/Res_withTRM_vs_noTRM_cells_DEG_sig_res.rds")
ent_deg <- unique(cocluster_sig_res$feature[grepl("E-", cocluster_sig_res$group)])
de_idx <- sapply(manual_interactions$Epithelium, function(x){
    length(intersect(strsplit(x, split="/")[[1]], ent_deg))>0
})
manual_interactions$Epi_DEG <- de_idx
manual_interactions$plot_name <- paste0("[", manual_interactions$Intraepithelial_Lymphocytes, ", ", manual_interactions$Epithelium, "]")
saveRDS(manual_interactions, file="Res_manually_curated_interaction_with_cellChatDB_interaction_name.rds")


# get the interactions where the epithelial components are DEGs
epi_de_interactions <- unique(unlist(strsplit(split=",", manual_interactions$db_name[which(manual_interactions$Epi_DEG)])))

# visualize the interaction probability
setwd("/home/yuq22/Bacteria_TRM/TRM/update_epi_anno_and_cellChat/")

cellchat <- readRDS("Res_cellchat_cocultured_epi_TRM.rds")
df.net <- subsetCommunication(cellchat)
combined <- readRDS("Dat_cocultured_epi_TRM_combined_seurat_object.rds")
# load the interaction probability between enterocytes and immune cells
raw_dat <- readRDS("Res_enterocyte_interaction_prob.rds")

# subse to the manually curated ones
# load manually curated interacting LR pairs
manual_interactions <- readRDS("Res_manually_curated_interaction_with_cellChatDB_interaction_name.rds")
# exclude the ones that are not in the cellchat database
manual_interactions <- manual_interactions[!is.na(manual_interactions$db_name),]

df <- unique(manual_interactions[,c("db_name", "plot_name")])
db_name_list <- lapply(seq(nrow(df)), function(i){
    strsplit(df$db_name[i], split=",")[[1]]
})
db_name_vec <- unlist(db_name_list)
plot_name_vec <- rep(df$plot_name, sapply(db_name_list, length))
df <- unique(data.frame('db_name'=db_name_vec, 'plot_name'=plot_name_vec, stringsAsFactors = F))

selected_interactions <- intersect(sort(unique(unlist(strsplit(manual_interactions$db_name, split=",")))), rownames(raw_dat))
input_interaction <-  selected_interactions[which(apply(raw_dat[selected_interactions,], 1, sd)>0)]
input <- raw_dat[input_interaction,]

res <- prepareTreeAndHeatmapInput(expr.mat=input, column.reorder=F)
saveRDS(res, file="Res_selected_enterocyte_interaction_heatmap_input_v3.rds")

interaction_names <- rownames(res$heatmap_input)
row_cols <- ifelse(interaction_names %in% epi_de_interactions, "#303030", "#d0d0d0")

used_plot_name <- sapply(interaction_names, function(x){
    paste(unique(df$plot_name[which(df$db_name==x)]), collapse=",")
})
rownames(res$heatmap_input) <- used_plot_name

pdf("Plot_heatmap_selected_enterocyte_interaction_scaled_probability_v4.pdf", height=5, width=5)
gplots::heatmap.2(res$heatmap_input, col = beach.col.heatmap, RowSideColors=row_cols, trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 1, cexCol=1, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
dev.off()

genes <- c("IFNGR1", "IFNGR2", "TNFRSF14", "ALCAM", "IL15", "MICA", "ULBP3", "CCL25", "CDH1", "HLA-E", "ULBP3")
res <- setNames(genes%in%ent_deg, genes)
# check expression of BTNL3/BTNL8 in epithelium +/- TRM, TRGV4 / TRGV2 / TRGV9 in TRM
epi <- readRDS("/home/yuq22/Bacteria_TRM/TRM/Dat_merged_filtered_ileum_epithelium_TRM_coculture_noImmune_noLowQuality_clustered_updated_v2.rds")
trm <- readRDS("/home/yuq22/Bacteria_TRM/TRM/TRM_from_coculture/cell_type_annotation/coarse_cell_type_annotation/Dat_TRM_with_updated_cell_type_anno.rds")

plot_list <- list()
for(g in c("BTNL3", "BTNL8")){
    plot_list[[g]] <- SCpubr::do_BoxPlot(epi, feature=g, group.by="Cell_type_per_condition")+NoLegend()
}
p1 <- cowplot::plot_grid(plotlist=plot_list, ncol=1)
ggsave(p1, filename="Plot_boxplot_BTNL3_BTNL8_expression_in_epi.pdf", width=5, height=10)
genes <- c("TRGV4", "TRGV2", "TRGV9")
plotFeature(seu.obj=trm, genes.to.plot=genes, plot.name="Plot_UMAP_TRGVs_in_TRM.png", col.num=3, nCols=beach.col)

# generate heatmap to show expression of selected genes in epithelium + TRM
# load epi and TRM seurat object
trm$group <- trm$Updated_cell_type_index
epi$group <- epi$Cell_type_per_condition
combined <- merge(x=epi, y=trm)
ave_expr <- getAveExpr(seu.obj=combined, feature.to.calc="group", colname.prefix=NULL)
genes <- intersect(c("BTNL3", "BTNL8", "BTNL1", "BTNL6", "TRGV4", "TRGV2", "TRGV9"), rownames(ave_expr))
input <- ave_expr[genes,]
res <- prepareTreeAndHeatmapInput(expr.mat=input, column.reorder=F)
saveRDS(res, file="Res_selected_BTNLs_TRGVs_expr_heatmap_input_2.rds")

pdf("Plot_heatmap_BTNLs_TRGVs_expr_2.pdf", height=5, width=7)
gplots::heatmap.2(res$heatmap_input, col = beach.col.heatmap, trace="none", scale="none", density.info = "none", margins = c(15,10), cexRow = 1, cexCol=1, Rowv = FALSE, Colv = FALSE, dendrogram = "none")
dev.off()

SCpubr::do_BoxPlot(epi, feature="BHLHE40", group.by="Cell_type_per_condition")


