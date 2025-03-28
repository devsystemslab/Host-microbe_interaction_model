setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/GRN_on_more_filtered_data/incoorperate_multiple_GRN/topDE_GRN")

library(Pando)
data(motif2tf)
pando_tf <- unique(motif2tf$tf)
# load gene module
module_combined <- readRDS("../Res_combined_GRN_modules_no_all_epi_cell_type_module.rds")

# load infection vs control DEG test result
combined_de_res <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/used_object/Res_Salmonella_infection_on_colon_epi_goblet_and_s2E_all_DEG_res.rds")
# get the fold change of infection - control
all_fc_data <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/NFKB_barrier_innate_sensing/combined_with_Ruben_gene_set/Dat_expr_logFC_infection_control.rds")
all_pval_data <- readRDS("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/NFKB_barrier_innate_sensing/combined_with_Ruben_gene_set/Dat_expr_pval_infection_control.rds")

# for each expressed genes, get the minimum p-value and the maximum absolute fold change
min_pval <- setNames(apply(all_pval_data, 1, min),rownames(all_pval_data))
max_abs_fc <- setNames(apply(all_fc_data, 1, function(vec){vec[which.max(abs(vec))]}), rownames(all_fc_data))
# get the cluster, then cell type with the maximum absolute fold change
max_abs_fc_cl <- setNames(apply(all_fc_data, 1, function(vec){colnames(all_fc_data)[which.max(abs(vec))]}), rownames(all_fc_data))
max_abs_fc_ct <- sapply(max_abs_fc_cl, function(x){

    vec <- strsplit(split="_", x)[[1]]
    ct <- paste(vec[-c((length(vec)-1):length(vec))], collapse="_")
    return(ct)
    
})
df_gene <- data.frame(
    "gene"=rownames(all_fc_data),
    "max_fc"=max_abs_fc,
    "abs_max_fc"=abs(max_abs_fc),
    "min_pval"=min_pval,
    "logp"=-log10(min_pval),
    "max_fc_cl"=max_abs_fc_cl,
    "max_fc_ct"=max_abs_fc_ct,
    "direction"=ifelse(max_abs_fc>0, 1, -1),
    "tf"=ifelse(rownames(all_fc_data) %in% pando_tf, 1, 0),
    stringsAsFactors = F
)
# group the cell type to a more coarse grain
ct <- c("Colonocyte_1", "Colonocyte_2", "Colonocyte_3", "Colonocyte_4", 
        "Colonocyte_2_infection_specific_C14", "Colonocyte_2_infection_specific_C17", "Colonocyte_2_infection_specific_C16", 
        "Colon_SC", "Goblet")
coarse_ct <- rep(c("Non_infection_specific_colonocyte", "Infection_specific_colonocyte_subtype_1", 
                    "Infection_specific_colonocyte_subtype_2", "Colon_SC", "Goblet"),
                    c(4,2,1,1,1))                        
names(coarse_ct) <- ct
df_gene$coarse_ct <- coarse_ct[df_gene$max_fc_ct]

# get the significant DEGs
combined_sig_res <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/used_object/Res_Salmonella_infection_on_colon_epi_goblet_and_s2E_DEG_sig_res.rds")
union_degs <- unique(combined_sig_res$feature)
# indicate whether the gene is a significant DEG
df_gene$deg <- df_gene$gene%in%union_degs
saveRDS(df_gene, file="Res_salmonella_infection_DEG_fc_pval_cluster_cell_type.rds")

# select top DEGs with the largest change magnitude
top_deg <- df_gene$gene[order(df_gene$abs_max_fc, decreasing=T)[1:500]]
intersect(top_deg, c("PLA2G2A", "DMBT1"))
# get the TFs among the top DEGs
top_deg_tf <- intersect(top_deg, pando_tf)
df_gene$top_deg <- df_gene$gene%in%top_deg
saveRDS(df_gene, file="Res_salmonella_infection_DEG_fc_pval_cluster_cell_type.rds")

# get the TF-target pairs in the significant modules where the target is in the top DEGs
sig_module <- unique(module_combined[which(module_combined$padj<0.05 & module_combined$target%in%top_deg), c("tf", "target")])
# load cluster average expression
ave_expr <- readRDS("../Data_ave_expr_per_cluster_group_per_condition.rds")
# calculate the correlation between the TFs and the target genes
cor_vec <- sapply(seq(nrow(sig_module)), function(i){
    tf <- sig_module$tf[i]
    target <- sig_module$target[i]
    cor(ave_expr[tf,], ave_expr[target,], method="pearson")
})
sig_module$cor <- cor_vec
# for each target gene, get the top N TF with the strongest correlation
N <- 3
tf_num_vec <- table(sig_module$target)
sig_module$input <- rep(F, nrow(sig_module))
# for target genes with at most N TFs, select all of them
sig_module[sig_module$target %in% names(tf_num_vec)[which(tf_num_vec<=N)], "input"] <- T 
# for target genes with more than N TFs, select the top N TFs with the strongest correlation
popular_target <- names(tf_num_vec)[which(tf_num_vec>N)]
for(g in popular_target){
    tf_idx <- which(sig_module$target==g)
    top_tf_idx <- tf_idx[order(abs(sig_module$cor[tf_idx]), decreasing=T)[1:N]]
    sig_module[top_tf_idx, "input"] <- T
}


# check whether all the top DE TFs are included in the input module
top_deg_tf_not_in_input <- setdiff(top_deg_tf, input_module$tf)
length(top_deg_tf_not_in_input) # 1
# get the most correlated target genes for the TF
module1 <- unique(module_combined[which(module_combined$padj<0.05 & module_combined$tf==top_deg_tf_not_in_input), c("tf", "target")])
cor_vec <- cor(t(ave_expr[module1$target,]), ave_expr[top_deg_tf_not_in_input,], method="pearson")[,1]
module1$cor <- cor_vec
module1$input <- F
module1$input[order(abs(module1$cor), decreasing=T)[1:N]] <- T

input_module <- rbind(sig_module[sig_module$input,], module1[module1$input,])
saveRDS(input_module, file="Res_salmonella_infection_GRN_input_module.rds")

input_gene <- union(input_module$tf, input_module$target)
length(input_gene)
# construct co-expression network of the genes
pcc_mat <- cor(t(ave_expr[input_gene,]), method="pearson")
input <- pcc_mat
pca_mat <- irlba::prcomp_irlba(input, n=20)$x
rownames(pca_mat) <- rownames(input)
umap_layout <- uwot::umap(as.matrix(pca_mat))
saveRDS(umap_layout, file="Res_salmonella_infection_GRN_umap_layout.rds")
nrow(umap_layout)


umap_layout <- readRDS("Res_salmonella_infection_GRN_umap_layout.rds")
umap_genes <- intersect(rownames(umap_layout), rownames(df_gene))
length(umap_genes)
node_df <- df_gene[umap_genes,]
node_df$X <- umap_layout[umap_genes,1]
node_df$Y <- umap_layout[umap_genes,2]
node_df$logp_transformed <- node_df$logp^0.25
# highlight the TFs with the largest logFC
tf_df <- node_df[node_df$tf==1,]
#tf_highlight <- tf_df[order(abs(tf_df$max_fc), decreasing=T)[1:20], "gene"]
tf_highlight <- tf_df[order(abs(tf_df$max_fc), decreasing=T)[1:35], "gene"]
node_df$highlight <- ifelse(node_df$gene %in% tf_highlight, 1, 0)
saveRDS(node_df, file="Res_salmonella_infection_GRN_coex_network_node_df.rds")

node_df <- readRDS("Res_salmonella_infection_GRN_coex_network_node_df.rds")

# get edge list
edge_df <- input_module[which(input_module$tf %in% node_df$gene & input_module$target %in% node_df$gene), ]
edge_df$direction <- ifelse(edge_df$cor>0, 1, -1)
saveRDS(edge_df, file="Res_salmonella_infection_GRN_coex_network_edge_df.rds")

edge_df <- readRDS("Res_salmonella_infection_GRN_coex_network_edge_df.rds")
cols <- readRDS("/projects/site/pred/ihb-g-deco/USERS/yuq22/Bacteria_TRM/used_object/cols.rds")
ct_cols <- cols$cell_type

# generate a graph to visualize the co-expression network
# color the node according to the cell type with maximal expression level
# Create a graph object using tidygraph
library(tidygraph)
library(ggraph)
# Create a graph object from the edge list
ggraph <- tbl_graph(nodes = node_df, edges = edge_df, directed = FALSE)
saveRDS(ggraph, file="Res_salmonella_infection_GRN_module_coex_graph.rds")

ggraph <- readRDS("Res_salmonella_infection_GRN_module_coex_graph.rds")
p1 <- ggraph(ggraph, x=X, y=Y) +
    # color the edges based on the direction of the correlation
    geom_edge_diagonal(aes(color=as.factor(direction)),
                     width=0.5,
                     alpha=0.1)+
    scale_edge_color_manual(values = c("1"="#67001F", "-1"="#053061")) +
    geom_node_point(aes(size = abs_max_fc, fill = coarse_ct),
                  shape = 21, color='darkgrey') +
    scale_fill_manual(values=ct_cols)+
    scale_size_continuous(range = c(1,10)) +
    geom_node_label(aes(label=gene, filter = node_df$highlight == 1, color=as.factor(direction)), repel=T, max.overlaps = 999) +
    scale_color_manual(values = c("1"="#67001F", "-1"="#053061")) +
    #labs(title=x) +
    theme_void()
p1  
ggsave(p1, filename="Plot_ggraph_regulome_coexpression_network_including_BACH1.pdf", height=10, width=10)

write.table(edge_df, file="Dat_inferred_TF_target_pairs.txt", sep="\t", row.names=F, quote=F)

setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/GRN_on_more_filtered_data/incoorperate_multiple_GRN/topDE_GRN/")
edge_df <- readRDS("Res_salmonella_infection_GRN_coex_network_edge_df.rds")
example <- edge_df[which(edge_df$tf%in%c("STAT3","BHLHE40")),]
example <- edge_df[which(edge_df$target%in%"PLA2G2A"),]
rownames(example) <- NULL

# show the expression of selected genes in colonocytes between conditions



