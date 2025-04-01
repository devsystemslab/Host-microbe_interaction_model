setwd("/home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/ATAC_on_more_filtered_data/DAR")
library(Seurat)
library(Signac)
library(ggplot2)

# load input data
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
x <- args[1]
input_list <- readRDS(paste0("Res_Salmonella_infection_vs_control_",x,"_input_list.rds"))
bi_mat <- input_list$bi_mat
test_vec <- input_list$test_vec
umi_vec <- input_list$umi_vec  

library(doParallel)
registerDoParallel(20)
# run test
dd_test_res <- foreach(k=seq(nrow(bi_mat)), .multicombine = T, .combine = 'rbind')%dopar%{
  e <- as.numeric(as.vector(bi_mat[k,]))
  m0 <- glm(e ~ umi_vec, family = binomial)
  m1 <- glm(e ~ umi_vec+test_vec, family = binomial)
  a0 <- anova(m0)
  a1 <- anova(m1)
  p_anova <- anova(m1,m0, test = "Chisq")$Pr[2]
  p_resi <- pf((a0[nrow(a0),"Resid. Dev"]/a0[nrow(a0),"Resid. Df"]) / (a1[nrow(a1),"Resid. Dev"]/a1[nrow(a1),"Resid. Df"]), 
               df1 = a0[nrow(a0),"Resid. Df"], df2 = a1[nrow(a1),"Resid. Df"], lower.tail = F)
  coef <- coef(m1)[length(coef(m1))]
  return(c(p_anova, p_resi ,coef))
}
stopImplicitCluster()
rownames(dd_test_res) <- rownames(bi_mat)
colnames(dd_test_res) <- c("p_ANOVA", "p_Resi", "Coef")
resi_p_adj <- p.adjust(dd_test_res[,2], method = "BH")
anova_p_adj <- p.adjust(dd_test_res[,1], method = "BH")
dd_test_res <- as.data.frame(dd_test_res)
expressed_prop <- sapply(sort(unique(test_vec)), function(x){
  idx <- which(test_vec==x)
  prop_vec <- rowMeans(bi_mat[,idx])
  return(prop_vec)
})
colnames(expressed_prop) <- sort(unique(test_vec))
expressed_prop <- as.data.frame(expressed_prop)
expressed_prop$Prop_diff <- expressed_prop$Salmonella-expressed_prop$Control
res_mat <- data.frame("Test_anova_P"=dd_test_res$p_ANOVA,
                      "Test_resi_P"=dd_test_res$p_Resi,
                      "Corrected_anova_P"=anova_p_adj,
                      "Corrected_resi_P"=resi_p_adj,
                      "Condition_coef"=dd_test_res$Coef,
                      "Infection_expressed_prop"=expressed_prop$Salmonella,
                      "Control_expressed_prop"=expressed_prop$Control,
                      "Infection_control_prop_diff"=expressed_prop$Prop_diff)
 
saveRDS(res_mat, file=paste0("Res_Samonella_infection_vs_control_",x,"_differential_detection_test_res.rds"))

 files <- list.files(pattern="Res_Samonella_infection_vs_control_.*_differential_detection_test_res.rds")
dar_list <- list()
for(file in files){
    ct <- sub("_differential_detection_test_res.rds", "", sub("Res_Samonella_infection_vs_control_", "", file))
    res <- readRDS(file)
    sal_up_region <- rownames(res)[res$Corrected_anova_P<0.05 & res$Infection_control_prop_diff>0.005] 
    sal_down_region <- rownames(res)[res$Corrected_anova_P<0.05 & res$Infection_control_prop_diff< -0.005] 
    dar_list[[paste0(ct,"_up")]] <- sal_up_region
    dar_list[[paste0(ct,"_down")]] <- sal_down_region
}
all_dar_regions <- unique(unlist(dar_list))
matrix(F, nrow=length(all_dar_regions), ncol=length(dar_list), dimnames=list(all_dar_regions, names(dar_list))) -> dar_matrix
for(ct in names(dar_list)){
    dar_matrix[dar_list[[ct]],ct] <- T
}
saveRDS(dar_matrix, file="Res_Salmonella_infection_vs_control_DAR_matrix.rds")

# plot the DAR numbers
size <- sqrt(colSums(dar_matrix))
size[grep("down", names(size))] <- -size[grep("down", names(size))]
ct <- sub("_up", "", sub("_down", "", names(size)))
cols <- readRDS("~/Bacteria_TRM/used_object/cols.rds")
ct_cols <- cols$cell_type
pdf("Plot_barplot_DAR_numbers_per_cell_type.pdf", width=5, height=7)
par(mar=c(10,5,3,3))
barplot(size[grep("up", names(size))], col=ct_cols[ct[grep("up", names(size))]], main="Salmonella infection induced DARs", las=2, ylim=c(min(size), max(size)), ylab="sqrt(#DAR)", names=ct[grep("up", names(size))])
barplot(size[grep("down", names(size))], col=ct_cols[ct[grep("down", names(size))]], las=2, add=T, names="", xaxt="n", yaxt="n")
dev.off()
 
