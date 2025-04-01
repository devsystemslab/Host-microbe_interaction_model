coverage_plot_with_evo_sig <- function(seu_obj, 
                                       peak_assay="peaks_merged", 
                                       rna_assay="RNA", 
                                       group="Coarse_cell_type_per_condition", 
                                       colors.use=NULL, 
                                       peaks, 
                                       evo_list_path="/projects/site/pred/ihb-intestine-evo/evo_signature/summary_Feb_2024/Dat_SNC_GWAS_HAR_PS_HAQER.rds", 
                                       peak_anno_path=NULL, #"Data_frame_peak_gene_pairs.rds"
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
                                       plot_suffix=""){

        require(Seurat)
        require(Signac)
        require(IRanges)
        require(GenomicRanges)
        require(ggplot2)
        require(cowplot)
        require(grid)
        require(gridExtra)

        # read evolutionary signature list
        evo_list <- readRDS(evo_list_path)
        # read peak annotation
        if(!is.null(peak_anno_path)){
          peak_anno <- readRDS(peak_anno_path)
        }

        DefaultAssay(seu_obj) <- peak_assay
        composite_plot_list <- list()
        for(p_human in peaks){
            print(p_human)
            if(gene_mode=="Specified"){
              g <- g_vec[p_human]
            }else if(gene_mode=="ChIPseeker"){
              g <- peak_anno[which(peak_anno[,peak_anno_peak_col]==p_human),peak_anno_gene_col]
            }

            p_list <- list()
            # 1. coverage data in hg38 coordinate
            cov_plot <- CoveragePlot(
              object = seu_obj,
              assay=peak_assay,
              group.by = group,
              region = p_human,
              extend.upstream = extend_up,
              extend.downstream = extend_down,
              region.highlight = StringToGRanges(p_human),
              annotation = FALSE,
              peaks = FALSE,
              tile = FALSE,
              links = FALSE,
              window = window_size
            )
            if(!is.null(colors.use)){
              cov_plot <- cov_plot + scale_fill_manual(values = colors.use)
            } 
            p_list[["coverage"]] <- cov_plot
            print("get coverage plot")

            # 2. genomic signature
            p_human_vec <- strsplit(p_human, split="-")[[1]]
            region <- GRanges(
              seqnames=p_human_vec[1],
              ranges=IRanges(
                start=as.numeric(p_human_vec[2])-extend_up,
                end=as.numeric(p_human_vec[3])+extend_down
              )
            )

            ## plot the overlapped sites
            for(x in names(evo_list)){
              # subset to covered range
              peak.intersect <- subsetByOverlaps(x = evo_list[[x]], ranges = region)
              peak.df <- as.data.frame(x = peak.intersect)
              start.pos <- start(x = region)
              end.pos <- end(x = region)
              chromosome <- seqnames(x = region)

              if (nrow(x = peak.df) > 0) {
              
                peak.df$start[peak.df$start < start.pos] <- start.pos
                peak.df$end[peak.df$end > end.pos] <- end.pos

                peak.plot <- ggplot(
                  data = peak.df
                )
                if(x %in% c("SNC","CMS","GWAS","All_human_SNC") ){
                  peak.plot <- ggplot(
                    data = peak.df
                  ) +
                  geom_segment(aes(x = start, y = 0, xend = end, yend = 1),
                               size =  0.5,
                               color = color,
                               data = peak.df)
                }else if(x %in% c("HAQER","HAR")){
                  peak.plot <- ggplot(
                  data = peak.df
                ) +
                  geom_rect(
                  data = peak.df,
                  aes_string(
                    xmin = "start",
                    xmax = "end",
                    ymin = 0,
                    ymax = 1),
                  fill = color
                  )
                }

              } else {
                # no peaks present in region, make empty panel
                peak.plot <- ggplot(data = peak.df)
              }
              peak.plot <- peak.plot + theme_classic() +
                ylab(label = x) +
                theme(axis.ticks.y = element_blank(),
                      axis.text.y = element_blank()) +
                xlab(label = paste0(chromosome, " position (bp)")) +
                xlim(c(start.pos, end.pos))

              # remove legend, change color
              peak.plot <- peak.plot +
                scale_color_manual(values = color) +
                theme(legend.position = "none")

              p_list[[x]] <- peak.plot
              }
            print("get signature overlapping plot")

            # 3. human peaks within the region
            peak_plot <- PeakPlot(
              object = seu_obj,
              region = p_human,
              extend.upstream = 10000,
              extend.downstream = 10000
            )
            p_list[["peak"]] <- peak_plot
            print("get peak plot")

            # 4. human annotation within the region
            gene_plot <- AnnotationPlot(
              object = seu_obj,
              region = p_human,
              extend.upstream = 10000,
              extend.downstream = 10000
            )
            p_list[["gene"]] <- gene_plot
            print("get annotation plot")

            # expression plot
            if(sum(g%in%rownames(seu_obj@assays[[rna_assay]]@data))>0){
                expr_plot <- ExpressionPlot(
                  object = seu_obj,
                  features = g,
                  assay = rna_assay,
                  group.by=group,
                )
                if(!is.null(colors.use)){
                  expr_plot <- expr_plot + scale_fill_manual(values = colors.use)
                } 
                print("get expression plot")
            }else{
                expr_plot <- ggplot()
                print("ChIPseeker annotated genes not detected in the scRNA-seq data")
            }
            p_list[["expr"]] <- expr_plot

            plot_combined <- CombineTracks(
                plotlist = p_list[-length(p_list)],
                expression.plot=expr_plot,
                heights = c(4, rep(0.3,length(p_list)-3), 0.5),
                widths = c(10, 2)
              )
            if(do_plot_individual){
              plot_name <- paste0("Plot_coveragePlot_for_hg38_",p_human,"_",g,plot_suffix,".pdf")
              ggsave(plot_combined, filename=plot_name, width=10, height=15)
            }
            composite_plot_list[[p_human]] <- plot_combined

          }
          if(do_plot_combined){
            ncol=min(c(4, length(peaks)))
            nrow=ceiling(length(peaks)/ncol)
            p <- plot_grid(plotlist = composite_plot_list, align = "hv", ncol = ncol)
            if(is.null(combined_plot_name)){
              plot_name <- "Plot_coveragePlot_for_hg38_selected_region.pdf"
            }
            ggsave(p, filename=plot_name, width=10*ncol, height=15*nrow)
          }
          return(composite_plot_list)
}
