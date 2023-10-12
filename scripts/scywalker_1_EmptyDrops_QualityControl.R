

p <- add_argument(p, "--percent_mito", help = "The cut-off for percentage of mitochondrial reads per cell", type = "numeric", nargs = '1', default = "2.5")
p <- add_argument(p, "--input", help = "Path to directory containing raw files in 10X format", type = "character", nargs = '1')
p <- add_argument(p, "--workdir", help = "Path to working directory", type = "character", nargs = '+')

arg <- parse_args(p)

input_df <- arg$input
i <- input_df$samplename


# 1. Load files
file_path <- paste(input_df$location)
barcodes <- paste(file_path, "barcodes.tsv.gz", sep = "")
features <- paste(file_path, "features.tsv.gz", sep = "")
matrix <- paste(file_path, "matrix.mtx.gz", sep = "")
mat <- Matrix::readMM(matrix)
features <- read.delim(features, header=F)
barcodes <- read.delim(barcodes, header=F)
colnames(mat) <- barcodes[,1]
rownames(mat) <- features[,2]

  
# 2. EmptyDrops
mito.genes <- grep(pattern = "^MT-", x = rownames(mat), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", x = rownames(mat), value = TRUE)
qc <- c(mito.genes, RPS.genes)
keep <- !row.names(mat) %in% unlist(qc)

br.out <- barcodeRanks(mat)
df <- data.frame(Rank = br.out$rank, Total = br.out$total, row.names = rownames(br.out))
df$AboveInflection <- df$Total > metadata(br.out)$inflection
e.out <- emptyDrops(mat[keep, ])
is.cell <- e.out$FDR <= 0.001
is.cell[is.na(is.cell)] <- FALSE
df_is.cell <- as.data.frame(is.cell)
rownames(df_is.cell) <- rownames(e.out)
df_is.cell <- merge(df_is.cell, df, by = "row.names")
  
p1 <- ggplot(df_is.cell, aes(x = Rank, y = Total, color = is.cell)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Rank", y = "Total") +
    theme_classic() +
    geom_hline(yintercept = metadata(br.out)$knee, color = "dodgerblue", linetype = "dotted") +
    geom_hline(yintercept = metadata(br.out)$inflection, color = "forestgreen", linetype = "dotted")+
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), 
                     guide = guide_legend(title = "Above Inflection"))

write.table(df_is.cell, paste0(arg$workdir, "/emptyDrops_barcode_classification.tsv"), row.names = FALSE, sep = "\t", quote = F)
filt_mat <- mat[, is.cell]
rm(mat, barcodes, features, e.out, is.cell)
  
SOB <- CreateSeuratObject(filt_mat, min.cells = 3)
rm(filt_mat)

 # 3. Quality metrics check
SOB[["percent.mt"]] <- PercentageFeatureSet(SOB, pattern = "^MT-")
t <- textGrob("can be filled")
p2 <- VlnPlot(SOB, features="nCount_RNA") + theme_classic() + theme(legend.position = 'none', axis.title.x=element_blank()) + scale_fill_manual(values = alpha("red", 0.5))
p3 <- VlnPlot(SOB, features= "nFeature_RNA") + theme_classic() + theme(legend.position = 'none', axis.title.x=element_blank()) + scale_fill_manual(values = alpha("red", 0.5))
p4 <- VlnPlot(SOB, features= "percent.mt") + theme_classic() + theme(legend.position = 'none', axis.title.x=element_blank()) + scale_fill_manual(values = alpha("red", 0.5))
p2$layers[[2]]$aes_params$alpha = 0.2
p3$layers[[2]]$aes_params$alpha = 0.2
p4$layers[[2]]$aes_params$alpha = 0.2
p5 <- ggplot(data = SOB@meta.data) +
    geom_point(mapping = aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), position = 'jitter', alpha = 0.3) +
    scale_color_viridis_b(option="magma", limits = c(0, 13)) +
    theme_classic() + theme(legend.position = 'none')
SOB_filtered <- subset(SOB, subset =  percent.mt < arg$percent_mt)
p6 <- VlnPlot(SOB_filtered, features= "percent.mt") + theme_classic() + theme(legend.position = 'none', axis.title.x=element_blank()) + scale_fill_manual(values = alpha("#5ec962", 0.5))
p6$layers[[2]]$aes_params$alpha = 0.2
p7 <- ggplot(data = SOB_filtered@meta.data) +
    geom_point(mapping = aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), position = 'jitter', alpha = 0.3) +
    scale_color_viridis_b(option="magma", limits = c(0, 13)) +
    theme_classic() + guides(shape = guide_legend(override.aes = list(size = 2)), 
                             color = guide_legend(override.aes = list(size = 2))) +
    theme(legend.title = element_text(size = 10), legend.text = element_text(size = 9))
  
title_grob <- textGrob(paste0("QC for ",i), gp = gpar(fontsize = 16, fontface = "bold"))
subtitle_grob <- textGrob(paste0(nrow(SOB_filtered), " features x ", ncol(SOB_filtered), " cells"), gp = gpar(fontsize = 12))
final_plot <- gridExtra::arrangeGrob(p1, t, p2, p3, p4, p6, p5, p7, ncol = 2, nrow = 4)
margin <- unit(0.5, "line")
updated_plot_combined_gof <-  grid.arrange(title_grob, 
                                           subtitle_grob, 
                                           final_plot,
                                           heights = unit.c(grobHeight(title_grob) + 1.2*margin, 
                                                            grobHeight(subtitle_grob) + margin,
                                                            unit(1,"null")))

png(filename = paste0(arg$workdir, "/outdir/QualityControl_graph_", i,".png"), width = 9.5, height = 5.5, units = "in", res = 300)
  plot(updated_plot_combined_gof)
dev.off()

write10xCounts(SOB_filtered, (arg$workdir, "/filtered_matrix_", i, "/")

