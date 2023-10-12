

p <- add_argument(p, "--marker_genes", help = "Path to cell type marker gene file", type = "character", nargs = '+')

arg <- parse_args(p)

SOB_filtered <- readRDS(SOB_filtered, (arg$workdir, "/outdir/", i, "_processed.rds")

# 6. Cell type annotation
markers <- fread(cellmarker_input)
unique_values <- unique(markers$Type)
num_colors <- length(unique_values)
color_palette <- c("mediumblue", "#762a83", "turquoise1", "#e31a1c", "#ff7f00","#33a02c", "#b2df8a", "yellow", "#c51b7d")
color_palette <- color_palette[1:num_colors]
color_list <- setNames(color_palette, unique_values)

varFeats <- FindVariableFeatures(SOB_filtered, selection.method = "vst", nfeatures = 5000)
topgenes <- head(VariableFeatures(varFeats), 5000)

expr <- GetAssayData(SOB_filtered, slot = "data", assay = "RNA")
topgene_filter <- rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1 # Filter genes with zero expression in more than 10% ot total cells
table(topgene_filter)
topgenes <- topgenes[topgene_filter]
  
picked_genes <- unique(c(markers$Marker, topgenes))
expr <- expr[rownames(expr) %in% picked_genes, ]
rts <- scSorter(expr, markers)
rm(expr, varFeats, topgenes, topgene_filter, picked_genes)
df <- data.frame(rbind(table(as.character(rts$Pred_Type))), check.names = FALSE) %>% gather(key = "CellType", value = "CellNumber")
df$Percentage <- round((df$CellNumber / sum(df$CellNumber)) * 100, 1)

p2 <- ggplot(df, aes(x="", y=CellNumber, fill=CellType)) +
    geom_bar(width = 1, stat = "identity") + 
    coord_polar("y") +
    geom_text(aes(x = 1.4, label = paste(CellType, "\n", Percentage, "%")), position = position_stack(vjust = 0.5), size = 1.8) +  # Adjust vjust to move labels closer to the outside
    scale_fill_manual(values = color_list) +
    theme_void() + theme(legend.position = 'none', axis.text.x=element_blank()) +
    labs(title = paste0("        % Cell Types in ", i))
    
SOB_filtered@meta.data$celltypes <- rts$Pred_Type
SOB_filtered@meta.data$cluster <- Idents(SOB_filtered)
Idents(SOB_filtered) <- rts$Pred_Type
  
p3 <- DimPlot(SOB_filtered, reduction = "umap", cols = color_list, pt.size = 0.5) + labs(title = paste0("UMAP plot ", i))  + theme_classic() + 
    guides(shape = guide_legend(override.aes = list(size = 2)), color = guide_legend(override.aes = list(size = 2))) +
    theme(legend.title = element_text(size = 10), legend.text = element_text(size = 9))
  
final_plot <- p2 + p3
final_plot <- final_plot + plot_layout(guides = "collect")

png(filename = paste0(arg$workdir, "/outdir/CellType_graph_", i,".png"), width = 9.5, height = 5.5, units = "in", res = 300)
  plot(final_plot)
dev.off()

saveRDS(SOB_filtered, (arg$workdir, "/outdir/", i, "processed_celltyped.rds")

