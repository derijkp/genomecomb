proc sc_celltyper_scsorter_job {args} {
	upvar job_logdir job_logdir
	set cellmarkerfile ""
	set tissue ""
	set extrarootname scsorter
	cg_options sc_celltyper_scsorter args {
		-cellmarkerfile {set cellmarkerfile [file_absolute $value]}
		-tissue {set tissue $value}
		-extrarootname {set extrarootname $value}
	} {scgenefile scisoformfile} 2 2
	#
	if {$cellmarkerfile eq ""} {
		error "celmarkerfile has to be given for sc_celltyper_scsorter"
	}
	set scgenefile [file_absolute $scgenefile]
	set scisoformfile [file_absolute $scisoformfile]
	set rootname [file_rootname $scgenefile]
	set dir [file dir $scgenefile]
	set scgenefile10x [file root [gzroot $scgenefile]].10x
	set groupfile $dir/sc_group-$extrarootname-$rootname.tsv
	set umappng $dir/sc_celltype_umap-$extrarootname-$rootname.png
	set outrds $dir/sc_celltype-$extrarootname-$rootname.rds
	set tempresult $groupfile.temp
	job sc_celltyper_scsorter_scgenefile10x-$extrarootname-$rootname -deps {
		$scgenefile
	} -targets {
		$scgenefile10x
	} -vars {
		scgenefile scgenefile10x
	} -code {
		cg tsv210x $scgenefile $scgenefile10x.temp
		file rename -force $scgenefile10x.temp $scgenefile10x
	}
	job sc_celltyper_$extrarootname-$rootname -mem 10G -deps {
		$scgenefile10x
	} -targets {
		$groupfile $umappng
	} -vars {
		scgenefile10x umappng cellmarkerfile tissue rootname outrds tempresult groupfile
	} -code {
		set out10x {}
		set metadatafile {}
		set outRfile ~/tmp/sc_filter.R
		set outRfile {}
		R -outRfile $outRfile -vars {
			scgenefile10x umappng cellmarkerfile tissue rootname outrds tempresult
		} {
			suppressMessages(library(Seurat))
			suppressMessages(library(tidyverse))
			suppressMessages(library(scSorter))
			suppressMessages(library(patchwork))
			ggplotColours <- function(n = 6, h = c(0, 360) + 15){
			  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
			  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
			}
			#
			# read markers (and flag errors if present)
			cellmarkers <- read.table(cellmarkerfile,sep='\t',header=TRUE)
			if (tissue != "") {
				if (length(cellmarkers$tissue) == 0) {
					stop("cellmarker file has no column named tissue")
				}
				cellmarkers = cellmarkers[cellmarkers$tissue == tissue]
			}
			markergenes=cellmarkers$marker
			if (length(markergenes) == 0) {
				markergenes=cellmarkers$gene
			}
			if (length(markergenes) == 0) {
				stop("cellmarker file has no column named marker (or gene)")
			}
			if (length(cellmarkers$celltype) == 0) {
				stop("cellmarker file has no column named celltype")
			}
			markers <- data.frame(Type=cellmarkers$celltype,Marker=markergenes)
			#
			# 1. Load sample
			mat=Read10X(scgenefile10x)
			mito.genes <- grep(pattern = "^MT-", x = rownames(mat), value = TRUE)
			RPS.genes <- grep(pattern = "^RPS", x = rownames(mat), value = TRUE)
			qc <- c(mito.genes, RPS.genes)
			keep <- !row.names(mat) %in% unlist(qc)
			filt_mat <- mat[keep, ]
			rm(mat)
			#
			# Seurat object		
			SOB <- CreateSeuratObject(filt_mat)
			rm(filt_mat)
			SOB[["percent.mt"]] <- PercentageFeatureSet(SOB, pattern = "^MT-")
			#
			# Mitochondrial cut-off and standard processing (scDblFinder does not prefer this processing step first, DoubletFinder does)
			cat("running Mitochondrial cutoff")
			SOB_filtered <- subset(SOB, subset = percent.mt < 5)
			SOB@assays$RNA_org <- CreateAssayObject(SOB@assays$RNA@counts)
			SOB@assays$RNA_org@key <- "newkey_"
			genes_to_remove <- rownames(SOB)[grep("^novelg-", rownames(SOB))]
			SOB <- SOB[!rownames(SOB) %in% genes_to_remove, ]
			SOB <- NormalizeData(SOB, normalization.method = "LogNormalize", scale.factor = 10000)
			SOB <- FindVariableFeatures(SOB, selection.method = "vst", nfeatures = 2000)
			SOB <- ScaleData(SOB)
			SOB <- RunPCA(SOB)
			SOB <- FindNeighbors(SOB, reduction = "pca", dims = 1:50)
			SOB <- FindClusters(SOB, resolution = 0.5)
			SOB <- RunUMAP(SOB, dims = 1:50)
			SOB <- RunTSNE(SOB, dims = 1:50)
			SOB_filtered <- SOB
			
			# 6. Cell type annotation
			if ("weight" %in% colnames(cellmarkers)) {
				markers$Weight = cellmarkers$weight
				markers$Weight[which(is.na(markers$Weight))] = 1
			}
			if ("markertype" %in% colnames(cellmarkers)) {
				markers = markers[cellmarkers$markertype == "up" | cellmarkers$markertype == "",]
			}
			#
			varFeats <- FindVariableFeatures(SOB_filtered, selection.method = "vst", nfeatures = 5000)
			topgenes <- head(VariableFeatures(varFeats), 5000)
			#
			expr <- GetAssayData(SOB_filtered, slot = "data", assay = "RNA")
			topgene_filter <- rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1 # Filter genes with zero expression in more than 10% ot total cells
			table(topgene_filter)
			topgenes <- topgenes[topgene_filter]
			#
			picked_genes <- unique(c(markers$Marker, topgenes))
			expr <- expr[rownames(expr) %in% picked_genes, ]
			markerspresent = markers$Marker %in% rownames(expr)
			if (!all(markerspresent)) {
				missing = markers$Marker[!markerspresent]
				warning(paste("skipping markers missing in data set:",paste(missing,collapse=",")))
				cmarkers = markers[markers$Marker %in% rownames(expr), ]
				cmarkers=droplevels(cmarkers)
				rownames(cmarkers) = 1:nrow(cmarkers)
			} else {
				cmarkers = markers
			}
			rts <- scSorter(expr, cmarkers)
			unique_values <- unique(cmarkers$Type)
			num_colors <- length(unique_values)
			color_palette <- c("mediumblue", "#762a83", "turquoise1", "#e31a1c", "#ff7f00","#33a02c", "#b2df8a", "yellow", "#c51b7d")
			if (num_colors <= length(color_palette)) {
				color_palette <- color_palette[1:num_colors]
			} else {
				color_palette <- ggplotColours(n=num_colors)
			}
			color_list <- setNames(color_palette, unique_values)
			rm(expr, varFeats, topgenes, topgene_filter, picked_genes)
			df <- data.frame(rbind(table(as.character(rts$Pred_Type))), check.names = FALSE) %>% gather(key = "CellType", value = "CellNumber")
			df$Percentage <- round((df$CellNumber / sum(df$CellNumber)) * 100, 1)
			#
			p2 <- ggplot(df, aes(x="", y=CellNumber, fill=CellType)) +
			    geom_bar(width = 1, stat = "identity") + 
			    coord_polar("y") +
			    geom_text(aes(x = 1.4, label = paste(CellType, "\n", Percentage, "%")), position = position_stack(vjust = 0.5), size = 1.8) +  # Adjust vjust to move labels closer to the outside
			    scale_fill_manual(values = color_list) +
			    theme_void() + theme(legend.position = 'none', axis.text.x=element_blank()) +
			    labs(title = paste0("        % Cell Types in ",rootname))
			    
			SOB_filtered@meta.data$celltypes <- rts$Pred_Type
			SOB_filtered@meta.data$cluster <- Idents(SOB_filtered)
			Idents(SOB_filtered) <- rts$Pred_Type
			#
			p3 <- DimPlot(SOB_filtered, reduction = "umap", cols = color_list, pt.size = 0.5) + 
				labs(title = paste0("UMAP plot ", rootname))  + theme_classic() + 
				guides(shape = guide_legend(override.aes = list(size = 2)), color = guide_legend(override.aes = list(size = 2))) +
				theme(legend.title = element_text(size = 10), legend.text = element_text(size = 9))
			final_plot <- p2 + p3
			final_plot <- final_plot + plot_layout(guides = "collect")
			#
			final_plot
			ggsave(umappng, width = 9.5, height = 5.5, units = "in", dpi = 300)
			#
			saveRDS(SOB_filtered, outrds)
			# output tsv
			qc_data <- FetchData(SOB, vars = c("ident", "UMAP_1", "UMAP_2"))
			out=data.frame(
				cell=rownames(SOB@meta.data),
				group=rts$Pred_Type,
				group_filtered=rts$Pred_Type,
				score=rep(NA,nrow(SOB@meta.data)),
				ncells=rep(NA,nrow(SOB@meta.data)),
				UMAP_1=qc_data$UMAP_1,
				UMAP_2=qc_data$UMAP_2,
				seurat_clusters=SOB@meta.data$seurat_clusters,
				nCount_RNA=SOB@meta.data$nCount_RNA,
				nFeature_RNA=SOB@meta.data$nFeature_RNA,
				percent.mt=SOB@meta.data$percent.mt,
				RNA_snn_res.0.5=SOB@meta.data$RNA_snn_res.0.5
			)
			for(type in unique(rts$Pred_Type)){
				out$ncells[out$group == type] = sum(rts$Pred_Type == type)
			}
			write.table(out,file=tempresult,sep="\t",quote=FALSE,row.names = FALSE)
		}
		file rename -force $tempresult $groupfile
	}
	sc_pseudobulk_job $scgenefile $scisoformfile $groupfile
}

proc cg_sc_celltyper_scsorter {args} {
	set args [job_init {*}$args]
	sc_celltyper_scsorter_job {*}$args
	job_wait
}
