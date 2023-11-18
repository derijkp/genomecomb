proc sc_filter_default_job {args} {
	upvar job_logdir job_logdir
	#
	cg_options sc_filter_default args {
	} {scgenefile scisoformfile expectedcells} 3 3
	if {$expectedcells eq ""} {
		error "failed sc_filter_default: please give expected number of cells (-sc_expectedcells) for single cell analysis"
	}
	#
	set scgenefile [file_absolute $scgenefile]
	set scisoformfile [file_absolute $scisoformfile]
	set rootname [file_rootname $scgenefile]
	set dir [file dir $scgenefile]
	set out10x $dir/sc_gene_counts_filtered-$rootname.10x
	set target $dir/sc_gene_counts_filtered-$rootname.tsv.zst
	set outrds $dir/sc_seurat_filtered-$rootname.rds
	set outkneeplot $dir/kneeplot-$rootname.png
	set scgenefile10x [file root [gzroot $scgenefile]].10x
	set scisoformfile10x [file root [gzroot $scisoformfile]].weighed_count.10x
	set outfile [file root [gzroot $scgenefile]].rds
	job sc_filter_default_scgenefile10x-$rootname -deps {
		$scgenefile
	} -targets {
		$scgenefile10x
	} -vars {
		scgenefile scgenefile10x
	} -code {
		cg tsv210x $scgenefile $scgenefile10x.temp
		file rename $scgenefile10x.temp $scgenefile10x
	}
	job sc_filter_default_scisoformfile10x-$rootname -deps {
		$scisoformfile
	} -targets {
		$scisoformfile10x
	} -vars {
		scisoformfile scisoformfile10x
	} -code {
		cg tsv210x -countfield counts_weighed $scisoformfile $scisoformfile10x.temp
		file rename $scisoformfile10x.temp $scisoformfile10x
	}
	set outinfofile $dir/sc_temp_info-$rootname.tsv
	set umappng $dir/sc_umap-$rootname.png
	set tsnepng $dir/sc_tsne-$rootname.png
	job sc_filter_default_r-$rootname -deps {
		$scgenefile10x
		$scisoformfile10x
	} -targets {
		$outinfofile $outrds $umappng $tsnepng
	} -vars {
		scgenefile10x scisoformfile10x expectedcells metadatafile outrds out10x outkneeplot outinfofile
		umappng tsnepng
	} -code {
		# clear output
		foreach out [list $out10x $outinfofile $outrds] {
			if {[file exists $out]} {
				catch {file delete -force $out.old}
				file rename $out $out.old
			}
		}
	
		set out10x {}
		set metadatafile {}
		set outRfile ~/tmp/sc_filter.R
		set outRfile {}
		R -outRfile $outRfile -vars {
			scgenefile10x scisoformfile10x expectedcells metadatafile outrds out10x outkneeplot outinfofile
			umappng tsnepng
		} {
			suppressMessages(library(tidyverse))
			suppressMessages(library(argparser))
			suppressMessages(library(DropletUtils))
			suppressMessages(library(DoubletFinder))
			suppressMessages(library(scDblFinder))
			suppressMessages(library(Seurat))
			# not yet
			usescrublet <- 0
			#
			# 1. Load sample
			mat=Read10X(scgenefile10x)
			#
			# kneeplot
			bcrank <- barcodeRanks(mat)
			inflection <- metadata(bcrank)$inflection
			knee <- metadata(bcrank)$knee
			if (outkneeplot != "") {
				png(file=outkneeplot)
				uniq <- !duplicated(bcrank$rank)
				plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy", xlab="Rank", ylab="Total UMI count", cex.lab=1.2)
				abline(h=inflection, col="darkgreen", lty=2)
				abline(h=knee, col="dodgerblue", lty=2)
				legend("bottomleft", legend=c(paste0("Inflection (",round(inflection),")"), paste0("Knee (",knee,")")), col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)
				dev.off()
			}
			#
			# 2. EmptyDrops
			cat("running EmptyDrops\n")
			mito.genes <- grep(pattern = "^MT-", x = rownames(mat), value = TRUE)
			RPS.genes <- grep(pattern = "^RPS", x = rownames(mat), value = TRUE)
			qc <- c(mito.genes, RPS.genes)
			keep <- !row.names(mat) %in% unlist(qc)
			#
			# not very clean, but on some data emptyDropsCellRanger gives an error (that I could not easily patch)
			# so use emptyDrops as a fallback where this is the case
			emptydrops.error=0
			rm(e.out) ; rm(is.cell)
			e.out <- tryCatch({
				emptyDropsCellRanger(mat[keep, ],n.expected.cells=expectedcells)
			}, error = function(err) {
				cat("warning: emptyDropsCellRanger failed, fallback to emptyDrops\n")
			})
			if (!is.null(e.out)) {
				is.cell <- e.out$FDR <= 0.001
				is.cell[is.na(is.cell)] <- FALSE
				cat("detected",sum(is.cell),"cells using emptyDropsCellRanger")
				is.cell.method="emptyDropsCellRanger"
				is.cell.lower=NA
			} else {
				e.out <- tryCatch({
					sortedcols = sort(colSums(mat),decreasing=TRUE)
					lower=sortedcols[expectedcells]/10
					if (lower < 100) {lower = 100}
					e.out <- emptyDrops(mat[keep,],niters=20000,lower=lower,test.ambient=TRUE)
				}, error = function(err) {
					cat("warning: emptyDrops failed, fallback to inflectionpoint kneeplot\n")
				})
				if (!is.null(e.out)) {
					is.cell <- e.out$FDR <= 0.001
					is.cell[is.na(is.cell)] <- FALSE
					cat("detected",sum(is.cell),"cells using emptyDrops")
					is.cell.method="emptyDrops"
					is.cell.lower=lower
				} else {
					emptydrops.error=1
					is.cell <- colSums(mat) >= knee
					cat("detected",sum(is.cell),"cells using knee")
					is.cell.method="knee"
					is.cell.lower=knee
				}
			}
			filt_mat <- mat[, is.cell]
			rm(mat, barcodes, features, e.out, is.cell)
			#
			# Seurat object		
			SOB <- CreateSeuratObject(filt_mat)
			rm(filt_mat)
			SOB[["percent.mt"]] <- PercentageFeatureSet(SOB, pattern = "^MT-")
			#
			# 3. Doublet identification
			# 3.1 scDblFinder
			cat("running scDblFinder\n")
			sce_tmp <- as.SingleCellExperiment(SOB)
			sce_tmp <- scDblFinder(sce_tmp, dbr.sd = 1)
			seurat_tmp <- as.Seurat(sce_tmp, counts = "counts", data = NULL)
			#
			SOB@meta.data <- seurat_tmp@meta.data
			rm(seurat_tmp, sce_tmp)
			SOB@meta.data$scDblFinder.class <- gsub("singlet", "Singlet", SOB@meta.data$scDblFinder.class)
			SOB@meta.data$scDblFinder.class <- gsub("doublet", "Doublet", SOB@meta.data$scDblFinder.class)
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
			ggsave(umappng,DimPlot(SOB))
			ggsave(tsnepng,DimPlot(SOB,reduction="tsne"))
			# plot=DimPlot(SOB)
			# htmlwidgets::saveWidget(ggplotly(plot),umapplot)
			# HoverLocator(plot = plot, information = FetchData(SOB, vars = c("ident", "nCount_RNA", "nFeature_RNA", "DoubletFinder")))
			#
			# 3.2 DoubletFinder
			cat("running DoubletFinder")
			sweep.res.list_SOB <- paramSweep_v3(SOB, PCs = 1:50, sct = FALSE)
			sweep.stats_SOB <- summarizeSweep(sweep.res.list_SOB, GT = FALSE)
			bcmvn_SOB <- find.pK(sweep.stats_SOB)
			optimal_pK <- as.numeric(as.character(bcmvn_SOB$pK[which.max(bcmvn_SOB$BCmetric)]))
			homotypic.prop <- modelHomotypic(SOB$seurat_clusters)
			expected_doublet_rate_per_10000_cells <- 0.075
			total_cells <- ncol(SOB)
			nExpDoublets <- (expected_doublet_rate_per_10000_cells * (total_cells))
			nExpDoublets.adj <- round(nExpDoublets*(1-homotypic.prop))
			#
			SOB <- doubletFinder_v3(SOB, pN = 0.25, pK = optimal_pK, nExp = nExpDoublets, PCs = 1:50)
			colnames(SOB@meta.data)[colnames(SOB@meta.data) == grep("^DF.classifications_", colnames(SOB@meta.data), value = TRUE)] <- "DoubletFinder"
			SOB <- doubletFinder_v3(SOB, pN = 0.25, pK = optimal_pK, nExp = nExpDoublets.adj, PCs = 1:50)
			colnames(SOB@meta.data)[colnames(SOB@meta.data) == grep("^DF.classifications_", colnames(SOB@meta.data), value = TRUE)] <- "DoubletFinder.adj"
			#
			if (usescrublet == 1)	{
				# 3.3 Scrublet
				cat("running Scrublet")
				reticulate::use_condaenv("~/miniconda3/envs/scrublet")
				scrub <- reticulate::import(module = "scrublet", convert = FALSE,delay_load = TRUE)
				counts_matrix = reticulate::r_to_py(Seurat::GetAssayData(object = SOB[["RNA"]], slot = "counts"))$T$tocsc()
				
				scr <- eval(rlang::expr(scrub$Scrublet(counts_matrix = counts_matrix, expected_doublet_rate = expected_doublet_rate_per_10000_cells))) 
				doublet_results <- eval(rlang::expr(reticulate::py_to_r(scr$scrub_doublets())))
				doublet_score <- doublet_results[[1]]
				names(doublet_score) <- colnames(SOB)
				doublet_prediction <- doublet_results[[2]]
				names(doublet_prediction) <- colnames(SOB)
				SOB <- Seurat::AddMetaData(SOB, metadata = doublet_score, col.name = "scrublet_doublet_score")
				SOB <- Seurat::AddMetaData(SOB, metadata = doublet_prediction, col.name = "scrublet_doublet_prediction")
				SOB@meta.data$scrublet_doublet_prediction <- gsub("FALSE", "Singlet", SOB@meta.data$scrublet_doublet_prediction)
				SOB@meta.data$scrublet_doublet_prediction <- gsub("TRUE", "Doublet", SOB@meta.data$scrublet_doublet_prediction)
				rm(counts_matrix, scr, doublet_results, doublet_score, doublet_prediction)
				
				# Doublet Score
				SOB@meta.data <- SOB@meta.data %>% mutate(doublet_score = ifelse(scDblFinder.class == "Doublet", 1,0) + ifelse(scrublet_doublet_prediction == "Doublet", 1, 0) + ifelse(DoubletFinder == "Doublet", 0.3, 0) + ifelse(DoubletFinder.adj == "Doublet", 0.7, 0))
			} else {
				SOB@meta.data <- SOB@meta.data %>% mutate(doublet_score = ifelse(scDblFinder.class == "Doublet", 1,0) + ifelse(DoubletFinder == "Doublet", 0.3, 0) + ifelse(DoubletFinder.adj == "Doublet", 0.7, 0))
			}
			cat(paste0("Quality control for ", scgenefile10x, "Completed!"), "\n")
			#
			if (metadatafile != "") {
				# 4. Add Metadata
				cat("Adding Metadata")
				sample_info <- fread(metadatafile)
				SOB@meta.data <- merge(SOB@meta.data, sample_info, by.x = "dnumber", by.y = "dnumber", all = TRUE)
				rm(sample_info)
			}
			#
			# 5. Add transcript data
			cat("Add transcript data")
			transcripts <- Read10X(scisoformfile10x, gene.column = 1)
			barcodes_to_keep <- colnames(transcripts) %in% colnames(SOB)
			transcripts <- transcripts[, barcodes_to_keep]
			SOB[["ISOFORM"]] <- CreateAssayObject(counts = transcripts)
			rm(transcripts)
			cat("Metadata and isoforms added", "\n")
			#
			if (outrds != "") {
				# 6. Save preprocessed object
				cat(paste0("Save rds in ",outrds,"\n"))
				# name <- paste("SOB_processed_", unique(SOB$Sample), sep = "")
				saveRDS(object = SOB, file = outrds)
			}
			#
			if (out10x != "") {
				cat("Writing 10x output to", out10x, "\n")
				write10xCounts(out10x,SOB_filtered@assays$RNA@counts,version="3")
			}
			#
			cat("Writing info file", outinfofile, "\n")
			info=transform(SOB@meta.data, cell=rownames(SOB@meta.data))
			write.table(info,file=outinfofile,sep="\t",quote=FALSE,row.names = FALSE)
			cat("is_cell_method\tis_cell_lower\n",is.cell.method,"\t",is.cell.lower,"\n",sep="",file=paste0(outinfofile,".analysisinfo"))
			cat(paste0("Object saved, next!"), "\n")
		}
	}
	#
	# make cell info file (sc_cellinfo_raw-$rootname.tsv.zst and sc_cellinfo_filtered-$rootname.tsv.zst)
	set reads_per_cell [jobgzfile $dir/reads_per_cell_raw-$rootname.tsv $dir/reads_per_cell-$rootname.tsv $dir/reads_per_cell.tsv]
	set umis_per_cell [jobgzfile $dir/umis_per_cell_raw-$rootname.tsv $dir/umis_per_cell-$rootname.tsv $dir/umis_per_cell.tsv]
	# -------------------
	job sc_filter_write_cellinfo-$rootname -deps {
		$outinfofile
		$reads_per_cell $umis_per_cell
	} -targets {
		$dir/sc_cellinfo_raw-$rootname.tsv.zst
		$dir/sc_cellinfo_filtered-$rootname.tsv.zst
	} -vars {
		dir outinfofile rootname reads_per_cell umis_per_cell
	} -code {
		catch {gzclose $o} ; catch {gzclose $f}
		set f [open $outinfofile]
		set header [tsv_open $f]
		set newfields [list cell {*}[list_remove $header orig.ident ident cell]]
		set poss [list_cor $header $newfields]
		set newfields [lrange $newfields 1 end]
		set empty [list_fill [llength $poss] {}]
		unset -nocomplain a
		while 1 {
			if {[gets $f line] == -1} break
			set line [list_sub [split $line \t] $poss]
			set cell [lindex $line 0]
			set a($cell) [lrange $line 1 end]
		}
		close $f
		unset -nocomplain umia
		set temp [file_read $umis_per_cell]
		regsub \n\$ $temp "" temp
		array set umia [split $temp \n\t]
		set f [open $reads_per_cell]
		set header [tsv_open $f]
		set o [wgzopen $dir/sc_cellinfo_raw-$rootname.tsv.zst]
		puts $o cell\treadcount\tumicount\tis_cell\t[join $newfields \t]
		while 1 {
			if {[gets $f line] == -1} break
			set line [split $line \t]
			foreach {cell readcount} $line break
			lappend line [get umia($cell) 0]
			if {[info exists a($cell)]} {
				lappend line 1 {*}$a($cell)
			} else {
				lappend line 0 {*}$empty
			}
			puts $o [join $line \t]
		}
		gzclose $o
		gzclose $f
		cg select -q {$is_cell == 1} $dir/sc_cellinfo_raw-$rootname.tsv.zst $dir/sc_cellinfo_filtered-$rootname.tsv.zst
	}
	# 
	# make filtered gene files (sc_gene_counts_filtered-$rootname.tsv.zst)
	# ------------------------
	# get cells to keep (in array a)
	job sc_counts_filtered-$rootname -deps {
		$scgenefile $scisoformfile
		$dir/sc_cellinfo_filtered-$rootname.tsv.zst
	} -targets {
		$dir/sc_gene_counts_filtered-$rootname.tsv.zst
		$dir/sc_isoform_counts_filtered-$rootname.tsv.zst
	} -vars {
		dir outinfofile rootname scgenefile scisoformfile
	} -code {
		set cells [cg select -sh /dev/null -f cell -q {$is_cell == 1} $dir/sc_cellinfo_filtered-$rootname.tsv.zst]
		unset -nocomplain a
		foreach cell $cells {
			set a($cell) 1
		}
		# write sc_gene_counts_filtered
		catch {gzclose $o} ; catch {gzclose $f}
		set o [wgzopen $dir/sc_gene_counts_filtered-$rootname.tsv.temp.zst]
		set f [gzopen $scgenefile]
		set header [tsv_open $f]
		set pos [lsearch $header cell]
		puts $o [join $header \t]
		while 1 {
			if {[gets $f line] == -1} break
			set cell [lindex [split $line \t] $pos]
			if {![info exists a($cell)]} continue
			puts $o $line
		}
		gzclose $f
		gzclose $o
		file rename -force $dir/sc_gene_counts_filtered-$rootname.tsv.temp.zst $dir/sc_gene_counts_filtered-$rootname.tsv.zst
		# make filtered isoform file (sc_isoform_counts_filtered-$rootname.tsv.zst)
		# ------------------------------------
		# write sc_isoform_counts_filtered
		catch {gzclose $o} ; catch {gzclose $f}
		set o [wgzopen $dir/sc_isoform_counts_filtered-$rootname.tsv.temp.zst]
		set f [gzopen $scisoformfile]
		set header [tsv_open $f]
		set pos [lsearch $header cell]
		puts $o [join $header \t]
		while 1 {
			if {[gets $f line] == -1} break
			set cell [lindex [split $line \t] $pos]
			if {![info exists a($cell)]} continue
			puts $o $line
		}
		gzclose $f
		gzclose $o
		file rename -force $dir/sc_isoform_counts_filtered-$rootname.tsv.temp.zst $dir/sc_isoform_counts_filtered-$rootname.tsv.zst
	}
	return $dir/sc_cellinfo_filtered-$rootname.tsv.zst
}

proc cg_sc_filter_default {scgenefile scisoformfile expectedcells} {
	set args [job_init {*}$args]
	sc_filter_default_job {*}$args
	job_wait
}
