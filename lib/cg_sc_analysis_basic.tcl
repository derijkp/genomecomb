proc sc_analysis_basic_job {args} {
	upvar job_logdir job_logdir
	#
	set organelles {}
	cg_options sc_analysis_basic args {
		-reftranscripts {
			set reftranscripts $value
		}
		-organelles {
			set organelles $value
		}
	} {scgenefile scisoformfile} 2 2
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
	job sc_analysis_basic_scgenefile10x-$rootname -deps {
		$scgenefile
	} -targets {
		$scgenefile10x
	} -vars {
		scgenefile scgenefile10x
	} -code {
		cg tsv210x $scgenefile $scgenefile10x.temp
		file rename $scgenefile10x.temp $scgenefile10x
	}
	job sc_analysis_basic_scisoformfile10x-$rootname -deps {
		$scisoformfile
	} -targets {
		$scisoformfile10x
	} -vars {
		scisoformfile scisoformfile10x
	} -code {
		set featurefields [findfields [cg select -h $scisoformfile] {gene transcript}]
		cg tsv210x -featurefields $featurefields -countfield counts_weighed $scisoformfile $scisoformfile10x.temp
		file rename $scisoformfile10x.temp $scisoformfile10x
	}
	set outinfofile $dir/sc_temp_info-$rootname.tsv
	set umappng $dir/sc_umap-$rootname.png
	set tsnepng $dir/sc_tsne-$rootname.png
	ref_transcripts_convert $reftranscripts tsvreftranscripts gtfreftranscripts
	set theader [cg select -h $tsvreftranscripts]
	set cfield [lindex $theader [tsv_basicfields $theader 1]]
	if {[llength $organelles] <= 1} {
		set mitchr [lindex $organelles 0]
	} else {
		set pos [list_find -regexp $organelles M]
		set mitchr [lindex $organelles [lindex $pos 0]]
	}
	set genefields [findfields [cg select -h $tsvreftranscripts] {gene}]
	set mitgenes [cg select -sh /dev/null -f $genefields -q "\$$cfield eq \"$mitchr\"" $tsvreftranscripts]
	job sc_analysis_basic_r-$rootname -deps {
		$scgenefile10x
		$scisoformfile10x
	} -targets {
		$outinfofile $outrds $umappng
	} -vars {
		scgenefile10x scisoformfile10x expectedcells metadatafile outrds out10x outkneeplot outinfofile
		umappng tsnepng mitgenes
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
		R -outRfile $outRfile -listvars {
			mitgenes
		} -vars {
			scgenefile10x scisoformfile10x expectedcells metadatafile outrds out10x outkneeplot outinfofile
			umappng tsnepng
		} {
			suppressMessages(library(tidyverse))
			suppressMessages(library(argparser))
			suppressMessages(library(DropletUtils))
			suppressMessages(library(Seurat))
			# not yet
			usescrublet <- 0
			#
			# 1. Load sample
			mat=Read10X(scgenefile10x)
			#
			# Seurat object
			# to avoid "Feature names of counts matrix cannot be empty" error, remove these rows
			emptynames=which(rownames(mat) == '')
			if (length(emptynames) > 0) {
				mat = mat[-emptynames,]
			}
			# create seurat object
			SOB <- CreateSeuratObject(mat)
			mitgenes=mitgenes[mitgenes %in% rownames(GetAssayData(SOB))]
			if (length(mitgenes) > 0) {
				SOB[["percent.mt"]] <- PercentageFeatureSet(SOB, features = mitgenes)
			} else {
				SOB[["percent.mt"]] <- 0
			}
			#
			# Mitochondrial cut-off and standard processing
			cat("running Mitochondrial cutoff")
			# error if no mitochondrial genes were found
			# SOB_filtered <- subset(SOB, subset = percent.mt < 5)
			SOB@assays$RNA_org <- CreateAssayObject(SOB@assays$RNA@counts)
			SOB@assays$RNA_org@key <- "newkey_"
			genes_to_remove <- rownames(SOB)[grep("^novelg-", rownames(SOB))]
			SOB <- SOB[!rownames(SOB) %in% genes_to_remove, ]
			SOB <- NormalizeData(SOB, normalization.method = "LogNormalize", scale.factor = 10000)
			SOB <- FindVariableFeatures(SOB, selection.method = "vst", nfeatures = 2000)
			SOB <- ScaleData(SOB)
			if (ncol(SOB) > 50) {npcs = 50} else {npcs=ncol(SOB)-1}
			SOB <- RunPCA(SOB, npcs = npcs)
			SOB <- FindNeighbors(SOB, reduction = "pca", dims = 1:npcs)
			SOB <- FindClusters(SOB, resolution = 0.5)
			SOB <- RunUMAP(SOB, dims = 1:npcs)
			ggsave(umappng,DimPlot(SOB))
			e.out <- tryCatch({
				SOB <- RunTSNE(SOB, dims = 1:npcs)
				ggsave(tsnepng,DimPlot(SOB,reduction="tsne"))
			}, error = function(err) {
				cat("warning: RunTSNE failed\n")
			})
			cat(paste0("Quality control for ", scgenefile10x, "Completed!"), "\n")
			#
		}
	}
}

proc cg_sc_analysis_basic {scgenefile scisoformfile expectedcells} {
	set args [job_init {*}$args]
	sc_analysis_basic_job {*}$args
	job_wait
}
