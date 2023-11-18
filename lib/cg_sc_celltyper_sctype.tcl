proc sc_celltyper_sctype_job {args} {
	upvar job_logdir job_logdir
	set tissue ""
	# tissue can be e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
	cg_options sc_celltyper_sctype args {
		-tissue {set tissue $value}
	} {scgenefile scisoformfile} 2 2
	#
	set scgenefile [file_absolute $scgenefile]
	set scisoformfile [file_absolute $scisoformfile]
	set rootname [file_rootname $scgenefile]
	set dir [file dir $scgenefile]
	set scgenefile10x [file root [gzroot $scgenefile]].10x
	set groupfile $dir/sc_group-sctype-$rootname.tsv
	set tempresult $groupfile.temp
	job sctype_scgenefile10x-$rootname -deps {
		$scgenefile
	} -targets {
		$scgenefile10x
	} -vars {
		scgenefile scgenefile10x
	} -code {
		cg tsv210x $scgenefile $scgenefile10x.temp
		file rename $scgenefile10x.temp $scgenefile10x
	}
	set umappng $dir/sc_celltype_umap-$rootname.png
	set R [file_resolve [findR]]
	set dirR [file dir $R]
	set sctypedir [lindex [bsort [glob $dirR/lib64/R/library/sc-type-*]] end]
	job sc_celltyper_sctype-$rootname -deps {
		$scgenefile10x
	} -targets {
		$groupfile $umappng
	} -vars {
		sctypedir scgenefile10x tissue umappng tempresult groupfile
	} -code {
		set out10x {}
		set metadatafile {}
		set outRfile ~/tmp/sc_celltyper.R
		set outRfile {}
		R -outRfile $outRfile -vars {
			scgenefile10x tissue umappng tempresult sctypedir
		} {
			library(Seurat)
			library(tidyverse)
			library(sleepwalk)
			theme_set(theme_bw(base_size = 14))
			
			data = Read10X(scgenefile10x)
			SOB <- CreateSeuratObject(data)
			SOB[["percent.mt"]] <- PercentageFeatureSet(SOB, pattern = "^MT-")
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
			# DimPlot(SOB, reduction = "umap")
			
			# load gene set preparation function
			lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
			source(paste0(sctypedir,"/R/gene_sets_prepare.R"))
			# load cell type annotation function
			source(paste0(sctypedir,"/R/sctype_score_.R"))
			# DB file
			db_ = paste0(sctypedir,"/ScTypeDB_full.xlsx");
			
			# tissue can be e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
			db_read = openxlsx::read.xlsx(db_); tissues_ = unique(db_read$tissueType); result_ = c()
			if (tissue == "") {
				for(tissue in tissues_){
					print(paste0("Checking...", tissue));
					# prepare gene sets
					gs_list = gene_sets_prepare(db_, tissue);
					# prepare obj
					obj = as.matrix(SOB[["RNA"]]@scale.data)
					es.max = sctype_score(scRNAseqData = obj, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative, 
						marker_sensitivity = gs_list$marker_sensitivity, verbose=!0);
					cL_results = do.call("rbind", lapply(unique(SOB@meta.data$seurat_clusters), function(cl){
						es.max.cl = sort(rowSums(es.max[ ,rownames(SOB@meta.data[SOB@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
						head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl), 10)
					}))
					dt_out = cL_results %>% group_by(cluster) %>% top_n(n = 1)
					# return mean score for tissue
					result_ = rbind(result_, data.frame(tissue = tissue, score = mean(dt_out$scores)))
				}
				# order by mean score
				result_ = result_[order(-result_$score),]
				tissue = result_[1,1]
			}
			
			if (!(tissue %in% tissues_)) {
				stop(paste0("tissue not supported by sctype, must be one of: ",paste0(tissues_,collapse=",")))
			}

			# prepare gene sets
			gs_list = gene_sets_prepare(db_, tissue)
			
			# get cell-type by cell matrix
			es.max = sctype_score(scRNAseqData = SOB[["RNA"]]@scale.data, scaled = TRUE, 
				gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
			
			# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
			# In case Seurat is used, it is either SOB[["RNA"]]@scale.data (default), SOB[["SCT"]]@scale.data, in case sctransform is used for normalization,
			# or SOB[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.
			
			# merge by cluster
			cL_results = do.call("rbind", lapply(unique(SOB@meta.data$seurat_clusters), function(cl){
				es.max.cl = sort(rowSums(es.max[ ,rownames(SOB@meta.data[SOB@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
				head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(SOB@meta.data$seurat_clusters==cl)), 10)
			}))
			sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
			
			# set low-confident (low ScType score) clusters to "unknown"
			# sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
			# print(sctype_scores[,1:3])
			
			SOB@meta.data$customclassif = ""
			for(j in unique(sctype_scores$cluster)){
				cl_type = sctype_scores[sctype_scores$cluster==j,][1,]; 
				SOB@meta.data$customclassif[SOB@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
				if (as.numeric(as.character(cl_type$scores)) < cl_type$ncells/4 ) {
					SOB@meta.data$customclassif[SOB@meta.data$seurat_clusters == j] = paste(as.character(cl_type$type[1]),"(Unknown)")
				}
			}
			DimPlot(SOB, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')
			ggsave(umappng)
			qc_data <- FetchData(SOB, vars = c("ident", "UMAP_1", "UMAP_2"))
			out=data.frame(
				cell=rownames(SOB@meta.data),
				group=SOB@meta.data$customclassif,
				group_filtered=SOB@meta.data$customclassif,
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
			for(j in unique(sctype_scores$cluster)){
				cl_type = sctype_scores[sctype_scores$cluster==j,][1,]; 
				out$score[out$seurat_clusters == j] = cl_type$scores[1]
				out$ncells[out$seurat_clusters == j] = cl_type$ncells[1]
				if (as.numeric(as.character(cl_type$scores)) < cl_type$ncells/4 ) {
					out$group_filtered[out$seurat_clusters == j] = "Unknown"
				}
			}
			write.table(out,file=tempresult,sep="\t",quote=FALSE,row.names = FALSE)
		}
		file rename $tempresult $groupfile
	}
	sc_pseudobulk_job $scgenefile $scisoformfile $groupfile
}

proc cg_sc_celltyper_sctype {args} {
	set args [job_init {*}$args]
	sc_celltyper_sctype_job {*}$args
	job_wait
}
