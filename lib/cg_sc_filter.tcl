suppressMessages(library(tidyverse))
suppressMessages(library(argparser))
suppressMessages(library(DropletUtils))
suppressMessages(library(DoubletFinder))
suppressMessages(library(scDblFinder))
suppressMessages(library(Seurat))

p <- arg_parser("sc_filter")
p <- add_argument(p, "--percent_mito", help = "The cut-off for percentage of mitochondrial reads per cell", type = "numeric", nargs = '1', default = "5")
p <- add_argument(p, "--input", help = "Path to directory containing raw files in 10X format", type = "character", nargs = '1')
p <- add_argument(p, "--transcriptsfile", help = "Path to directory containing transcripts data in 10X format", type = "character", nargs = '1')
p <- add_argument(p, "--metadatafile", help = "file with metadata (future)", type = "character", nargs = '1',default = "")
p <- add_argument(p, "--out10x", help = "name of (output) filtered gene file", type = "character", nargs = '+')
p <- add_argument(p, "--outinfofile", help = "name of (output) cell info file", type = "character", nargs = '+')
p <- add_argument(p, "--outrds", help = "name of (output) seurat rds file", type = "character", nargs = '+', default = "")
p <- add_argument(p, "--outkneeplot", help = "name of (output) kneeplot file", type = "character", nargs = '+', default = "")
p <- add_argument(p, "--nexpected", help = "number of expected cells", type = "numeric", nargs = '+', default = "3000")


arg <- parse_args(p)

input <- arg$input
workdir <- arg$workdir
transcriptsfile <- arg$transcriptsfile
metadatafile <- arg$metadatafile
out10x <- arg$out10x
outrds <- arg$outrds
outinfofile <- arg$outinfofile
outkneeplot <- arg$outkneeplot
nexpected <- arg$nexpected
# not yet
usescrublet <- 0

# 1. Load sample
mat=Read10X(input)

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

# 2. EmptyDrops
cat("running EmptyDrops\n")
mito.genes <- grep(pattern = "^MT-", x = rownames(mat), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", x = rownames(mat), value = TRUE)
qc <- c(mito.genes, RPS.genes)
keep <- !row.names(mat) %in% unlist(qc)

# not very clean, but on some data emptyDropsCellRanger gives an error (that I could not easily patch)
# so use emptyDrops as a fallback where this is the case
emptydrops.error=0
rm(e.out) ; rm(is.cell)
e.out <- tryCatch({
	emptyDropsCellRanger(mat[keep, ],n.expected.cells=nexpected)
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
		lower=sortedcols[nexpected]/10
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

SOB <- CreateSeuratObject(filt_mat)
rm(filt_mat)

SOB[["percent.mt"]] <- PercentageFeatureSet(SOB, pattern = "^MT-")

# 3. Doublet identification
# 3.1 scDblFinder
cat("running scDblFinder\n")
sce_tmp <- as.SingleCellExperiment(SOB)
sce_tmp <- scDblFinder(sce_tmp, dbr.sd = 1)
seurat_tmp <- as.Seurat(sce_tmp, counts = "counts", data = NULL)

SOB@meta.data <- seurat_tmp@meta.data
rm(seurat_tmp, sce_tmp)
SOB@meta.data$scDblFinder.class <- gsub("singlet", "Singlet", SOB@meta.data$scDblFinder.class)
SOB@meta.data$scDblFinder.class <- gsub("doublet", "Doublet", SOB@meta.data$scDblFinder.class)

# Intermezzo. Mitochondrial cut-off and standard processing (scDblFinder does not prefer this processing step first, DoubletFinder does)
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

# 3.2 DoubletFinder
cat("running DoubletFinder")
sweep.res.list_SOB <- paramSweep_v3(SOB, PCs = 1:50, sct = FALSE)
sweep.stats_SOB <- summarizeSweep(sweep.res.list_SOB, GT = FALSE)
bcmvn_SOB <- find.pK(sweep.stats_SOB)
optimal_pK <- as.numeric(as.character(bcmvn_SOB$pK[which.max(bcmvn_SOB$BCmetric)]))
homotypic.prop <- modelHomotypic(SOB$seurat_clusters)
# doesn't work -> which package?
# nExpDoublets <- estimate_doublet_rate(SOB)
expected_doublet_rate_per_10000_cells <- 0.075
total_cells <- ncol(SOB)
nExpDoublets <- (expected_doublet_rate_per_10000_cells * (total_cells))
nExpDoublets.adj <- round(nExpDoublets*(1-homotypic.prop))

SOB <- doubletFinder_v3(SOB, pN = 0.25, pK = optimal_pK, nExp = nExpDoublets, PCs = 1:50)
colnames(SOB@meta.data)[colnames(SOB@meta.data) == grep("^DF.classifications_", colnames(SOB@meta.data), value = TRUE)] <- "DoubletFinder"
SOB <- doubletFinder_v3(SOB, pN = 0.25, pK = optimal_pK, nExp = nExpDoublets.adj, PCs = 1:50)
colnames(SOB@meta.data)[colnames(SOB@meta.data) == grep("^DF.classifications_", colnames(SOB@meta.data), value = TRUE)] <- "DoubletFinder.adj"

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


cat(paste0("Quality control for ", input, "Completed!"), "\n")

if (metadatafile != "") {
	# 4. Add Metadata
	cat("Adding Metadata")
	sample_info <- fread(metadatafile)
	SOB@meta.data <- merge(SOB@meta.data, sample_info, by.x = "dnumber", by.y = "dnumber", all = TRUE)
	rm(sample_info)
}

# 5. Add transcript data
cat("Add transcript data")
transcripts <- Read10X(transcriptsfile, gene.column = 1)
barcodes_to_keep <- colnames(transcripts) %in% colnames(SOB)
transcripts <- transcripts[, barcodes_to_keep]
SOB[["ISOFORM"]] <- CreateAssayObject(counts = transcripts)
rm(transcripts)

cat("Metadata and isoforms added", "\n")

if (outrds != "") {
	# 6. Save preprocessed object
	cat(paste0("Save rds in ",outrds,"\n"))
	# name <- paste("SOB_processed_", unique(SOB$Sample), sep = "")
	saveRDS(object = SOB, file = outrds)
}

if (out10x != "") {
	cat("Writing 10x output to", out10x, "\n")
	write10xCounts(out10x,SOB_filtered@assays$RNA@counts,version="3")
}

cat("Writing info file", outinfofile, "\n")
info=transform(SOB@meta.data, cell=rownames(SOB@meta.data))
write.table(info,file=outinfofile,sep="\t",quote=FALSE,row.names = FALSE)
cat("is_cell_method\tis_cell_lower\n",is.cell.method,"\t",is.cell.lower,"\n",sep="",file=paste0(outinfofile,".analysisinfo"))

cat(paste0("Object saved, next!"), "\n")
		      

