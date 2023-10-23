

p <- add_argument(p, "--percent_mito", help = "The cut-off for percentage of mitochondrial reads per cell", type = "numeric", nargs = '1', default = "5")
p <- add_argument(p, "--input", help = "Path to directory containing raw files in 10X format", type = "character", nargs = '1')
p <- add_argument(p, "--workdir", help = "Path to working directory", type = "character", nargs = '+')

arg <- parse_args(p)

input_df <- arg$input
i <- input_df$samplename


# 1. Load sample
cat(paste0("Loading Sample ", i), "\n")
file_path <- paste0("~/cgsc_samples_converted_to_10x/", i, "/")
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

e.out <- emptyDropsCellRanger(mat[keep, ])
is.cell <- e.out$FDR <= 0.001
is.cell[is.na(is.cell)] <- FALSE
filt_mat <- mat[, is.cell]
rm(mat, barcodes, features, e.out, is.cell)

SOB <- CreateSeuratObject(filt_mat)
rm(filt_mat)

SOB[["percent.mt"]] <- PercentageFeatureSet(SOB, pattern = "^MT-")

# 3. Doublet identification
# 3.1 scDblFinder
sce_tmp <- as.SingleCellExperiment(SOB)
sce_tmp <- scDblFinder(sce_tmp, dbr.sd = 1)
seurat_tmp <- as.Seurat(sce_tmp, counts = "counts", data = NULL)

SOB@meta.data <- seurat_tmp@meta.data
rm(seurat_tmp, sce_tmp)
SOB@meta.data$scDblFinder.class <- gsub("singlet", "Singlet", SOB@meta.data$scDblFinder.class)
SOB@meta.data$scDblFinder.class <- gsub("doublet", "Doublet", SOB@meta.data$scDblFinder.class)

# Intermezzo. Mitochondrial cut-off and standard processing (scDblFinder does not prefer this processing step first, DoubletFinder does)
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
sweep.res.list_SOB <- paramSweep_v3(SOB, PCs = 1:50, sct = FALSE)
sweep.stats_SOB <- summarizeSweep(sweep.res.list_SOB, GT = FALSE)
bcmvn_SOB <- find.pK(sweep.stats_SOB)
optimal_pK <- as.numeric(as.character(bcmvn_SOB$pK[which.max(bcmvn_SOB$BCmetric)]))
homotypic.prop <- modelHomotypic(SOB$seurat_clusters)
nExpDoublets <- estimate_doublet_rate(SOB)
nExpDoublets.adj <- round(nExpDoublets*(1-homotypic.prop))

SOB <- doubletFinder_v3(SOB, pN = 0.25, pK = optimal_pK, nExp = nExpDoublets, PCs = 1:50)
colnames(SOB@meta.data)[colnames(SOB@meta.data) == grep("^DF.classifications_", colnames(SOB@meta.data), value = TRUE)] <- "DoubletFinder"
SOB <- doubletFinder_v3(SOB, pN = 0.25, pK = optimal_pK, nExp = nExpDoublets.adj, PCs = 1:50)
colnames(SOB@meta.data)[colnames(SOB@meta.data) == grep("^DF.classifications_", colnames(SOB@meta.data), value = TRUE)] <- "DoubletFinder.adj"
						      
# 3.3 Scrublet
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

cat(paste0("Quality control for ", i, "Completed!"), "\n")

# 4. Add Metadata
sample_info <- fread("/location/to/metadata.tsv")
SOB@meta.data <- merge(SOB@meta.data, sample_info, by.x = "dnumber", by.y = "dnumber", all = TRUE)
rm(sample_info)

# 5. Add transcript data
transcripts <- Read10X(paste0(sub("/$", "", file_path), "_counts_weighed/"), gene.column = 1)
barcodes_to_keep <- colnames(transcripts) %in% colnames(SOB)
transcripts <- transcripts[, barcodes_to_keep]
SOB[["ISOFORM"]] <- CreateAssayObject(counts = transcripts)
rm(transcripts)

cat(paste0("Metadata and isoforms added!"), "\n")

# 6. Save preprocessed object
name <- paste("SOB_processed_", unique(SOB$Sample), sep = "")
saveRDS(object = SOB, file = paste0(file_path, name, ".rds"))

write10xCounts(SOB_filtered, (arg$workdir, "/filtered_matrix_", i, "/")

cat(paste0("Object saved, next!"), "\n")
				      
