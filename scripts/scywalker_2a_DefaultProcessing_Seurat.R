

p <- add_argument(p, "--normMethod", help = "Method to perform normalization and scaling. Choose between log-normalization and SCT", type = "character", nargs = '1')
p <- add_argument(p, "--workdir", help = "Path to working directory", type = "character", nargs = '+')

arg <- parse_args(p)


mat_filtered <- Read10X(arg$workdir, "/filtered_matrix_", i, "/")
SOB_filtered <- CreateSeuratObject(mat_filtered)

SOB_filtered <- NormalizeData(SOB_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
SOB_filtered <- FindVariableFeatures(SOB_filtered, selection.method = "vst", nfeatures = 2000)
SOB_filtered <- ScaleData(SOB_filtered)

SOB_filtered <- RunPCA(SOB_filtered, verbose = FALSE)
SOB_filtered <- FindNeighbors(SOB_filtered, reduction = "pca", dims = 1:50)
SOB_filtered <- FindClusters(SOB_filtered, resolution = 0.5)
SOB_filtered <- RunUMAP(SOB_filtered, dims = 1:50)

saveRDS(SOB_filtered, (arg$workdir, "/outdir/", i, "_processed.rds")

