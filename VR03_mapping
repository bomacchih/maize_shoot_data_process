library(Seurat)
library(Matrix)

# pick the assay that contains raw counts (usually "RNA")
DefaultAssay(pbmc) <- "RNA"

# 1) Filter genes detected in >= 3 cells  (equivalent to min.cells = 3)
cnt <- GetAssayData(pbmc, slot = "counts")
keep_genes <- Matrix::rowSums(cnt > 0) >= 3
pbmc <- pbmc[keep_genes, ]

# 2) Make sure per-cell QC columns exist
pbmc[["nFeature_RNA"]] <- Matrix::colSums(cnt[keep_genes, ] > 0)
pbmc[["nCount_RNA"]]   <- Matrix::colSums(cnt[keep_genes, ])

# 3) Filter cells with >= 200 detected genes (equivalent to min.features = 200)
pbmc <- subset(pbmc, subset = nFeature_RNA >= 200)

library(Seurat)

# 0) Make sure we're using the assay that has raw counts
DefaultAssay(pbmc) <- "RNA"   # change if your counts live elsewhere

# 1) Your mito gene list
mt.genes <- c("MT-ND1","MT-ND2","MT-CO1","MT-CO2","MT-ATP6")  # example

# 2) Keep only genes that exist in this object
mt.use <- intersect(mt.genes, rownames(pbmc))

if (length(mt.use) == 0) {
  stop("None of the provided mitochondrial genes were found in rownames(pbmc). 
       Check naming (symbols vs Ensembl IDs; case; species).")
}

# 3) Compute %MT
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, features = mt.use)

