library(Seurat)
library(SingleCellExperiment)

# seurat_sc: your Seurat v5 object with cell-type labels in seurat_sc$celltype
DefaultAssay(seurat_sc) <- "RNA"  # or your RNA-like assay

sce_ref <- as.SingleCellExperiment(seurat_sc)  # keeps counts & colData

# make sure cell-type labels are present for SPOTlight grouping
if (is.null(colData(sce_ref)$celltype)) {
  colData(sce_ref)$celltype <- seurat_sc$celltype
}

library(SpatialExperiment)
library(Seurat)

# seurat_vis: your Seurat v5 spatial object
#  - For Visium v1/v2 (spots): assay is often "RNA"
#  - For Visium HD: multiple binnings (e.g., "bins8", "bins16"); pick one

assay_name <- if ("RNA" %in% Assays(seurat_vis)) "RNA" else DefaultAssay(seurat_vis)
counts_mtx <- GetAssayData(seurat_vis, assay = assay_name, slot = "counts")

# Get 2D spot/bin coordinates
coords <- GetTissueCoordinates(seurat_vis)  # returns imagerow/imagecol by image
# If multiple images/capture areas, this already includes image context in rownames
# Make sure ordering matches the count matrix columns
coords <- coords[colnames(counts_mtx), , drop = FALSE]

# Build SPE
spe_vis <- SpatialExperiment(
  assays = list(counts = counts_mtx),
  colData = DataFrame(seurat_vis@meta.data[colnames(counts_mtx), , drop = FALSE]),
  spatialCoords = as.matrix(coords[, c("imagerow", "imagecol")])
)

sce_ref <- logNormCounts(sce_ref)

# Get all rownames
all_genes <- rownames(sce)

# Exclude ribosomal (RpL or RpS), mitochondrial, and chlorophyll genes
genes<- !(grepl("^Rp[l|s]", all_genes) |
                   all_genes %in% MT_genes_list |
                   all_genes %in% Pltd_genes_list)

dec <- modelGeneVar(sce_ref, subset.row = genes)
plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
curve(metadata(dec)$trend(x), col = "blue", add = TRUE)

hvg<-getTopHVGs(dec, n=3000)

colLabels(sce) <- colData(sce)$free_annotation

# Compute marker genes
mgs <- scoreMarkers(sce, subset.row = genes)

#Then we want to keep only those genes that are relevant for each cell identity:

mgs_fil <- lapply(names(mgs), function(i) {
    x <- mgs[[i]]
    # Filter and keep relevant marker genes, those with AUC > 0.8
    x <- x[x$mean.AUC > 0.8, ]
    # Sort the genes from highest to lowest weight
    x <- x[order(x$mean.AUC, decreasing = TRUE), ]
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i
    data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)
# split cell indices by identity
idx <- split(seq(ncol(sce)), sce$free_annotation)
# downsample to at most 20 per identity & subset
# We are using 5 here to speed up the process but set to 75-100 for your real
# life analysis
n_cells <- 100
cs_keep <- lapply(idx, function(i) {
    n <- length(i)
    if (n < n_cells)
        n_cells <- n
    sample(i, n_cells)
})
sc_ref <- sce_ref[, unlist(cs_keep)]

res <- SPOTlight(
    x = sce,
    y = spe,
    groups = as.character(sce$free_annotation),
    mgs = mgs_df,
    hvg = hvg,
    weight_id = "mean.AUC",
    group_id = "cluster",
    gene_id = "gene")

# Extract deconvolution matrix
head(mat <- res$mat)[, seq_len(3)]

# Extract NMF model fit
mod <- res$NMF


plotTopicProfiles(
    x = mod,
    y = sce$free_annotation,
    facet = FALSE,
    min_prop = 0.01,
    ncol = 1) +
    theme(aspect.ratio = 1)


plotTopicProfiles(
    x = mod,
    y = sce$free_annotation,
    facet = TRUE,
    min_prop = 0.01,
    ncol = 6)

library(Seurat)
# seurat_vis: your Seurat v5 spatial object
# mat: rows = barcodes/spots, cols = cell types (proportions 0–1)

# 1) align rows of mat to Seurat cells
common <- intersect(rownames(mat), colnames(seurat_vis))
stopifnot(length(common) > 0)
mat_aln <- mat[common, , drop = FALSE]
mat_aln <- mat_aln[match(colnames(seurat_vis), rownames(mat_aln)), , drop = FALSE]  # reorder to cells

# 2) make safe column names and add to meta.data
new_names <- paste0("SPOT_", make.names(colnames(mat_aln)))
colnames(mat_aln) <- new_names
seurat_vis <- AddMetaData(seurat_vis, metadata = mat_aln)

# 3) useful summaries
library(matrixStats)
seurat_vis$SPOT_top_type <- colnames(mat_aln)[max.col(mat_aln, ties.method = "first")]
seurat_vis$SPOT_top_prop <- rowMaxs(as.matrix(mat_aln))
seurat_vis$SPOT_entropy  <- -rowSums(mat_aln * log(mat_aln + 1e-9))

# 4) plot any proportion as a “feature”
SpatialFeaturePlot(seurat_vis, features = "SPOT_Mesophyll")  # replace with one of new_names



