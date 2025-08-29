library(Seurat)
library(patchwork)

gene <- "Zm00001eb251820"

# pick an assay that actually contains the gene
assay_to_use <- NA
if ("RNA" %in% Assays(obj) && gene %in% rownames(obj[["RNA"]])) assay_to_use <- "RNA"
if (is.na(assay_to_use) && "SCT" %in% Assays(obj) && gene %in% rownames(obj[["SCT"]])) assay_to_use <- "SCT"
if (is.na(assay_to_use)) stop(paste("Gene", gene, "not found in RNA or SCT assays."))

DefaultAssay(obj) <- assay_to_use

# global cutoffs (same scale for every sample)
expr_vec <- as.numeric(GetAssayData(obj, assay = assay_to_use, slot = "data")[gene, ])
vmin <- as.numeric(quantile(expr_vec, 0.02, na.rm = TRUE))
vmax <- as.numeric(quantile(expr_vec, 0.98, na.rm = TRUE))

# pick the first 14 spatial images (adjust order if needed)
imgs <- Images(obj)
if (length(imgs) < 14) stop("Object has fewer than 14 spatial images.")
imgs <- imgs[1:14]

# one plot per image
plots <- lapply(imgs, function(im) {
  SpatialFeaturePlot(
    object = obj,
    features = gene,
    images = im,
    slot = "data",
    min.cutoff = vmin,
    max.cutoff = vmax,
    pt.size.factor = 1.6
  ) + ggtitle(im)
})

# pad to 16 to make a 4x4 grid
if (length(plots) < 16) {
  plots <- c(plots, replicate(16 - length(plots), patchwork::plot_spacer(), simplify = FALSE))
}

wrap_plots(plots, ncol = 4) +
  plot_annotation(title = paste0(gene, " (GRF6) across 14 samples"))
