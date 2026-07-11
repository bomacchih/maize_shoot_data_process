library(Seurat)
library(ggplot2)
library(patchwork)

obj <- XGE2021XGE202122_S5_subset_embleaf_harmony_join

## 1) Make a clean 2-level factor in meta.data
obj$HighV_label <- ifelse(obj$SPOT_Vascular_tissue > 0.7, "High", "Low")
obj$HighV_label <- factor(obj$HighV_label, levels = c("Low", "High"))
cols_named <- c("Low" = "grey80", "High" = "red")

## 2) List images (slices)
imgs <- Images(obj)
if (length(imgs) == 0) stop("No spatial images found in the object.")

## 3) Helper: cells for one image (via tissue coordinates)
get_cells_for_image <- function(object, image) {
  tc <- GetTissueCoordinates(object = object, image = image)
  rownames(tc)  # cell/spot names for that image
}

## 4) Build plots per image with safe color mapping
plots <- lapply(imgs, function(im) {
  cells_im <- get_cells_for_image(obj, im)
  # Present levels in this image
  levs_present <- levels(droplevels(obj$HighV_label[cells_im]))
  cols_use <- cols_named[levs_present]  # only colors for present levels

  # Counts for title (useful QC)
  tab_im <- table(droplevels(obj$HighV_label[cells_im]))
  tab_str <- paste(paste0(names(tab_im), "=", as.integer(tab_im)), collapse = "  ")

  p <- SpatialDimPlot(
    obj,
    images = im,
    group.by = "HighV_label",
    pt.size.factor = 1.6,
    combine = TRUE
  ) +
    scale_color_manual(values = cols_use, drop = TRUE) +
    ggtitle(paste0("High SPOT_Vascular_tissue (>0.7): ", im, "   [", tab_str, "]"))

  p
})

## 5) Arrange in 4×4 grids; paginate if >16
wrap_and_print <- function(pl, ncol = 4, nrow = 4) {
  n_per <- ncol * nrow
  n_pages <- ceiling(length(pl) / n_per)
  for (pg in seq_len(n_pages)) {
    idx <- (((pg - 1) * n_per) + 1):min(pg * n_per, length(pl))
    print(
      wrap_plots(pl[idx], ncol = ncol, nrow = nrow) +
        plot_annotation(title = paste0(
          "High vascular spots across slices (page ", pg, "/", n_pages, ")"
        ))
    )
  }
}

wrap_and_print(plots, ncol = 4, nrow = 4)

## 6) (Optional) save one plot per page to a PDF
# pdf("High_SPOT_Vascular_tissue_overlays_by_image.pdf", width = 7.5, height = 7.5, useDingbats = FALSE)
# for (p in plots) print(p)
# dev.off()
