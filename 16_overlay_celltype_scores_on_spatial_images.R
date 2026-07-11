library(Seurat)
library(ggplot2)
library(patchwork)

obj <- XGE202122_S5_subset_embleaf_harmony_join

# 1) Label High vs Low at threshold 0.1
obj$HighV_label <- ifelse(obj$SPOT_Vascular_tissue > 0.1, "High", "Low")
obj$HighV_label <- factor(obj$HighV_label, levels = c("Low","High"))

# 2) List images (slices)
imgs <- Images(obj)
stopifnot(length(imgs) > 0)

# 3) Helper: cells for one image
get_cells_for_image <- function(object, image) {
  rownames(GetTissueCoordinates(object = object, image = image))
}

# 4) Build per-image plots (plot only High cells; keep counts in title)
plots <- lapply(imgs, function(im) {
  cells_im <- get_cells_for_image(obj, im)

  # counts per image (ensure both levels shown even if zero)
  tab_im <- table(droplevels(obj$HighV_label[cells_im]))
  low_n  <- if ("Low"  %in% names(tab_im)) as.integer(tab_im["Low"])  else 0
  high_n <- if ("High" %in% names(tab_im)) as.integer(tab_im["High"]) else 0
  tot_n  <- length(cells_im)

  ttl <- paste0("SPOT_Vascular_tissue > 0.1 : ", im,
                "   [High = ", high_n, " / Low = ", low_n, " / Total = ", tot_n, "]")

  # plot only High cells
  high_im <- intersect(cells_im, colnames(obj)[obj$HighV_label == "High"])

  if (length(high_im) > 0) {
    p <- SpatialDimPlot(
      obj,
      images = im,
      cells  = high_im,
      cols   = "red",           # only high shown
      pt.size.factor = 1.6,
      combine = TRUE
    ) + ggtitle(ttl)
  } else {
    # no High: show background image only, with counts
    p <- ImagePlot(obj, images = im) + ggtitle(ttl)
  }
  p
})

# 5) Arrange 4×4 grids; paginate if >16
wrap_and_print <- function(pl, ncol = 4, nrow = 4) {
  n_per <- ncol * nrow
  n_pages <- ceiling(length(pl) / n_per)
  for (pg in seq_len(n_pages)) {
    idx <- (((pg - 1) * n_per) + 1):min(pg * n_per, length(pl))
    print(
      wrap_plots(pl[idx], ncol = ncol, nrow = nrow) +
        plot_annotation(title = paste0(
          "High-only overlays (threshold 0.1) — page ", pg, " of ", n_pages
        ))
    )
  }
}

wrap_and_print(plots, ncol = 4, nrow = 4)

# 6) (Optional) save one plot per page to a PDF
# pdf("High_only_SPOT_Vascular_overlays_by_image.pdf", width = 7.5, height = 7.5, useDingbats = FALSE)
# for (p in plots) print(p)
# dev.off()



library(Seurat)
library(ggplot2)
library(patchwork)

obj <- XGE202122_S5_subset_embleaf_harmony_join

# 1) Label High vs Low at threshold 0.01 (SAM)
obj$HighSAM_label <- ifelse(obj$SPOT_shoot_apical_meristem > 0.01, "High", "Low")
obj$HighSAM_label <- factor(obj$HighSAM_label, levels = c("Low","High"))

# 2) List images (slices)
imgs <- Images(obj)
stopifnot(length(imgs) > 0)

# 3) Helper: cells for one image
get_cells_for_image <- function(object, image) {
  rownames(GetTissueCoordinates(object = object, image = image))
}

# 4) Build per-image plots (plot only High cells; keep counts in title)
plots <- lapply(imgs, function(im) {
  cells_im <- get_cells_for_image(obj, im)

  # counts per image (ensure both levels shown even if zero)
  tab_im <- table(droplevels(obj$HighSAM_label[cells_im]))
  low_n  <- if ("Low"  %in% names(tab_im)) as.integer(tab_im["Low"])  else 0
  high_n <- if ("High" %in% names(tab_im)) as.integer(tab_im["High"]) else 0
  tot_n  <- length(cells_im)

  ttl <- paste0("SPOT_shoot_apical_meristem > 0.01 : ", im,
                "   [High = ", high_n, " / Low = ", low_n, " / Total = ", tot_n, "]")

  # plot only High cells
  high_im <- intersect(cells_im, colnames(obj)[obj$HighSAM_label == "High"])

  if (length(high_im) > 0) {
    p <- SpatialDimPlot(
      obj,
      images = im,
      cells  = high_im,
      cols   = "red",           # only high shown
      pt.size.factor = 1.6,
      combine = TRUE
    ) + ggtitle(ttl)
  } else {
    # no High: show background image only, with counts
    p <- ImagePlot(obj, images = im) + ggtitle(ttl)
  }
  p
})

# 5) Arrange 4×4 grids; paginate if >16
wrap_and_print <- function(pl, ncol = 4, nrow = 4) {
  n_per <- ncol * nrow
  n_pages <- ceiling(length(pl) / n_per)
  for (pg in seq_len(n_pages)) {
    idx <- (((pg - 1) * n_per) + 1):min(pg * n_per, length(pl))
    print(
      wrap_plots(pl[idx], ncol = ncol, nrow = nrow) +
        plot_annotation(title = paste0(
          "High-only SAM overlays (threshold 0.01) — page ", pg, " of ", n_pages
        ))
    )
  }
}

wrap_and_print(plots, ncol = 4, nrow = 4)

# 6) (Optional) save each slice on its own PDF page
# pdf("High_only_SPOT_SAM_overlays_by_image.pdf", width = 7.5, height = 7.5, useDingbats = FALSE)
# for (p in plots) print(p)
# dev.off()
