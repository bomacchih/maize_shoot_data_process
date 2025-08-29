library(Seurat)
library(patchwork)
library(dplyr)

# ---- genes to plot ----
genes_named <- c(
  ZmGRF6  = "Zm00001eb251820",
  ZmGRF8  = "Zm00001eb071830",
  ZmGRF10 = "Zm00001eb193180",
  ZmARF20 = "Zm00001eb232120",
  ZmARF28 = "Zm00001eb408800"
)

# ---- domains live in obj$domains ----
dom_order <- c("SAM","P1_P2","P3","P4","P5")
obj$domains <- factor(as.character(obj$domains), levels = dom_order)

# ---- helpers ----
pick_assay_for_gene <- function(object, gene) {
  if ("RNA" %in% Assays(object) && gene %in% rownames(object[["RNA"]])) return("RNA")
  if ("SCT" %in% Assays(object) && gene %in% rownames(object[["SCT"]])) return("SCT")
  stop(sprintf("Gene %s not found in RNA or SCT assays.", gene))
}

feat_plot_one <- function(object, gene_id, title_txt) {
  a <- pick_assay_for_gene(object, gene_id)
  DefaultAssay(object) <- a
  FeaturePlot(
    object, features = gene_id, reduction = "umap",
    slot = "data", cols = c("grey90","red"),
    min.cutoff = "q02", max.cutoff = "q98",
    order = TRUE, pt.size = 0.5
  ) + ggtitle(title_txt) + theme(legend.position = "none")
}

vln_plot_one <- function(object, gene_id, title_txt) {
  a <- pick_assay_for_gene(object, gene_id)
  DefaultAssay(object) <- a
  Idents(object) <- "domains"
  VlnPlot(
    object, features = gene_id, group.by = "domains", pt.size = 0,
    cols = rep("grey70", length(levels(obj$domains)))
  ) +
    ggtitle(title_txt) + xlab("") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1)
    )
}

# ---- build and arrange 2×5 ----
umaps <- mapply(
  function(nm, gid) feat_plot_one(obj, gid, paste0(nm, " (", gid, ")")),
  names(genes_named), genes_named, SIMPLIFY = FALSE
)
violins <- mapply(
  function(nm, gid) vln_plot_one(obj, gid, nm),
  names(genes_named), genes_named, SIMPLIFY = FALSE
)

wrap_plots(c(umaps, violins), nrow = 2, ncol = 5) +
  plot_annotation(title = "UMAP (top, red) and domain-ordered violin (bottom)")
