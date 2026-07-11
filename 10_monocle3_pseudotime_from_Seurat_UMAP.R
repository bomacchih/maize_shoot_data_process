library(Seurat)
library(monocle3)

# Example Seurat object
seurat_obj <- YourSeuratObject

# Extract raw counts
counts <- GetAssayData(seurat_obj, assay="RNA", layer= "counts")

# Extract metadata
cell_metadata <- seurat_obj@meta.data

# Extract gene (feature) information
gene_metadata <- data.frame(gene_short_name = rownames(counts))
rownames(gene_metadata) <- rownames(counts)

# Extract UMAP embeddings
umap_coordinates <- Embeddings(seurat_obj, reduction = "umap")

# Create a Monocle3 cell_data_set object
cds <- new_cell_data_set(
  expression_data = counts,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

# Ensure UMAP embeddings and cell names match
if (!all(rownames(umap_coordinates) %in% colnames(cds))) {
  stop("Cell names in Seurat and Monocle3 do not match!")
}

# Add UMAP coordinates to the Monocle3 dataset
reducedDims(cds)$UMAP <- umap_coordinates

cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

plot_cells(cds, color_cells_by = "new.cluster.ids")

cds<-cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds)

plot_cells(cds,
           color_cells_by = "new.cluster.ids",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

cds<-order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=3, trajectory_graph_segment_size = 2)
