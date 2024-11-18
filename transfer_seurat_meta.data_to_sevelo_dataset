# Extract metadata in R
metadata <- seurat_obj@meta.data

# Select the desired property (e.g., "cell_type")
property <- metadata$cell_type

# Ensure rownames of metadata match cell barcodes
property_data <- data.frame(Barcode = rownames(metadata), Cell_Type = property)

# Save to CSV
write.csv(property_data, "seurat_metadata.csv", row.names = FALSE)


# in Python  
# Load the Seurat metadata
cell_type = pd.read_csv("seurat_metadata.csv")
cell_type_filter=cell_type.loc[inter_cell]
 x.obs['cell_type']=cell_type_filter.loc[x.obs_names, 'Cell_Type']
