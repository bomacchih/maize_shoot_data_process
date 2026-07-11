# Maize shoot spatial transcriptomics data analysis

This repository contains example scripts used for processing and analyzing serial 10x Visium spatial transcriptomes of maize seedling shoot apices and developing leaf primordia. The scripts accompany the Bio-protocol manuscript:

**Serial Cryosection for the Spatial Transcriptomes of Maize Seedling Leaf**

The workflow was developed for serial sections of maize shoot apices generated using 10x Visium V1 Gene Expression slides. The analysis supports section splitting, spatial-domain annotation, pseudobulk comparison, developmental trend analysis, RNA velocity, pseudotime inference, scRNA-seq integration, cell-type mapping, and 3D visualization.

> **Note**  
> Most scripts are provided as analysis examples rather than fully automated command-line pipelines. File paths, object names, sample names, gene annotations, and plotting parameters should be adjusted for each dataset. The scripts were tested mainly with Visium V1 data. Application to Visium HD data may require modification, especially for bin-level data structures, image handling, and coordinate processing.

---

## Recommended analysis workflow

1. Process raw sequencing data using **Space Ranger**.
2. Import Space Ranger outputs into **R** using **STUtility** and **Seurat**.
3. Split multiple serial sections from the same capture area.
4. Perform quality control and remove low-quality spots or genes.
5. Normalize, integrate, and cluster spatial transcriptome datasets.
6. Annotate anatomical domains.
7. Identify marker genes and tissue supergroups.
8. Generate pseudobulk profiles for structural domains.
9. Analyze developmental gene-expression trends and GO enrichment.
10. Perform RNA velocity analysis.
11. Perform Monocle 3 pseudotime analysis.
12. Construct domain-ordered transcription factor co-expression networks.
13. Integrate Visium data with maize scRNA-seq reference data.
14. Assign scRNA-seq cell identities using SCINA.
15. Map scRNA-seq cell-type information onto Visium spots using Seurat anchors and SPOTlight.
16. Reconstruct serial sections in 3D.

---

## Repository guide and suggested pipeline titles

| Current file name | Suggested pipeline title | Suggested new file name | Main purpose |
|---|---|---|---|
| `general process and 3D` | Import, section splitting, STUtility processing, Harmony integration, and 3D visualization | `01_import_crop_integrate_3D_STUtility_Seurat.R` | Imports Space Ranger output, manually annotates multiple sections in one capture area, crops sections, masks and aligns images, performs QC, SCTransform, Harmony integration, clustering, marker detection, and 3D visualization. |
| `intergration` | Seurat normalization and dataset integration | `02_seurat_normalization_integration_clustering.R` | Performs SCTransform, PCA, integration, clustering, domain/sample factor ordering, and UMAP plotting. |
| `FindAllMarker_Seurat` | Marker-gene identification | `03_find_marker_genes_Seurat.R` | Identifies cluster marker genes using `FindAllMarkers()` and extracts significant genes. |
| `export_counts_by_domains` | Export domain-level expression matrices | `04_export_domain_expression_matrices.R` | Exports or prepares expression values grouped by anatomical domain for downstream analyses. |
| `pseudo-bulk-comparsion` | Pseudobulk comparison between spatial and scRNA-seq data | `05_pseudobulk_domain_scRNA_comparison.R` | Aggregates expression profiles and calculates Spearman correlations between Visium and scRNA-seq pseudobulk profiles. |
| `clustering gene expression` | Developmental gene-expression trend clustering | `06_developmental_expression_trend_clustering.R` | Orders domains, scales gene expression, clusters genes by developmental expression trends, plots median trends, and selects representative genes. |
| `WGCNA_domains` | Domain-level co-expression analysis | `07_domain_coexpression_WGCNA_input.R` | Generates domain-averaged expression matrices and hierarchical clustering for co-expression or TO-GCN-related analysis. |
| `scvelo_data_computing` | RNA velocity analysis using velocyto and scVelo | `08_RNA_velocity_velocyto_scVelo.R` | Runs `velocyto`, reads loom files, merges spliced/unspliced matrices with Seurat-derived metadata, and performs RNA velocity analysis using scVelo. |
| `RNA_velocity_celldancer` | Alternative RNA velocity analysis using cellDancer | `09_RNA_velocity_cellDancer.R` | Provides an alternative RNA velocity workflow using cellDancer-related processing. |
| `monocle3 with preprocessed UMPA` | Monocle 3 pseudotime analysis using Seurat UMAP | `10_monocle3_pseudotime_from_Seurat_UMAP.R` | Converts a Seurat object to a Monocle 3 `cell_data_set`, transfers UMAP coordinates, learns a trajectory graph, orders cells/spots, and plots pseudotime. |
| `transfer_seurat_metadata_to_monocle3` | Transfer Seurat metadata to Monocle 3 | `11_transfer_Seurat_metadata_to_Monocle3.R` | Transfers Seurat UMAP coordinates and metadata, such as cell type or domain labels, into a Monocle 3 object. |
| `maize_scRNA_integration` | Maize scRNA-seq reference integration | `12_maize_scRNA_reference_integration_Seurat.R` | Processes and integrates maize scRNA-seq libraries using Seurat workflows, including PCA, Harmony/CCA-style integration, clustering, UMAP, and layer joining. |
| `maize scRNA celltype_assign` | SCINA-based scRNA-seq cell-type assignment | `13_scRNA_celltype_annotation_SCINA.R` | Builds marker lists, runs SCINA, adds predicted labels and posterior probabilities to Seurat metadata, and visualizes cell-type annotations. |
| `SPOTlight` | SPOTlight cell-type deconvolution | `14_SPOTlight_celltype_deconvolution.R` | Converts Seurat objects to SingleCellExperiment/SpatialExperiment objects, selects marker genes, runs SPOTlight, adds estimated cell-type proportions to Visium metadata, and visualizes spatial proportions. |
| `Plot_specific_celltype_on_visium_UMAP` | Plot selected cell types on Visium UMAP | `15_plot_selected_celltypes_on_Visium_UMAP.R` | Visualizes selected predicted or deconvolved cell types on the Visium UMAP. |
| `SPOT_overlay_image` | Overlay SPOTlight results on tissue images | `16_overlay_celltype_scores_on_spatial_images.R` | Plots predicted cell-type proportions or SPOTlight scores on spatial tissue images. |
| `gene_spatialFeaturePlot_overlay_images` | Overlay gene expression on spatial images | `17_overlay_gene_expression_on_spatial_images.R` | Generates spatial feature plots for selected genes overlaid on tissue images. |
| `VR03_mapping` | Example mapping for one representative Visium sample | `18_example_VR03_spatial_mapping.R` | Provides a sample-specific example for mapping or visualization using the VR03 dataset. |
| `SAM_V_UMAP_Violin_plot` | SAM and vascular marker visualization | `19_SAM_vein_marker_UMAP_violin_plots.R` | Generates UMAP and violin plots for selected SAM- or vascular-associated marker genes. |
| `along middle to margin of p5` | P5 medial-to-margin expression analysis | `20_P5_middle_to_margin_expression_analysis.R` | Examines gene-expression patterns across the middle-to-margin axis of P5. |
| `RNAread_normalized` | RNA read normalization | `21_RNA_read_normalization.R` | Normalizes read-count or expression matrices for downstream comparison. |
| `HDL` | Additional RNA velocity or high-dimensional analysis notes | `22_additional_velocity_or_high_dimensional_workflow.R` | Contains additional processing commands related to velocity or high-dimensional data handling. Rename after confirming final purpose. |
| `mis` | Miscellaneous exploratory analysis | `99_miscellaneous_exploratory_analysis.R` | Contains exploratory or temporary analyses. Consider splitting useful parts into the numbered pipelines above. |
| `Sv_phenotype_ANOVA` | Setaria phenotype statistical analysis | `validation_Setaria_phenotype_ANOVA.R` | Performs statistical analysis for Setaria mutant phenotype validation. This is related to biological validation rather than the core Visium processing workflow. |

---

## Suggested folder organization

A cleaner repository structure could be:

```text
maize_shoot_data_process/
├── README.md
├── scripts/
│   ├── 01_import_crop_integrate_3D_STUtility_Seurat.R
│   ├── 02_seurat_normalization_integration_clustering.R
│   ├── 03_find_marker_genes_Seurat.R
│   ├── 04_export_domain_expression_matrices.R
│   ├── 05_pseudobulk_domain_scRNA_comparison.R
│   ├── 06_developmental_expression_trend_clustering.R
│   ├── 07_domain_coexpression_WGCNA_input.R
│   ├── 08_RNA_velocity_velocyto_scVelo.R
│   ├── 09_RNA_velocity_cellDancer.R
│   ├── 10_monocle3_pseudotime_from_Seurat_UMAP.R
│   ├── 11_transfer_Seurat_metadata_to_Monocle3.R
│   ├── 12_maize_scRNA_reference_integration_Seurat.R
│   ├── 13_scRNA_celltype_annotation_SCINA.R
│   ├── 14_SPOTlight_celltype_deconvolution.R
│   ├── 15_plot_selected_celltypes_on_Visium_UMAP.R
│   ├── 16_overlay_celltype_scores_on_spatial_images.R
│   ├── 17_overlay_gene_expression_on_spatial_images.R
│   ├── 18_example_VR03_spatial_mapping.R
│   ├── 19_SAM_vein_marker_UMAP_violin_plots.R
│   ├── 20_P5_middle_to_margin_expression_analysis.R
│   └── 21_RNA_read_normalization.R
├── validation/
│   └── validation_Setaria_phenotype_ANOVA.R
├── metadata/
│   ├── sample_table.csv
│   ├── domain_annotation.csv
│   └── marker_gene_lists/
├── references/
│   ├── maize_B73_v5_reference/
│   └── transcription_factor_gene_lists/
└── results/
    ├── figures/
    ├── tables/
    └── intermediate_objects/
```

---

## Input data expected by the workflow

### Space Ranger output

Each Visium capture area should contain standard Space Ranger outputs, including:

```text
spaceranger_count_output/
├── filtered_feature_bc_matrix.h5
├── spatial/
│   ├── tissue_hires_image.png
│   ├── tissue_lowres_image.png
│   ├── scalefactors_json.json
│   └── tissue_positions.csv or tissue_positions_list.csv
├── possorted_genome_bam.bam
└── web_summary.html
```

### Metadata

A sample metadata table is recommended. Example columns:

| Column | Description |
|---|---|
| `sample_id` | Biological replicate or library ID |
| `capture_area` | Visium capture area ID |
| `section_id` | Section number after section splitting |
| `section_order` | Serial order along the proximal–distal axis |
| `domain` | Anatomical domain label |
| `image_path` | Path to tissue image |
| `matrix_path` | Path to feature-barcode matrix |
| `scalefactor_path` | Path to scalefactors JSON file |

---

## Main analysis modules

### 1. Import, crop, and reconstruct serial sections

Use `01_import_crop_integrate_3D_STUtility_Seurat.R`.

This script imports Space Ranger output using STUtility, manually labels multiple serial sections within one capture area, defines crop windows, crops each section, removes spots from other sections, merges cropped datasets, masks images, aligns serial sections, and prepares a 3D stack.

Main functions used:

- `InputFromTable()`
- `ManualAnnotation()`
- `CropImages()`
- `SubsetSTData()`
- `MergeSTData()`
- `MaskImages()`
- `ManualAlignImages()`
- `FeaturePlot3D()`

### 2. Quality control, normalization, and integration

Use `02_seurat_normalization_integration_clustering.R`.

This script performs quality control, SCTransform normalization, PCA, integration, UMAP visualization, and graph-based clustering.

Recommended checks:

- Number of detected genes per spot
- Total UMI counts per spot
- Mitochondrial read percentage
- Chloroplast read percentage
- Section-specific or replicate-specific batch effects

### 3. Structural-domain annotation

Structural-domain labels should be assigned using histology together with gene-expression patterns. In the maize shoot apex dataset, the main domains are:

- `SAM`
- `P1_P2`
- `P3`
- `P4`
- `P5`
- `coleoptile`
- `co_v`

Spots overlapping more than one anatomical domain should be excluded from domain-level analyses.

### 4. Marker genes and tissue supergroups

Use `03_find_marker_genes_Seurat.R`.

This script identifies cluster marker genes using `FindAllMarkers()`. Marker genes can be used to validate domain annotations and define broader tissue supergroups.

### 5. Pseudobulk comparison

Use `05_pseudobulk_domain_scRNA_comparison.R`.

This script aggregates expression by domain or cell type and compares pseudobulk profiles. It can be used to compare Visium domains with scRNA-seq cell-type profiles.

### 6. Developmental trend analysis

Use `06_developmental_expression_trend_clustering.R`.

This script orders domains according to developmental progression, scales gene expression, clusters genes by expression trajectory, plots median trends, and identifies representative genes.

Example developmental order:

```text
SAM → P1_P2 → P3 → P4 → P5
```

For broader comparisons, coleoptile and coleoptile vein can also be included.

### 7. RNA velocity

Use `08_RNA_velocity_velocyto_scVelo.R`.

This workflow runs `velocyto` on Space Ranger output, loads loom files, matches barcodes and gene names with Seurat-derived metadata, merges spliced/unspliced matrices, and performs velocity analysis with scVelo.

Main tools:

- `velocyto run10x`
- `scanpy`
- `scvelo`
- `cellrank`

### 8. Pseudotime analysis

Use `10_monocle3_pseudotime_from_Seurat_UMAP.R`.

This workflow converts a Seurat object into a Monocle 3 `cell_data_set`, transfers Seurat UMAP embeddings, learns a trajectory graph, and assigns pseudotime. SAM spots should be used as the root state for maize shoot developmental ordering.

### 9. scRNA-seq reference integration

Use `12_maize_scRNA_reference_integration_Seurat.R`.

This workflow processes maize shoot scRNA-seq data and generates a reference object for label transfer or deconvolution.

Recommended steps:

- Create Seurat objects from raw count matrices
- Remove low-quality cells and doublets
- Normalize with SCTransform
- Integrate libraries
- Cluster and annotate cell types

### 10. SCINA cell-type assignment

Use `13_scRNA_celltype_annotation_SCINA.R`.

This script uses marker gene lists to assign cell identities in the scRNA-seq reference. It removes duplicated markers, filters marker lists, runs SCINA, adds posterior probabilities, and visualizes annotated cell types.

### 11. SPOTlight cell-type deconvolution

Use `14_SPOTlight_celltype_deconvolution.R`.

This workflow estimates cell-type proportions in Visium spots using the annotated scRNA-seq object as the reference. Estimated proportions are added back to the Visium Seurat metadata and can be visualized by `SpatialFeaturePlot()` or UMAP.

### 12. 3D visualization

Use `01_import_crop_integrate_3D_STUtility_Seurat.R` and downstream visualization scripts.

Serial sections are aligned and stacked along the z-axis. Gene-expression patterns can be visualized in 3D using STUtility-derived functions or exported for Napari-based visualization.

---

## Software

The main analysis uses:

- Space Ranger
- R
- Seurat
- STUtility
- Harmony
- Monocle 3
- SCINA
- SPOTlight
- SingleCellExperiment
- SpatialExperiment
- Python
- scanpy
- scVelo
- cellrank
- velocyto
- napari

---

## Notes for users

- Check and revise hard-coded file paths before running scripts.
- Keep object names consistent across scripts.
- Confirm that barcode naming is consistent between Seurat objects, Space Ranger output, loom files, and exported metadata.
- For multi-section capture areas, section splitting is required before section-specific analysis and 3D reconstruction.
- Domain annotation should not rely only on clustering. Histological position and known marker genes should also be used.
- Visium HD data structures differ from Visium V1 data; additional testing and modification are required before applying these scripts to Visium HD.

---

## Citation

If you use this workflow, please cite:

Wu, C.-C. et al. Serial Spatial Transcriptomes Reveal Regulatory Transitions in Maize Leaf Development. *Plant Biotechnology Journal* (2026). DOI: 10.1111/pbi.70515

Also cite the accompanying Bio-protocol manuscript when available.
