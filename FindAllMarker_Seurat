markers_SCT_all <- FindAllMarkers(
    object = XGE202122_S5_leaf_SCT_harmony, 
    only.pos = TRUE,               # Only return positive markers
    min.pct = 0.25,                # Minimum fraction of cells expressing the gene
    logfc.threshold = 0.25,        # Minimum log fold-change threshold
    test.use = "wilcox"   )   


significant_markers_RNA_all <- markers_RNA_all %>% filter(p_val_adj<0.05) %>% pull(gene) %>% unique()
