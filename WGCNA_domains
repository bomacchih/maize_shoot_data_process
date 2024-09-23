DefaultAssay(XGE202122_S5_leaf_SCT_harmony)<-"SCT"
XGE202122_leaf_SCT_averaged_domains<-AverageExpression(XGE202122_S5_leaf_SCT_harmony, assay="SCT", slot="data")
data_matrix_SCT_averaged_domsins<-as.data.frame(XGE202122_leaf_SCT_averaged_domains$SCT)
data_matrix_SCT_averaged_domsins<-data_matrix_SCT_averaged_domsins[rowSums(data_matrix_SCT_averaged_domsins) != 0, ]
WGCNA_Tree = hclust(dist(data_matrix_SCT_averaged_domsins), method = "complete")
