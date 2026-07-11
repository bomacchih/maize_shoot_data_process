# Seurat v5
rna.pb   <- AggregateExpression(sc_merged, return.seurat = FALSE)$RNA
vis.pb   <- AverageExpression(XGE202122_merged_emleaf)$RNA
cor.mat  <- cor(rna.pb , vis.pb , method = "spearman")
pheatmap::pheatmap(cor.mat)
