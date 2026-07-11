# sum of cells by domains
table(XGE202122_S5_subset_leaf_joinlayers$domains)
sum(table(XGE202122_S5_subset_leaf_joinlayers$domains)))

# sum of genes by domains
domains<-unique(XGE202122_S5_subset_leaf$domains)
> sapply(X=domains, function(c) {
+  cells.c<-WhichCells(XGE202122_S5_subset_leaf_joinlayers, idents=c)
+  nFeatures.c <-sum(rowSums(XGE202122_S5_subset_leaf_joinlayers[['RNA']]$counts[, cells.c]) !=0)
+  return(nFeatures.c)
+ }
+ )
