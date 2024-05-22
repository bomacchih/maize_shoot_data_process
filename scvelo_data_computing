#!/usr/bin/bash
velocyto run10x /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/XGE202122_short_MTCL/XGE20_UL01_B73V5S210_short_MTCL /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/Zea_mays.Zm-B73-REFERENCE_NAM-5.0.51.gtf
velocyto run10x /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/XGE202122_short_MTCL/XGE20_UL02_B73V5S210_short_MTCL /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/Zea_mays.Zm-B73-REFERENCE_NAM-5.0.51.gtf
velocyto run10x /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/XGE202122_short_MTCL/XGE20_UL03_B73V5S210_short_MTCL /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/Zea_mays.Zm-B73-REFERENCE_NAM-5.0.51.gtf
velocyto run10x /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/XGE202122_short_MTCL/XGE20_UL04_B73V5S210_short_MTCL /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/Zea_mays.Zm-B73-REFERENCE_NAM-5.0.51.gtf
velocyto run10x /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/XGE202122_short_MTCL/XGE21_VR01_B73V5S210_short_MTCL /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/Zea_mays.Zm-B73-REFERENCE_NAM-5.0.51.gtf
velocyto run10x /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/XGE202122_short_MTCL/XGE21_VR02_B73V5S210_short_MTCL /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/Zea_mays.Zm-B73-REFERENCE_NAM-5.0.51.gtf
velocyto run10x /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/XGE202122_short_MTCL/XGE21_VR03_B73V5S210_short_MTCL /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/Zea_mays.Zm-B73-REFERENCE_NAM-5.0.51.gtf
velocyto run10x /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/XGE202122_short_MTCL/XGE21_VR04_B73V5S210_short_MTCL /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/Zea_mays.Zm-B73-REFERENCE_NAM-5.0.51.gtf
velocyto run10x /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/XGE202122_short_MTCL/XGE22_DQ01_B73V5S210_short_MTCL /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/Zea_mays.Zm-B73-REFERENCE_NAM-5.0.51.gtf
velocyto run10x /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/XGE202122_short_MTCL/XGE22_DQ02_B73V5S210_short_MTCL /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/Zea_mays.Zm-B73-REFERENCE_NAM-5.0.51.gtf
velocyto run10x /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/XGE202122_short_MTCL/XGE22_DQ03_B73V5S210_short_MTCL /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/Zea_mays.Zm-B73-REFERENCE_NAM-5.0.51.gtf
velocyto run10x /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/XGE202122_short_MTCL/XGE22_DQ04_B73V5S210_short_MTCL /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/Zea_mays.Zm-B73-REFERENCE_NAM-5.0.51.gtf
velocyto run10x /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/XGE202122_short_MTCL/XGE22_DQ05_B73V5S210_short_MTCL /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/Zea_mays.Zm-B73-REFERENCE_NAM-5.0.51.gtf
velocyto run10x /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/XGE202122_short_MTCL/XGE22_DQ06_B73V5S210_short_MTCL /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/Zea_mays.Zm-B73-REFERENCE_NAM-5.0.51.gtf
velocyto run10x /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/XGE202122_short_MTCL/XGE22_DQ07_B73V5S210_short_MTCL /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/Zea_mays.Zm-B73-REFERENCE_NAM-5.0.51.gtf
velocyto run10x /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/XGE202122_short_MTCL/XGE22_DQ08_B73V5S210_short_MTCL /home/macchihlee/ST/spaceranger/spaceranger-2.1.0/Zea_mays.Zm-B73-REFERENCE_NAM-5.0.51.gtf


## in Python command
import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd
import scvelo as scv
import cellrank as cr
import anndata as ad



# load sparse matrix
X_XGE22=io.mmread("/home/macchihlee/R/XGE20_21_22/XGE22_combined.counts_matrix.mtx")

# create anndata object
adata_XGE22=anndata.AnnData(X=X_XGE22.transpose().tocsr())

# create cell metadata
cell_meta_XGE22=pd.read_csv("/home/macchihlee/R/XGE20_21_22/XGE22_combined.metadata.csv")

# load gene names
with open ("/home/macchihlee/R/XGE20_21_22/XGE22_combined.gene_names.csv", 'r') as f:gene_names=f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata_XGE22.obs=cell_meta_XGE22
adata_XGE22.obs.index=adata_XGE22.obs['barcode']
adata_XGE22.var.index=gene_names

# load dimensional reduction
pca_XGE22=pd.read_csv("/home/macchihlee/R/XGE20_21_22/XGE22_combined.pca.csv")
pca_XGE22.index=adata_XGE22.obs.index

# set pca and umap
adata_XGE22.obsm['X_pca']=pca_XGE22.to_numpy()
adata_XGE22.obsm['X_umap']=np.vstack((adata_XGE22.obs['UMAP_1'].to_numpy(), adata_XGE22.obs['UMAP_2'].to_numpy())).T

# plot a UMAP colored by domains to test
sc.pl.umap(adata_XGE22, color=['domains'], frameon=False, save=True)

scv.settings.verbosity=3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity=2

$ load loom files for spliced/unspliced matrices for each sample
ldata_UL01=scv.read('UL01.loom' cache=True)
ldata_UL02=scv.read('UL02.loom', cache=True)
ldata_UL03=scv.read('UL03.loom', cache=True)
ldata_UL04=scv.read('UL04.loom', cache=True)
ldata_VR01=scv.read('VR01.loom', cache=True)
ldata_VR02=scv.read('VR02.loom', cache=True)
ldata_VR03=scv.read('VR03.loom', cache=True)
ldata_VR04=scv.read('VR04.loom', cache=True)
ldata_DQ01=scv.read('DQ01t.loom', cache=True)
ldata_DQ02=scv.read('DQ02t.loom', cache=True)
ldata_DQ03=scv.read('DQ03t.loom', cache=True)
ldata_DQ04=scv.read('DQ04t.loom', cache=True)
ldata_DQ06=scv.read('DQ06t.loom', cache=True)
ldata_DQ07=scv.read('DQ07t.loom', cache=True)
ldata_DQ08=scv.read('DQ08t.loom', cache=True)

# check the pattern of 'barcode'
adata_XGE202122.obs['barcode']

# change the pattern of 'barcode' in ldata
barcodes=[bc.split(':')[1] for bc in ldata_UL01.obs.index.tolist()]
barcodes=[bc[0:len(bc)-1]+'-1_1_1' for bc in barcodes]
ldata_DQ01.obs.index=barcodes
# repeat for all ldata

# concatenate the looms
 ldata_XGE202122=ldata_UL01.concatenate([ldata_UL02, ldata_UL04, ldata_VR01, ldata_VR02, ldata_VR03, ldata_VR04, ldata_DQ01, ldata_DQ02, ldata_DQ03, ldata_DQ04, ldata_DQ06, ldata_DQ07, ldata_DQ08])

# rewrite obs index of merged loom
ldata_cell=ldata_XGE22.obs.index.tolist()
ldata_cell=[bc[0:len(bc)-2]for bc in ldata_cell]

# or cell2 = [c.rsplit('-', 1)[0] for c in ldata_cell] 
# if the serial number is not all two digit

ldata_XGE202122.obs.index=ldata_cell2
adata_cell=adata_XGE202122.obs.index.tolist()
adata_gene=adata_XGE202122.var.index.tolist()

# replace specific string if needed
adata_gene=[bc.replace('-XLOC-', '_XLOC_') for bc in adata_XGE202122_PB.var.index.tolist()]

# subset ldata or adata for merge (here subset ldata for adata)
adata_cell=adata_XGE22.obs.index.tolist()
adata_gene=adata_XGE22.var.index.tolist()


# subset ldata for merge
ldata_XGE202122_sub=ldata_XGE202122[ldata_XGE202122.obs.index.isin(adata_cell)]
ldata_XGE202122_sub=ldata_XGE202122_sub[:,ldata_XGE202122_sub.var.index.isin(adata_gene)]
adata_XGE202122_sub=adata_XGE202122[:, adata_XGE202122.var.index.isin(ldata_gene)]


# merge matrices into the original adata object
adata_XGE202122_merge=scv.utils.merge(adata_XGE202122_sub, ldata_XGE202122_sub)
adata_XGE202122_merge.write('adata_XGE202122_merged.h5ad')

# subset by domains
adata_XGE22_subcomined=adata_XGE22_combined[adata_XGE22_combined.obs['domains'].isin(['SAM','P1_P2','P3','P4','P5'])]

#computing velocity
scv.pl.proportions(adata_XGE202122_merge, groupby='domains')
scv.pp.filter_and_normalize(adata_XGE202122_merge)
scv.pp.moments(adata_XGE202122_merge)
scv.tl.velocity(adata_XGE202122_merge, mode='stochastic')
scv.tl.velocity_graph(adata_XGE202122_merge)
scv.pl.velocity_embedding(adata_XGE202122_merge, basis='umap', frameon=False, save='velocity_embedding.png')
scv.pl.velocity_embedding_grid(adata_XGE202122_merge, basis='umap', color='domains', save='velocity_embedding_grid.png', title='', scale=0.25)
scv.pl.velocity_embedding_stream(adata_XGE202122_merge, basis='umap', color=['domains'], save='velocity_embedding_stream.png', fitle='')
scv.tl.recover_dynamics(adata_XGE202122_merge, mode='dynamical')
scv.pl.velocity_embedding_stream(adata_XGE202122_merge, basis='umap', color=['domains'], save='velocity_embedding_stream_dynamic.png', title='')
scv.pl.velocity_embedding_stream(adata_XGE202122_merge, basis='pca', color=['domains'], save='velocity_embedding_stream_pca_dynamic.png', title='')
adata_XGE202122_sub.write('adata_XGE202122_sub.h5ad')
ldata_XGE202122_sub.write('ldata__XGE202122_sub.h5ad')

scv.pl.velocity_embedding_stream(adata_XGE202122_merge, basis='umap', color=['harmony_clusters'], save='velocity_embedding_stream_dynamic_clusters.png', title='', density=2, max_length=1)
scv.pl.velocity_embedding_stream(adata_XGE202122_merge, basis='umap', color=['domains','harmony_clusters'], save='velocity_embedding_stream_dynamic_domains_clusters.png', title='', density=2, max_length=1)
