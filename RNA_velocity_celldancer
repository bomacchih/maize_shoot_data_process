# https://guangyuwanglab2021.github.io/cellDancer_website/data_preprocessing.html
import os
import sys
import glob
import pandas as pd
import math
import matplotlib.pyplot as plt
import celldancer as cd
import celldancer.cdplt as cdplt
from celldancer.cdplt import colormap
import scvelo as scv
adata_XGE202122_merged=scv.read_h5ad('/home/macchihlee/XGE202122_short_2024/2024_short_loom/adata_XGE202122_merged.h5ad')
scv.pp.filter_and_normalize(adata_XGE202122_merged, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata_XGE202122_merged, n_pcs=30, n_neighbors=30)
cdutil.adata_to_df_with_embed(adata_XGE202122_merged,
                              us_para=['Mu','Ms'],
                              cell_type_para='celltype',
                              embed_para='X_umap',
                              save_path='cell_type_u_s.csv')
cell_type_u_s=pd.read_csv('/home/macchihlee/R/XGE202122_short_2024/cell_type_u_s.csv')
cellDancer_df=cd.velocity(cell_type_u_s, permutation_ratio=0.125, n_jobs=8)
ax=plt.subplots(figsize=(10,10))
cdplt.scatter_cell(ax, cellDancer_df, alpha=0.5, s=10, velocity=True, min_mass=15, arrow_grid=(20,20))
ax.axis('off')
cellDancer_df=cd.compute_cell_velocity(cellDancer_df=cellDancer_df, projection_neighbor_choice='gene', expression_scale='power10', projection_neighbor_size=10, speed_up=(100, 100))
loss.df, cellDancer_df=cd.velocity(cell_type_u_s, permutation_ratio=0.1, norm_u_s=False, False, n_jobs=8)
gene_list=list(set(cell_type_u_s.gene_name))
loss_df, cellDancer_df=cd.velocity(cell_type_u_s, gene_list=gene_list, permutation_ratio=0.125, n_jobs=8)
cellDancer_df=cd.compute_cell_velocity(cellDancer_df=cellDancer_df, projection_neighbor_size=100)
fig, ax=plt.subplot(figsize=(10,10))
ax=plt.subplot(figsize=(10,10))
domains_color={'SAM': '#ff4800', 'P1_P2': '#B81136', 'P3': 'green', 'P4': '#03B3B0'; 'P5': '#3361A5'}
im=cdplt.scatter_cell(ax, cellDancer_df, colors=domains_color, alpha=0.5, s=10, velocity=True, legend='on', min_mass=5, arrow_grid=(20,20))
ax.axis('off')
plt.show()
fig, ax=plt.subplots(figsize=(10,10))
im=cdplt.scatter_cell(ax, cellDancer_df, colors=domains_color, alpha=0.5, s=20, velocity=True, legend='on', min_mass=5, arrow_grid=(10,10))
ax.axis('off')
plt.show()

cellDancer_df=cd.compute_cell_velocity(cellDancer_df=cellDancer_df, projection_neighbor_size=100, projection_neighbor_choice='embedding', expression_scale='power10')
fig, ax=plt.subplots(figsize=(10,10))
im=cdplt.scatter_cell(ax, cellDancer_df, colors=domains_color, alpha=0.5, s=20, velocity=True, legend='on', min_mass=5, arrow_grid=(15,15), legend_marker_sieze=10)
ax.axis('off')
plt.show()

# set parameters
dt = 0.05
t_total = {dt:int(10/dt)}
n_repeats = 10
cellDancer_df=cd.pseudo_time(cellDancer_df=cellDancer_df, grid=(30,30), dt=dt, t_total=t_total[dt], n_repeats=n_repeats, speed_up=(100,100), n_paths=3, psrng_seeds_diffusion=[i for i in range(n_repeats)], n_jobs=8)
im=cdplt.scatter_cell(ax, cellDancer_df, colors='pseudotime', alpha=0.5, s=70, velocity=False, colorbar='on')
ax.axis('off')
plt.show()
