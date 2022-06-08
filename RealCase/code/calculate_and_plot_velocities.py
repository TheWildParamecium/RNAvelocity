import scvelo as scv
import sys
import scipy
import subprocess
import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

#Setting directories

file_path = sys.argv[1]
file_path = file_path.split("/")

file_name = file_path[-1]
file_name = file_name.split(".")[0]

file_dir = file_path[:-1]
file_dir = "/".join(file_dir)
file_dir = file_dir + "/"

results_dir = file_dir.replace("datasets", "results")
results_dir_full = results_dir + file_name + "/"


adata = scv.read(file_dir+file_name+".h5ad")
scv.pp.remove_duplicate_cells(adata)
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
sc.tl.tsne(adata)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)


######################################## VELOCYTO STREAM PLOTS ########################################
scv.tl.velocity(adata, n_jobs=7, mode="deterministic")
scv.tl.velocity_graph(adata, n_jobs=7)
scv.tl.velocity_confidence(adata)

#VELOCYTO + UMAP and TSNE STREAMPLOT
scv.pl.velocity_embedding_stream(adata, basis='umap', color="cluster_names", show=False, save=results_dir_full+"velocyto_types_umap.png", legend_loc="right margin") 
scv.pl.velocity_embedding_stream(adata, basis='umap', color="seurat_clusters",show=False, save=results_dir_full+"velocyto_clusters_umap.png")

scv.pl.velocity_embedding_stream(adata, basis='tsne', color="cluster_names", show=False, save=results_dir_full+"velocyto_types_tsne.png", legend_loc="right margin") 

#VELOCYTO + PCA STREAMPLOT
scv.tl.velocity_embedding(adata, basis="pca", direct_pca_projection=True)
scv.pl.velocity_embedding_stream(adata, basis='pca', color="cluster_names",show=False, save=results_dir_full + "velocyto_pca.png", legend_loc="right margin")


#CONFIDENCE AND LENGTH PLOTS
scv.tl.velocity_confidence(adata)
scaled = adata.obs["velocity_length"] / np.max(adata.obs["velocity_length"])
adata.obs["velocity_length"] = scaled

scv.pl.scatter(adata, show=False, c='velocity_confidence', cmap='coolwarm', perc=[5, 95], save=results_dir_full+"velocyto_umap_confidence.png")
scv.pl.scatter(adata, show=False, c='velocity_length', cmap='coolwarm', perc=[5, 95], save=results_dir_full+"velocyto_umap_velolength.png")

sns.distplot(adata.obs["velocity_length"])
plt.savefig(results_dir_full+'velocyto_length_dist.png')

sns.displot(adata.obs["velocity_confidence"])
plt.savefig(results_dir_full+'velocyto_confidence_dist.png')

######################################## SCVELO STREAM PLOTS ########################################
scv.tl.recover_dynamics(adata, n_jobs=7)
scv.tl.velocity(adata, n_jobs=7, mode ="dynamical")
scv.tl.velocity_graph(adata, n_jobs=7)
scv.tl.velocity_confidence(adata)
scv.tl.latent_time(adata)


#SCVELO + UMAP STREAMPLOT
scv.pl.velocity_embedding_stream(adata, basis='umap', color="cluster_names", show=False, save=results_dir_full+"scvelo_types.png", legend_loc="right margin") 
scv.pl.velocity_embedding_stream(adata, basis='umap', color="seurat_clusters",show=False, save=results_dir_full+"scvelo_clusters.png")

scv.pl.velocity_embedding_stream(adata, basis='tsne', color="cluster_names", show=False, save=results_dir_full+"velocyto_types_tsne.png", legend_loc="right margin") 

#SCVELO + PCA STREAMPLOT
scv.tl.velocity_embedding(adata, basis="pca", direct_pca_projection=True)
scv.pl.velocity_embedding_stream(adata, basis='pca', color="cluster_names", show=False, save=results_dir_full + "scvelo_pca.png", legend_loc="right margin")

#CONFIDENCE AND LENGTH PLOTS
#CONFIDENCE AND LENGTH PLOTS
scv.tl.velocity_confidence(adata)
scaled = adata.obs["velocity_length"] / np.max(adata.obs["velocity_length"])
adata.obs["velocity_length"] = scaled

scv.pl.scatter(adata, show=False, c='velocity_confidence', cmap='coolwarm', perc=[5, 95], save=results_dir_full+"scvelo_umap_confidence.png")
scv.pl.scatter(adata, show=False, c='velocity_length', cmap='coolwarm', perc=[5, 95], save=results_dir_full+"scvelo_umap_velolength.png")

sns.distplot(adata.obs["velocity_length"])
plt.savefig(results_dir_full+'scvelo_length_dist.png')

sns.displot(adata.obs["velocity_confidence"])
plt.savefig(results_dir_full+'scvelo_confidence_dist.png')

############################################# LATENT TIME PLOT  #############################################
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, show=False, save=results_dir_full+"scvelo_umap_latent_time.png")

top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, show=False,var_names=top_genes, sortby='latent_time', col_color='seurat_clusters', n_convolve=100, save=results_dir_full + "heatmap_clusters.png")
scv.pl.heatmap(adata, show=False,var_names=top_genes, sortby='latent_time', col_color='cluster_names', n_convolve=100, save=results_dir_full + "heatmap_typess.png")


############################################ TOP LIKELIHOOD GENES #############################################

#Top likelihood genes
np.save(file_dir + "global_likelihood_genes.csv", adata.var['fit_likelihood'].sort_values(ascending=False).index[0:50])