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

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

#FOR TESTING PURPOSES
debugging = True

#Setting the needed paths

file_path = sys.argv[1]
file_path = file_path.split("/")

file_name = file_path[-1]
file_name = file_name.split(".")[0]

file_dir = file_path[:-1]
file_dir = "/".join(file_dir)
file_dir = file_dir + "/"

results_dir = file_dir.replace("datasets", "results")
results_dir_full = results_dir + file_name + "/"


#Creating the results directory
bashCommand = f"mkdir -p {results_dir_full}"
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()


#Reading the different matrices and instanciating the anndata object
unspliced = pd.read_csv(f"{file_dir}/{file_name}_unspliced.csv")
spliced = pd.read_csv(f"{file_dir}/{file_name}_spliced.csv")
metadata = pd.read_csv(f"{file_dir}/{file_name}_metadata.csv")
real_velos = pd.read_csv(f"{file_dir}/{file_name}_velocities.csv")
real_velos = real_velos.transpose()


adata = anndata.AnnData(
X=spliced.transpose(),
obs=metadata,
layers=dict(
    unspliced=unspliced.transpose(),
    spliced=spliced.transpose()
))

#Preprocessing and clustering the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=20)
sc.tl.leiden(adata)
sc.tl.paga(adata)
sc.pl.paga(adata, plot=False)
sc.tl.umap(adata, init_pos='paga')
sc.tl.tsne(adata, n_pcs=20)

if debugging:
    print("Saving files at:", results_dir_full)

#Preprocessing data for scvelo
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata, n_pcs=20, n_neighbors=20)

if sys.argv[2] == "pseudotimeFalse":
    
    colors = "simulation_state"
    #colors = "clusters" if "01" in results_dir_full else "simulation_state"
    
    ############### FIRST CASE: VELOCYTO + PCA ##################
    scv.tl.velocity(adata, mode="deterministic")
    n_cells = adata.layers["velocity"].shape[0]
    n_genes = adata.layers["velocity"].shape[1]

    # Saving velocities and correlations for later evaluation metrics
    velocities1 = adata.layers["velocity"]
    velocities1[np.isnan(velocities1)] = 0
    np.savetxt(file_dir+file_name+"velocyto_vels.csv", velocities1)

    cell_cors = list(np.zeros(n_cells))
    gen_cors = list(np.zeros(n_genes))

    for i in range(n_genes):
        gen_cors[i] = scv.utils.vcorrcoef(real_velos.iloc[:,i],adata.layers["velocity"][:,i],mode="spearman")
    for i in range(n_cells):
        cell_cors[i] = scv.utils.vcorrcoef(real_velos.iloc[i,:],adata.layers["velocity"][i,:],mode="spearman")
    np.savetxt(file_dir+file_name+"velocyto_gene_cors.csv", gen_cors)
    np.savetxt(file_dir+file_name+"velocyto_cell_cors.csv", cell_cors)


    # Plotting the different embeddings
    scv.tl.velocity_embedding(adata, basis="pca", direct_pca_projection=True)
    scv.pl.velocity_embedding_stream(adata, show=False, basis='pca', color=colors, save=results_dir_full+"velocyto_pca.png", legend_loc="right margin", legend_fontweight="black", legend_fontsize=16, title="")    
    #scv.pl.velocity_embedding_grid(adata, show=False, basis='pca', color=colors, save=results_dir_full+"velocyto_pca_grid.png", legend_loc="right margin")
    #scv.pl.velocity_embedding(adata, show=False, basis='pca', color=colors, save=results_dir_full+"velocyto_pca_normal.png", legend_loc="right margin")

    ############### SECOND CASE: VELOCYTO + UMAP and tSNE ##################
    scv.tl.velocity_graph(adata, n_jobs=7)
    scv.tl.velocity_confidence(adata)


    scv.pl.velocity_embedding_stream(adata, show=False, basis='umap', color=colors, save=results_dir_full+"velocyto_umap2.png", legend_loc="right margin", legend_fontweight="black", legend_fontsize=16, title="")
    #scv.pl.velocity_embedding_grid(adata, show=False, basis='umap', color=colors, save=results_dir_full+"velocyto_umap2_grid.png", legend_loc="right margin")
    #scv.pl.velocity_embedding(adata, show=False, basis='umap', color=colors, save=results_dir_full+"velocyto_umap2_normal.png", legend_loc="right margin")

    scv.pl.velocity_embedding_stream(adata, show=False, basis='tsne', color=colors, save=results_dir_full+"velocyto_tsne.png", legend_loc="right margin", legend_fontweight="black", legend_fontsize=16, title="") 
    #scv.pl.velocity_embedding_grid(adata, show=False, basis='tsne', color=colors, save=results_dir_full+"velocyto_tsne_grid.png", legend_loc="right margin")
    #scv.pl.velocity_embedding(adata, show=False, basis='tsne', color=colors, save=results_dir_full+"velocyto_tsne_normal.png", legend_loc="right margin")


    #Calculating, plotting and saving velocity confidences and lengths
    scv.tl.velocity_confidence(adata)
    scaled = adata.obs["velocity_length"] / np.max(adata.obs["velocity_length"])
    adata.obs["velocity_length"] = scaled

    scv.pl.scatter(adata, show=False, c='velocity_confidence', cmap='coolwarm', perc=[5, 95], save=results_dir_full+"velocyto_umap_confidence.png")
    scv.pl.scatter(adata, show=False, c='velocity_length', cmap='coolwarm', perc=[5, 95], save=results_dir_full+"velocyto_umap_velolength.png")

    np.savetxt(file_dir+file_name+"lengths_velocyto.csv", scaled)
    np.savetxt(file_dir+file_name+"confidences_velocyto.csv", adata.obs["velocity_confidence"])

    ################ THIRD CASE: SCVELO + PCA ##################
    scv.tl.recover_dynamics(adata, n_jobs=7)
    scv.tl.velocity(adata, mode="dynamical")

    # Saving velocities and correlation for evaluation metrics
    velocities2 = adata.layers["velocity"]
    velocities2[np.isnan(velocities2)] = 0
    np.savetxt(file_dir+file_name+"scvelo_vels.csv", velocities2)

    cell_cors = list(np.zeros(n_cells))
    gen_cors = list(np.zeros(n_genes))

    for i in range(n_genes):
        gen_cors[i] = scv.utils.vcorrcoef(real_velos.iloc[:,i],adata.layers["velocity"][:,i],mode="spearman")
    for i in range(n_cells):
        cell_cors[i] = scv.utils.vcorrcoef(real_velos.iloc[i,:],adata.layers["velocity"][i,:],mode="spearman")
    np.savetxt(file_dir+file_name+"scvelo_gene_cors.csv", gen_cors)
    np.savetxt(file_dir+file_name+"scvelo_cell_cors.csv", cell_cors)


    # Plotting the embeddings
    scv.tl.velocity_embedding(adata, basis="pca", direct_pca_projection=True)
    scv.pl.velocity_embedding_stream(adata, show=False, basis='pca', color=colors, save=results_dir_full+"scvelo_pca.png", legend_loc="right margin", legend_fontweight="black", legend_fontsize=16, title="")
    #scv.pl.velocity_embedding_grid(adata, show=False, basis='pca', color=colors, save=results_dir_full+"scvelo_pca_grid.png", legend_loc="right margin")
    #scv.pl.velocity_embedding(adata, show=False, basis='pca', color=colors, save=results_dir_full+"scvelo_pca_normal.png", legend_loc="right margin")

    ################# FOURTH CASE: SCVELO + UMAP and tSNE ##################
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_confidence(adata)

    # Plotting the embeddings
    scv.pl.velocity_embedding_stream(adata, show=False, basis='umap', color=colors, save=results_dir_full+"scvelo_umap2.png", legend_loc="right margin", legend_fontweight="black", legend_fontsize=16, title="")
    #scv.pl.velocity_embedding_grid(adata, show=False, basis='umap', color=colors, save=results_dir_full+"scvelo_umap2_grid.png", legend_loc="right margin")
    #scv.pl.velocity_embedding(adata, show=False, basis='umap', color=colors, save=results_dir_full+"scvelo_umap2_normal.png", legend_loc="right margin")
    
    scv.pl.velocity_embedding_stream(adata, show=False, basis='tsne', color=colors, save=results_dir_full+"scvelo_tsne.png", legend_loc="right margin", legend_fontweight="black", legend_fontsize=16, title="")
    #scv.pl.velocity_embedding_grid(adata, show=False, basis='tsne', color=colors, save=results_dir_full+"scvelo_tsne_grid.png", legend_loc="right margin")
    #scv.pl.velocity_embedding(adata, show=False, basis='tsne', color=colors, save=results_dir_full+"scvelo_tsne_normal.png", legend_loc="right margin")
    

    #Calculating, plotting and saving velocity confidences and lengths
    scv.tl.velocity_confidence(adata)
    scaled = adata.obs["velocity_length"] / np.max(adata.obs["velocity_length"])
    adata.obs["velocity_length"] = scaled

    scv.pl.scatter(adata, show=False, c='velocity_confidence', cmap='coolwarm', perc=[5, 95], save=results_dir_full+"scvelo_umap_confidence.png")
    scv.pl.scatter(adata, show=False, c='velocity_length', cmap='coolwarm', perc=[5, 95], save=results_dir_full+"scvelo_umap_velolength.png")

    np.savetxt(file_dir+file_name+"lengths_scvelo.csv", scaled)
    np.savetxt(file_dir+file_name+"confidences_scvelo.csv", adata.obs["velocity_confidence"])
        

    ##################### LATENT TIME AND HEATMAP PLOTS ##############################
    scv.tl.latent_time(adata)
    scv.pl.scatter(adata, show=False, color='latent_time', color_map='gnuplot', size=80, save=results_dir_full+"latent_time.png")
    
    top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:100]
    #scv.pl.heatmap(adata, show=False, var_names=top_genes, sortby='latent_time', col_color='leiden', n_convolve=100, save=results_dir_full + "heatmap_clusters.png")
    scv.pl.heatmap(adata, show=False, var_names=top_genes, sortby='latent_time', col_color=colors, n_convolve=100, save=results_dir_full + "heatmap_typess.png")


    ##################### PLOTTING SOME PHASE PORTRAIT PLOTS ##########################
    
    ########## GLOBAL "BEST" AND "WORST" LIKELIHOOD RATED DYNAMICAL GENES
    top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index
    #WORST
    scv.pl.scatter(adata, show=False, save=results_dir_full + "phase_worst_all1.png", basis=top_genes[-8:-4], ncols=4, color="simulation_state", add_outline="starting, ending", outline_width=(0.03,0.01), wspace=0.5)
    scv.pl.scatter(adata, show=False, save=results_dir_full + "phase_worst_all2.png", basis=top_genes[-4:], ncols=4, color="simulation_state", add_outline="starting, ending", outline_width=(0.03,0.01), wspace=0.5)
    #BEST
    scv.pl.scatter(adata, show=False, save=results_dir_full + "phase_best_all1.png", basis=top_genes[0:4], ncols=4, color="simulation_state", add_outline="starting, ending", outline_width=(0.03,0.01), wspace=0.5)
    scv.pl.scatter(adata, show=False, save=results_dir_full + "phase_best_all2.png", basis=top_genes[4:8], ncols=4, color="simulation_state", add_outline="starting, ending", outline_width=(0.03,0.01), wspace=0.5)

    ########## BEST PART SPECIFIC LIKELIHOOD RATED DYNAMICAL GENES
    scv.tl.rank_dynamical_genes(adata, groupby="simulation_state")
    df = scv.get_df(adata, "rank_dynamical_genes/names")
    for cluster in pd.unique(adata.obs["simulation_state"]):
        scv.pl.scatter(adata, df[cluster][:5], show=False, save=results_dir_full+"best"+cluster+".png", ylabel="unspliced", color="simulation_state", add_outline="starting, ending", outline_width=(0.03,0.01), wspace=0.5)



##################### THIS IS ONLY NEED FOR CALCULATING PSEUDOTIME METRICS #####################
elif sys.argv[2] == "pseudotimeTrue":
    simulated_time = adata.obs["simulated_time"]

    ############## PSEUDOTIME VELOCYTO ##############
    scv.tl.velocity(adata, mode="deterministic")

    # Saving pseudotimes for evaluation metrics
    scv.tl.velocity_graph(adata, n_jobs=7)
    scv.tl.velocity_pseudotime(adata)
    pseudotime = adata.obs["velocity_pseudotime"]
    velocyto_correlation = np.round_(scipy.stats.spearmanr(pseudotime, simulated_time)[0],4)

    ############### PSEUDOTIME SCVELO ##############
    adata = scv.read(file_dir+file_name+".h5ad")
    scv.tl.recover_dynamics(adata, n_jobs=7)
    scv.tl.velocity(adata, mode="dynamical")

    # Saving pseudotimes for evaluation metrics
    scv.tl.velocity_graph(adata, n_jobs=7)
    scv.tl.latent_time(adata)
    latent_time = adata.obs["latent_time"]
    scvelo_correlation = np.round_(scipy.stats.spearmanr(latent_time, simulated_time)[0], 4)

    #SAVING RESULTS
    with open(results_dir+"velos_pseudotime_cors.csv", mode="a") as f:
        f.write(f"{velocyto_correlation},{scvelo_correlation}\n")



    #Deleting h5ad to save memory"
    bashCommand = f"rm {file_dir+file_name}.h5ad"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
