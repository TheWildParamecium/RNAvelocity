.libPaths( c( .libPaths(), "/home/ubuntu/workspace/lib/R/site-library") )
system("sudo chmod 666 /var/run/docker.sock")
library(VeloSim)
library(foreach)
library(doParallel)
library(ape)
library(Seurat)             # General workflow
library(SeuratDisk)         # Loading datasets
library(SeuratWrappers)     # Wrappers for Seurat (libraries developed by other adapted to work with Seurat)
library(velocyto.R)
library(scRNAseq)
library(SingleCellExperiment)
library(velociraptor)
library(scuttle)
library(scran)
library(scater)
library(lsa)
library(dyno)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(anndata)
library(SingleCellExperiment)
library(Seurat)
library(dyngen)
library(Matrix)

registerDoParallel(detectCores()-3)

args = commandArgs(trailingOnly=TRUE)

option = args[1]
iterations <- as.numeric(args[2])
capture.rate = as.numeric(args[3])

phyla <- ape::read.tree(text="((A:1,B:1):1.5);")

ngenes <- 500
nevf <- 20
vary <- "all"
Sigma <- 0.4
evf_center <- 1
gene_effect_prob <- 0.3
geffect_mean <- 0
gene_effects_sd <- 1
param_realdata <- "zeisel.imputed"
bimod <- 0
scale_s <- 1

generate_dataset <- function(seed){
    filename <- "dataset"
    n_de_evf_dynamic <- 10
    n_de_evf_static <- 0

    if(args[1] == "1"){
        system("mkdir -p ./datasets/10/")
        dataset_dir <- "./datasets/10/"

        n_dynamic <- 2000
        n_static_1 <- 0
        n_static_2 <- 0
        n_static_3 <- 0 
    }else if(args[1] == "2"){
        system("mkdir -p ./datasets/01/")
        dataset_dir <- "./datasets/01/"
        
        n_dynamic <- 0
        n_static_1 <- 2000
        n_static_2 <- 0
        n_static_3 <- 0
    }else if (args[1] == "3"){
        system("mkdir -p ./datasets/11/")
        dataset_dir <- "./datasets/11/"

        n_dynamic <- 1000
        n_static_1 <- 1000
        n_static_2 <- 0
        n_static_3 <- 0
    }else if (args[1] == "4"){
        system("mkdir -p ./datasets/14/")
        dataset_dir <- "./datasets/14/"

        n_dynamic <- 1000
        n_static_1 <- 600
        n_static_2 <- 300
        n_static_3 <- 100
    }
    print(paste("Generating dataset:", dataset_dir,filename,seed, sep=""))
    print("With parameters:")
    print(list(n_dynamic = n_dynamic, n_static_1 = n_static_1, n_static_2 = n_static_2, n_static_3 = n_static_3))

    randseed <- seed
    if(n_dynamic > 0){
        dynamic <- technicalNoise(SimulateVeloTree(ncells_total=n_dynamic,ngenes=ngenes, evf_center=1,nevf=nevf,
                            phyla=phyla, randseed=randseed, n_de_evf=n_de_evf_dynamic, vary=vary,Sigma=Sigma,geffect_mean=geffect_mean,
                            gene_effects_sd=gene_effects_sd,gene_effect_prob=gene_effect_prob,
                            bimod=bimod,param_realdata=param_realdata,scale_s=scale_s,
                            prop_hge=0.08, mean_hge=5, n_unstable=0, plot = FALSE), capture.rate = capture.rate)    
    }

    randseed <- seed + 10000 + 1
    if(n_static_1 > 0){
    static1 <- technicalNoise(SimulateVeloTree(ncells_total=n_static_1,ngenes=ngenes, evf_center=1,nevf=nevf,
                            phyla=phyla, randseed=randseed, n_de_evf=n_de_evf_static,vary=vary,Sigma=Sigma,geffect_mean=geffect_mean,
                            gene_effects_sd=gene_effects_sd,gene_effect_prob=gene_effect_prob,
                            bimod=bimod,param_realdata=param_realdata,scale_s=scale_s,
                            prop_hge=0.08, mean_hge=5, n_unstable=0, plot = FALSE), capture.rate = capture.rate)
    }

    randseed <- seed + 10000 + 2
    if(n_static_2 > 0){
    static2 <- technicalNoise(SimulateVeloTree(ncells_total=n_static_2,ngenes=ngenes, evf_center=1,nevf=nevf,
                            phyla=phyla, randseed=randseed, n_de_evf=n_de_evf_static,vary=vary,Sigma=Sigma,geffect_mean=geffect_mean,
                            gene_effects_sd=gene_effects_sd,gene_effect_prob=gene_effect_prob,
                            bimod=bimod,param_realdata=param_realdata,scale_s=scale_s,
                            prop_hge=0.08, mean_hge=5, n_unstable=0, plot = FALSE), capture.rate = capture.rate)
    }

    randseed <- seed + 10000 + 3
    if(n_static_3 > 0){
    static3 <- technicalNoise(SimulateVeloTree(ncells_total=n_static_3,ngenes=ngenes, evf_center=1,nevf=nevf,
                            phyla=phyla, randseed=randseed, n_de_evf=n_de_evf_static,vary=vary,Sigma=Sigma,geffect_mean=geffect_mean,
                            gene_effects_sd=gene_effects_sd,gene_effect_prob=gene_effect_prob,
                            bimod=bimod,param_realdata=param_realdata,scale_s=scale_s,
                            prop_hge=0.08, mean_hge=5, n_unstable=0, plot = FALSE), capture.rate = capture.rate)
    }
    
    ################ SAVE DATASET VALUES ##############################
     
    cell_counter <- 1
    pseudotime_static = 50000 #Arbitrary value to tag non dynamic cells
    all_pseudotimes <- c()

    ################# DYNAMIC ##############################
    if(n_dynamic > 0){
    pseudotime_dynamic = dynamic$cell_time
    all_pseudotimes <- c(all_pseudotimes, pseudotime_dynamic)

    dynamic_unspliced <- dynamic$counts_u
    dynamic_spliced <- dynamic$counts_s
    dynamic_velocities <- dynamic$velocity

    colnames(dynamic_unspliced) <- paste("cell", cell_counter:(ncol(dynamic_spliced)), sep="")
    rownames(dynamic_unspliced) <- paste("gen", 1:nrow(dynamic_unspliced), sep="")

    colnames(dynamic_spliced) <- paste("cell", cell_counter:(ncol(dynamic_spliced)), sep="")
    rownames(dynamic_spliced) <- paste("gen", 1:nrow(dynamic_spliced), sep="")

    colnames(dynamic_velocities) <- paste("cell", cell_counter:(ncol(dynamic_velocities)), sep="")
    rownames(dynamic_velocities) <- paste("gen", 1:nrow(dynamic_velocities), sep="")

    cell_counter <- cell_counter + ncol(dynamic_spliced)
    }

    ################ Static 1 ##############################
    if(n_static_1 > 0){
    all_pseudotimes <- c(all_pseudotimes, rep(pseudotime_static, n_static_1-1))

    static1_spliced <- static1$counts_s
    static1_unspliced <- static1$counts_u
    static1_velocities <- static1$velocity

    colnames(static1_spliced) <- paste("cell", cell_counter:(cell_counter+ncol(static1_spliced)-1), sep="")
    rownames(static1_spliced) <- paste("gen", 1:nrow(static1_spliced), sep="")

    colnames(static1_unspliced) <- paste("cell", cell_counter:(cell_counter+ncol(static1_unspliced)-1), sep="")
    rownames(static1_unspliced) <- paste("gen", 1:nrow(static1_unspliced), sep="")

    colnames(static1_velocities) <- paste("cell", cell_counter:(cell_counter+ncol(static1_velocities)-1), sep="")
    rownames(static1_velocities) <- paste("gen", 1:nrow(static1_velocities), sep="")

    cell_counter <- cell_counter + ncol(static1_spliced)
    }

    ################ Static 2 ##############################
    if(n_static_2 > 0){
    all_pseudotimes <- c(all_pseudotimes, rep(pseudotime_static, n_static_2-1))
    
    static2_spliced <- static2$counts_s
    static2_unspliced <- static2$counts_u
    static2_velocities <- static2$velocity

    colnames(static2_spliced) <- paste("cell", cell_counter:(cell_counter+ncol(static2_spliced)-1), sep="")
    rownames(static2_spliced) <- paste("gen", 1:nrow(static2_spliced), sep="")

    colnames(static2_unspliced) <- paste("cell", cell_counter:(cell_counter+ncol(static2_unspliced)-1), sep="")
    rownames(static2_unspliced) <- paste("gen", 1:nrow(static2_unspliced), sep="")

    colnames(static2_velocities) <- paste("cell", cell_counter:(cell_counter+ncol(static2_velocities)-1), sep="")
    rownames(static2_velocities) <- paste("gen", 1:nrow(static2_velocities), sep="")

    cell_counter <- cell_counter + ncol(static2_spliced)
    }
    ################ Static 3 ##############################
    if(n_static_3 > 0){
    all_pseudotimes <- c(all_pseudotimes, rep(pseudotime_static, n_static_3-1))

    static3_spliced <- static3$counts_s
    static3_unspliced <- static3$counts_u
    static3_velocities <- static3$velocity

    colnames(static3_spliced) <- paste("cell", cell_counter:(cell_counter+ncol(static3_spliced)-1), sep="")
    rownames(static3_spliced) <- paste("gen", 1:nrow(static3_spliced), sep="")

    colnames(static3_unspliced) <- paste("cell", cell_counter:(cell_counter+ncol(static3_spliced)-1), sep="")
    rownames(static3_unspliced) <- paste("gen", 1:nrow(static3_unspliced), sep="")

    colnames(static3_velocities) <- paste("cell", cell_counter:(cell_counter+ncol(static3_spliced)-1), sep="")
    rownames(static3_velocities) <- paste("gen", 1:nrow(static3_velocities), sep="")
    }


    ###################### Concatenate all cell clusters ##############################
    all_spliced = c()
    all_unspliced = c()
    all_velocities = c()

    if(n_dynamic > 0){
    all_spliced <- cbind(all_spliced, dynamic_spliced)
    all_unspliced <- cbind(all_unspliced, dynamic_unspliced)
    all_velocities <- cbind(all_velocities, dynamic_velocities)
    }
    if(n_static_1 > 0){
    all_spliced <- cbind(all_spliced, static1_spliced)
    all_unspliced <- cbind(all_unspliced, static1_unspliced)
    all_velocities <- cbind(all_velocities, static1_velocities)
    }
    if(n_static_2 > 0){
    all_spliced <- cbind(all_spliced, static2_spliced)
    all_unspliced <- cbind(all_unspliced, static2_unspliced)
    all_velocities <- cbind(all_velocities, static2_velocities)
    }
    if(n_static_3 > 0){
    all_spliced <- cbind(all_spliced, static3_spliced)
    all_unspliced <- cbind(all_unspliced, static3_unspliced)
    all_velocities <- cbind(all_velocities, static3_velocities)
    }


    #Creating some metadata with Seurat before saving the object
    ldat2 <- list(unspliced=all_unspliced, spliced=all_spliced)
    bm <- as.Seurat(x = ldat2)
    bm[["simulated_time"]] <- all_pseudotimes
    bm[["simulation_state"]] <- rep("dummy", length(all_pseudotimes))
    if(n_dynamic > 0){
    bm[["simulation_state"]] <- ifelse(bm[["simulated_time"]] > quantile(pseudotime_dynamic, 0.9), "ending", "intermediate" )
    bm[["simulation_state"]][pseudotime_dynamic < quantile(pseudotime_dynamic, 0.1),] <- "starting"
    }
    
    bm[["simulation_state"]][bm[["simulated_time"]] > 10000,] <- "non-dynamical"


    #Saving matrices and metadata for scvelo
    write.table(all_spliced, paste(dataset_dir,filename,seed,"_spliced.csv",sep=""), row.names = TRUE, col.names = TRUE, sep = ",")
    write.table(all_unspliced, paste(dataset_dir,filename,seed,"_unspliced.csv",sep=""), row.names = TRUE, col.names = TRUE, sep = ",")
    write.table(all_velocities, paste(dataset_dir,filename,seed,"_velocities.csv",sep=""), row.names = TRUE, col.names = TRUE, sep = ",")
    write.table(bm[[]], paste(dataset_dir,filename,seed,"_metadata.csv",sep=""), row.names = TRUE, col.names = TRUE, sep = ",")
}


foreach (i=1:iterations) %dopar% {
    generate_dataset(i)
}