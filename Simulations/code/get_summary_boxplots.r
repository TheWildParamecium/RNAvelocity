library(tidyverse)
library(dplyr)
library(ggplot2)
library(stringr)

debugging = TRUE

#Setting needed directories
args = commandArgs(trailingOnly=TRUE)
file_dir <- args[1]

results_dir<- gsub("datasets", "results", file_dir)

#Reading gen and cell-wise velocity correlations, length and confidences
total_gene_correlations <- data.frame(values=numeric(), replicates=character(), tool=character(), stringsAsFactors=FALSE)
total_cell_correlations <- data.frame(values=numeric(), replicates=character(), tool=character(), stringsAsFactors=FALSE)
total_lengths <- data.frame(values=numeric(), replicates=character(), tool=character(), state=character(), stringsAsFactors=FALSE)
total_confidences <- data.frame(values=numeric(), replicates=character(), tool=character(), state=character(), stringsAsFactors=FALSE)


for (muestra in paste("dataset", 1:10, sep="")){
    #Reading the required files
    scvelo_cells_cors = read.csv2(paste(file_dir, muestra, "scvelo_cell_cors.csv", sep=""), sep=",", header = FALSE)
    velocyto_cells_cors = read.csv2(paste(file_dir, muestra,"velocyto_cell_cors.csv", sep=""), sep=",", header = FALSE)
    
    scvelo_gene_cors = read.csv2(paste(file_dir, muestra,"scvelo_gene_cors.csv", sep=""), sep=",", header = FALSE)
    velocyto_gene_cors = read.csv2(paste(file_dir, muestra,"velocyto_gene_cors.csv", sep=""), sep=",", header = FALSE)

    scvelo_lengths = read.csv2(paste(file_dir, muestra,"lengths_scvelo.csv", sep=""), sep=",", header = FALSE)
    velocyto_lengths = read.csv2(paste(file_dir, muestra,"lengths_velocyto.csv", sep=""), sep=",", header = FALSE)

    scvelo_confidences = read.csv2(paste(file_dir, muestra,"confidences_scvelo.csv", sep=""), sep=",", header = FALSE)
    velocyto_confidences = read.csv2(paste(file_dir, muestra,"confidences_velocyto.csv", sep=""), sep=",", header = FALSE)

    metadata = read.csv2(paste(file_dir, muestra,"_metadata.csv", sep=""), sep=",", header = TRUE)

    #Creating empty dataframes with required columns for saving the metrics of each sample
    sample_gene_correlations <- data.frame(values=numeric(length(scvelo_gene_cors[,1])*2), replicates=character(length(scvelo_gene_cors[,1])*2), tool=character(length(scvelo_gene_cors[,1])*2), stringsAsFactors=FALSE)
    sample_cell_correlations <- data.frame(values=numeric(length(scvelo_cells_cors[,1])*2), replicates=character(length(scvelo_cells_cors[,1])*2), tool=character(length(scvelo_cells_cors[,1])*2), stringsAsFactors=FALSE)
    sample_lengths <- data.frame(values=numeric(2*length(scvelo_lengths[,1])), replicates=character(2*length(scvelo_lengths[,1])), tool=character(2*length(scvelo_lengths[,1])), state=character(2*length(scvelo_lengths[,1])), stringsAsFactors=FALSE)
    sample_confidences <- data.frame(values=numeric(2*length(scvelo_confidences[,1])), replicates=character(2*length(scvelo_confidences[,1])), tool=character(2*length(scvelo_confidences[,1])), state=character(2*length(scvelo_confidences[,1])), stringsAsFactors=FALSE)

    #Saving required metrics inside datasets
    sample_gene_correlations["values"] = as.numeric(c(scvelo_gene_cors[,1], velocyto_gene_cors[,1]))
    sample_gene_correlations["replicates"] = rep(muestra, (2*length(scvelo_gene_cors[,1])))
    sample_gene_correlations["tool"] = c(rep("scvelo",length(scvelo_gene_cors[,1])), rep("velocyto", length(velocyto_gene_cors[,1])))

    sample_cell_correlations["values"] = as.numeric(c(scvelo_cells_cors[,1], velocyto_cells_cors[,1]))
    sample_cell_correlations["replicates"] = rep(muestra, 2*length(scvelo_cells_cors[,1]))
    sample_cell_correlations["tool"] = c(rep("scvelo",length(scvelo_cells_cors[,1])), rep("velocyto", length(velocyto_cells_cors[,1])))

    sample_lengths["values"] = as.numeric(c(scvelo_lengths[,1], velocyto_lengths[,1]))
    sample_lengths["replicates"] = rep(muestra, 2*length(scvelo_lengths[,1]))
    sample_lengths["tool"] = c(rep("scvelo",length(scvelo_lengths[,1])), rep("velocyto", length(velocyto_lengths[,1])))
    sample_lengths["state"] = rep(ifelse(metadata$simulation_state == "non-dynamical", "non-dynamic", "dynamic"), 2)

    sample_confidences["values"] = as.numeric(c(scvelo_confidences[,1], velocyto_confidences[,1]))
    sample_confidences["replicates"] = rep(muestra, 2*length(scvelo_confidences[,1]))
    sample_confidences["tool"] = c(rep("scvelo",length(scvelo_confidences[,1])), rep("velocyto", length(velocyto_confidences[,1])))
    sample_confidences["state"] = rep(ifelse(metadata$simulation_state == "non-dynamical", "non-dynamic", "dynamic"), 2)

    #Saving sample values inside global dataframes
    total_gene_correlations = rbind(total_gene_correlations, sample_gene_correlations)
    total_cell_correlations = rbind(total_cell_correlations, sample_cell_correlations)
    total_lengths = rbind(total_lengths, sample_lengths)
    total_confidences = rbind(total_confidences, sample_confidences)
}

#Transforming possible NAs values to 0
total_gene_correlations[is.na(total_gene_correlations)] <- 0
total_cell_correlations[is.na(total_cell_correlations)] <- 0
total_lengths[is.na(total_lengths)] <- 0
total_confidences[is.na(total_confidences)] <- 0

#Converting replicates string column to a factor one
total_gene_correlations$replicates <- factor(total_gene_correlations$replicates, levels=paste("dataset", 1:10, sep=""))
total_cell_correlations$replicates <- factor(total_cell_correlations$replicates, levels=paste("dataset", 1:10, sep=""))
total_lengths$replicates <- factor(total_lengths$replicates, levels=paste("dataset", 1:10, sep=""))
total_confidences$replicates <- factor(total_confidences$replicates, levels=paste("dataset", 1:10, sep=""))

#Renaming the levels to sample names
levels(total_gene_correlations$replicates) <- paste("sample", 1:10, sep="")
levels(total_cell_correlations$replicates) <- paste("sample", 1:10, sep="")
levels(total_lengths$replicates) <- paste("sample", 1:10, sep="")
levels(total_confidences$replicates) <- paste("sample", 1:10, sep="")


if(debugging){
    print("Printing total_gene_correlations boxplot")
}

#Boxplot of genwise correlations
png(filename = paste0(results_dir,"genwise_correlations.png", sep=""), width = 2000, height = 1000)
ggplot(total_gene_correlations, aes(x=replicates, y=values, fill=tool)) +
    geom_boxplot(alpha=0.3) +  
    xlab("Sample") + ylab("Correlation") +
    theme(text = element_text(size = 50), axis.text.x = element_text(size = 30, angle = 15))
dev.off()

if(debugging){
    print("Printing total_cell_correlations boxplot")
}

#Boxplot of cellwise correlations
png(filename = paste0(results_dir,"cellwise_correlations.png", sep=""), width = 2000, height = 1000)
ggplot(total_cell_correlations, aes(x=replicates, y=values, fill=tool)) +
    geom_boxplot(alpha=0.3) +  
    xlab("Sample") + ylab("Correlation") +
    theme(text = element_text(size = 50), axis.text.x = element_text(size = 30, angle = 15))
dev.off()

if(debugging){
    print("Printing lengths boxplot")
}

#Boxplot of lengths
png(filename = paste0(results_dir,"lengths.png", sep=""), width = 2000, height = 1000)
ggplot(total_lengths, aes(x=replicates, y=values, fill=state)) +
    geom_boxplot(alpha=0.3) + 
    facet_wrap(~tool, nrow = 2, ncol=1) +
    xlab("Sample") + ylab("Scaled length") +
    theme(text = element_text(size = 50), axis.text.x = element_text(size = 30, angle = 15))
dev.off()

if(debugging){
    print("Printing confidences boxplot")
}

#Boxplot of confidences
png(filename = paste0(results_dir,"confidences.png", sep=""), width = 2000, height = 1000)
ggplot(total_confidences, aes(x=replicates, y=values, fill=state)) +
    geom_boxplot(alpha=0.3) +  
    facet_wrap(~tool, nrow = 2, ncol=1) +
    xlab("Sample") + ylab("Confidence") + 
    theme(text = element_text(size = 50), axis.text.x = element_text(size = 30, angle = 15))
dev.off()

if(debugging){
    print("FINISHED. RESULTS AT:")
    print(results_dir)
}