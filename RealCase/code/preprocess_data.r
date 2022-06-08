#Loading needed libraries
library(Seurat)             # General workflow
library(SeuratDisk)         # Loading datasets
library(SeuratWrappers)     # Wrappers for Seurat (libraries developed by other adapted to work with Seurat)
library(clusterProfiler)    # Functional significance analysis (by over-representation of GO terms)
library(celldex)            # bulk and single cell RNA-seq datasets for performing cell annotation
library(scRNAseq)           # Single-cell RNA-seq datasets
library(SingleR)            # Perform cell annotation comparing your query dataset against a (usually manually curated) reference dataset
library(scuttle)            # Processing single-cell object to make SingleR work
library(slingshot)          # Perform pseudotrajectory calculation
library(monocle3)           # Perform pseudotrajectory calculation
library(Nebulosa)           # Plots density graphs of two genes in UMAP # (useful for manual annotation of cell types)
library(dplyr)              # Suitable data rearranging
library(ggplot2)            # Improved plotting features
library(patchwork)
library(magrittr)
library(grDevices)
library(RColorBrewer)
library(velocyto.R)
library(SeuratWrappers)
library(BiocFileCache)
library(stringr)

library(future)             # Parallelization
library(parallel)           # Parallelization
library(BiocParallel)       # Parallelization
cores <- detectCores(all.tests = FALSE, logical = TRUE) - 1
plan("multicore", workers = cores)
options(future.globals.maxSize = 8000 * 1024^2)
set.seed(1)

#Setting directories
args = commandArgs(trailingOnly=TRUE)
file_dir <- args[1]
dataset_dir <- str_split(file_dir, "/", simplify = TRUE)
dataset_dir <- paste(dataset_dir[1, 1:(length(dataset_dir)-1)], collapse="/")

splitted_results_dir_full <- gsub("datasets", "results", file_dir)
splitted_results_dir_full <- str_split(splitted_results_dir_full, "/", simplify = TRUE)

results_dir <- paste(splitted_results_dir_full[1, 1:(length(splitted_results_dir_full)-1)], collapse="/")
filename <- splitted_results_dir_full[1,length(splitted_results_dir_full)]
filename <- str_split(filename, ".loom", simplify = TRUE)[1,1]

results_dir_full <- paste(results_dir, "/",filename,"/", sep="")

system(paste("mkdir -p", results_dir_full))

#Loading and preprocessing the data
full_path <- paste(dataset_dir, "/",filename, ".loom", sep="")
ldat <- ReadVelocity(file = full_path)
seuratObject <- as.Seurat(x = ldat)
seuratObject[["RNA"]] <- seuratObject[["spliced"]]
seuratObject[["percent.mt"]] <- PercentageFeatureSet(seuratObject, pattern = "^MT-")
seuratObject <- subset(seuratObject, subset = nFeature_RNA < 2500 & percent.mt < 5)
write.table(c(paste("ncells:", ncol(seuratObject)), paste("ntranscripts:", nrow(seuratObject))), file=paste(results_dir_full,"nCellsAfterSeuratProccessing.csv",sep=""))

seuratObject <- NormalizeData(seuratObject, normalization.method = "LogNormalize", scale.factor = 10000)
seuratObject <- FindVariableFeatures(seuratObject, selection.method = "vst", nfeatures = 2000)
seuratObject <- ScaleData(seuratObject)
seuratObject <- RunPCA(seuratObject)
dimensions <- 20
seuratObject <- FindNeighbors(object = seuratObject, dims=1:dimensions)
findClustersPromise <- future( FindClusters(object = seuratObject, resolution = 3), seed=TRUE)
seuratObject <- value(findClustersPromise)
seuratObject <- RunUMAP(object = seuratObject, dims=1:dimensions, n.neighbors = 80L, min.dist = 0.5, spread=1)

#Annotating the data with cell types with SingleR and Human Cell Atlas
ref <- celldex::HumanPrimaryCellAtlasData()
raw.seuratObject <- as.SingleCellExperiment(seuratObject)
pred.seuratObject <- SingleR(test = raw.seuratObject, ref = ref,
                     labels = ref$label.main, clusters = Idents(seuratObject),
                     assay.type.test = 1, BPPARAM = MulticoreParam(8))

#Adding predicted cell annotation to each cell
cluster_names <- pred.seuratObject$labels[Idents(seuratObject)]
names(cluster_names) <- colnames(seuratObject)
seuratObject <- AddMetaData( object = seuratObject, metadata = cluster_names, col.name = "cluster_names")

png(file = paste(results_dir_full, filename, "_clusters.png", sep=""), width = 1800, height = 900)
DimPlot(seuratObject, label=TRUE, reduction = "umap", pt.size = 0.5) + NoLegend() + DimPlot(seuratObject, reduction = "umap", group.by="cluster_names", pt.size = 0.5)
dev.off()

#Saving the data
saveRDS(seuratObject, file = paste(dataset_dir,"/",filename,".rds", sep=""))
DefaultAssay(seuratObject) <- "RNA"

SaveH5Seurat(seuratObject, filename = paste(dataset_dir,"/", filename, ".h5Seurat", sep=""), overwrite = TRUE)
Convert(paste(dataset_dir,"/", filename, ".h5Seurat", sep=""), dest = "h5ad", overwrite = TRUE)
