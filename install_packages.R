.libPaths( c( .libPaths(), "/home/ubuntu/workspace/lib/R/site-library") )
options(install.packages.compile.from.source = "always")
install.packages(c("BiocManager", "tidyverse", "dyngen","doParallel","devtools","Rcpp",
                    "foreach","expm","gee","dplyr", "igraph","dyngen","remotes",
                    "lsa", "phytools", "rgdal", 'Cairo', "rgeos"),
                 dependencies=TRUE, 
                 repos='http://cran.rstudio.com/',
                 lib="/home/ubuntu/workspace/lib/R/site-library")

library(BiocManager)
library(devtools)
library(remotes)
BiocManager::install(c("multtest", "scater", "scran", "scuttle","velociraptor",
                        "SingleCellExperiment", "scRNAseq"), 
                    lib="/home/ubuntu/workspace/lib/R/site-library")

install.packages(c("mutoss", "qqconf", "metap", "Seurat"),
                 dependencies=TRUE, 
                 repos='http://cran.rstudio.com/',
                 lib="/home/ubuntu/workspace/lib/R/site-library")
                    
devtools::install_github("velocyto-team/velocyto.R",
                    lib="/home/ubuntu/workspace/lib/R/site-library")
devtools::install_github("PeterZZQ/VeloSim",
                    lib="/home/ubuntu/workspace/lib/R/site-library")
remotes::install_github("dynverse/dyno",
                    lib="/home/ubuntu/workspace/lib/R/site-library")