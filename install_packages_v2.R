.libPaths( c( .libPaths(), "/home/ubuntu/workspace/lib/R/site-library") )
options(install.packages.compile.from.source = "always")
install.packages("config",
                 dependencies=TRUE, 
                 repos='http://cran.rstudio.com/',
                 lib="/home/ubuntu/workspace/lib/R/site-library")
remotes::install_github("satijalab/seurat-wrappers", lib="/home/ubuntu/workspace/lib/R/site-library")