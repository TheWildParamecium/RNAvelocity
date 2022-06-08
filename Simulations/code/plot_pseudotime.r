library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
file_dir <- args[1]
results_dir <- gsub("datasets", "results", file_dir)

print(results_dir)
velo_cors = read.table(paste(file_dir, "/velos_pseudotime_cors.csv", sep=""), header = TRUE, sep=",")
n_velometrics <- nrow(velo_cors)

ti_cors = read.table(paste(file_dir, "/ti_pseudotime_cors.csv", sep=""), header = TRUE, sep=",")
n_TImetrics <- nrow(ti_cors)


total_cors = data.frame(velocyto = velo_cors[,"velocyto"], 
                        scvelo = velo_cors[,"scvelo"],
                        slingshot = ti_cors[,"slingshot"], 
                        paga = ti_cors[,"paga"])

ncols = ncol(total_cors)
nrows = nrow(total_cors)


reformated = data.frame(
    correlation = c(velocyto = velo_cors[,"velocyto"], scvelo = velo_cors[,"scvelo"],
                    slingshot = ti_cors[,"slingshot"], paga = ti_cors[,"paga"]),
    method = c(rep("velocyto",n_velometrics), rep("scvelo",n_velometrics), 
                rep("slingshot",n_TImetrics), rep("paga",n_TImetrics))
)

#Boxplot of genwise correlations
title = "Correlation between simulated and predicted pseudotime"

png(filename = paste0(results_dir,"/pseudotime_correlations.png", sep=""), width = 1000, height = 1000)
ggplot(reformated, aes(x=correlation, y=method, fill=method)) +
    geom_boxplot(alpha=0.3) + 
    xlab("Correlation") + ylab("Method") +
    ggtitle(paste(strwrap(title, width = 70), collapse = "\n")) + 
    theme(text = element_text(size = 20), plot.title = element_text(hjust =0.5)) +
    coord_flip()
dev.off()


#Plotting the pseudotime correlation
