library("Seurat")
library("argparser")
library(dplyr)
library(Matrix)
library(stringi)
library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)

p <- arg_parser("permutation")
p <- add_argument(p, "percent", help="matrix count name")
p <- add_argument(p, "n_clusters", help="how many clusters BayeSpace should search")
p <- add_argument(p, "pcaDimensions", help="PCA dimension for seurat first number")
p <- add_argument(p, "index", help="Clulstering method: SIMLR tSne Griph")

argv <- parse_args(p)
cat(system("pwd"))
n_clusters=as.integer(argv$n_clusters)
pcaDimensions=as.integer(argv$pcaDimensions)
index=argv$index
percent=as.numeric(argv$percent)

source("./../../../home/functions.R")
sce <- readVisium("/scratch/")
sce = sce[,sort(colnames(sce))]
sce = sce[rowSums(counts(sce)) > 10,]
killedCell <- sample(ncol(sce),(ncol(sce)*percent/100))
sce = sce[,-killedCell]
sce <- spatialPreprocess(sce, platform="ST", n.PCs=pcaDimensions, n.HVGs=2000, log.normalize=TRUE)
sce <- spatialCluster(sce, q=n_clusters, platform="ST", d=pcaDimensions,
                        init.method="mclust", model="t", gamma=2,
                        nrep=10000, burn.in=100,
                        save.chain=TRUE)
suffix = stri_rand_strings(length=5,n=1)
png(paste('/scratch/bayespaceout/bayespace_',suffix,".png",sep=""))
clusterPlot(sce)
dev.off()
mainVector = colData(sce)$spatial.cluster
mainVector = as.numeric(mainVector)
jumping_clusters = sort(unique(mainVector))
for(i in 1:length(jumping_clusters)){
    mainVector[mainVector==jumping_clusters[i]] = i
}

write.table(mainVector,paste("./Permutation/clusterB_",index,".","txt",sep=""),sep="\t")
write.table(killedCell,paste("./Permutation/killC_",index,".","txt",sep=""),sep="\t")
rm(list=setdiff(ls(),"index"))
dir.create("./memory")
system(paste("cat /proc/meminfo >  ./memory/",index,".txt",sep=""))