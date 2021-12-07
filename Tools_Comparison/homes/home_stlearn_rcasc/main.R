setwd("/home")

source("functions.R")
library(Seurat)
library("argparser")
library(dplyr)
library("vioplot")
library("Publish")
 
p <- arg_parser("permutation")
p <- add_argument(p, "matrixName", help="matrix count name")
p <- add_argument(p, "lResolution", help="double for the resolution of Louvain algorithm")
p <- add_argument(p, "nPerm", help="Permutation number for bootstrap algorithm ")
p <- add_argument(p, "permAtTime", help="Number of permutation in parallel")
p <- add_argument(p, "percent", help="Percentage of cell removed for bootstrap algorithm ")
p <- add_argument(p, "pcaDimensions", help="PCA dimension for seurat first number")
p <- add_argument(p, "seed", help="Seed necessary for the reproducibility")

argv <- parse_args(p)

options(bitmapType='cairo')
Sys.setenv("DISPLAY"=":0.0")
matrixName=argv$matrixName
lResolution=argv$lResolution
nPerm=as.numeric(argv$nPerm)
permAtTime=as.numeric(argv$permAtTime)
percent=as.numeric(argv$percent)
seed=as.numeric(argv$seed)
pcaDimensions=as.numeric(argv$pcaDimensions)
set.seed(seed)
matrixNameBis = strsplit(matrixName,".",fixed = TRUE)[[1]][1]
dir.create(paste("./../scratch/",matrixNameBis,sep=""))

setwd(paste("./../scratch/",matrixNameBis,"/",sep=""))
nCluster=clustering(matrixName,lResolution,nPerm,permAtTime,percent,nCluster=0,pcaDimensions)

setwd("./../../../home")
setwd(paste("./../scratch/",matrixNameBis,"/",sep=""))
silhouettePlot(matrixNameBis,nCluster,"txt","\t")

setwd("./../..")

system("chmod -R 777 ./scratch") 
system("chmod -R 777 ./data")
