setwd("/home")

source("functions.R")
library(Seurat)
library("argparser")
library(dplyr)
library("vioplot")
library("Publish")
 
p <- arg_parser("permutation")
p <- add_argument(p, "matrixName", help="matrix count name")
p <- add_argument(p, "tissuePositionsName", help="matrix count name")
p <- add_argument(p, "imageName", help="matrix count name")
p <- add_argument(p, "use_histology", help="bool, use histological info or not")
p <- add_argument(p, "lResolution", help="double for the resolution of Louvain algorithm")
p <- add_argument(p, "pcaDimensions", help="PCA dimension for seurat first number")
p <- add_argument(p, "p", help="real, percentage of total expression contributed by neighborhoods")
p <- add_argument(p, "nPerm", help="Permutation number for bootstrap algorithm ")
p <- add_argument(p, "permAtTime", help="Number of permutation in parallel")
p <- add_argument(p, "percent", help="Percentage of cell removed for bootstrap algorithm ")
p <- add_argument(p, "seed", help="Seed necessary for the reproducibility")


argv <- parse_args(p)


options(bitmapType='cairo')
Sys.setenv("DISPLAY"=":0.0")
matrixName=argv$matrixName
tissuePositionsName=argv$tissuePositionsName
imageName=argv$imageName
use_histology=argv$use_histology
lResolution=argv$lResolution
pcaDimensions=as.numeric(argv$pcaDimensions)
p=argv$p
nPerm=as.numeric(argv$nPerm)
permAtTime=as.numeric(argv$permAtTime)
percent=as.numeric(argv$percent)
seed=as.numeric(argv$seed)
set.seed(seed)
matrixNameBis = strsplit(matrixName,".",fixed = TRUE)[[1]][1]
dir.create(paste("./../scratch/",matrixNameBis,sep=""))
 


setwd(paste("./../scratch/",matrixNameBis,"/",sep=""))
nCluster=clustering(matrixName,tissuePositionsName,imageName,use_histology,lResolution,pcaDimensions,p,nPerm,permAtTime,percent,nCluster=0)


setwd("./../../../home")
setwd(paste("./../scratch/",matrixNameBis,"/",sep=""))
silhouettePlot(matrixNameBis,nCluster,"txt","\t")

  
#dir.create("./../../data/Results")

#system("cp -r ./../* ./../../data/Results")
setwd("./../..")
#system("rm -r ./scratch/*")





system("chmod -R 777 ./scratch") 
system("chmod -R 777 ./data")
