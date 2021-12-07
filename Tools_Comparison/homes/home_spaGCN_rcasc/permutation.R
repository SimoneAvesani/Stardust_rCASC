library("Seurat")
library("argparser")
library(dplyr)
library(Matrix)
library(stringi)

p <- arg_parser("permutation")
p <- add_argument(p, "percent", help="matrix count name")
p <- add_argument(p, "matrixName", help="Percentage of cell removed for bootstrap algorithm ")
p <- add_argument(p, "tissuePositionsName", help="file with spot coordinates")
p <- add_argument(p, "imageName", help="parameter of minkowski distance bw transcriptional profiles")
p <- add_argument(p, "lResolution", help="double for the resolution of Louvain algorithm")
p <- add_argument(p, "pcaDimensions", help="PCA dimension for seurat first number")
p <- add_argument(p, "use_histology", help="parameter of minkowski distance bw spot position")
p <- add_argument(p, "p", help="real, percentage of total expression contributed by neighborhoods")
p <- add_argument(p, "index", help="Clulstering method: SIMLR tSne Griph")

argv <- parse_args(p)
cat(system("pwd"))
matrixName=argv$matrixName
tissuePositionsName=argv$tissuePositionsName
imageName=argv$imageName
lResolution=argv$lResolution
pcaDimensions=as.numeric(argv$pcaDimensions)
use_histology=argv$use_histology
p=argv$p
index=argv$index
percent=as.numeric(argv$percent)

source("./../../../home/functions.R")
countMatrix = Read10X_h5(paste("/scratch/",matrixName,sep=""))
suffix = stri_rand_strings(length=5,n=1)
killedCellFile = paste("/scratch/killedCell",suffix,".txt",sep="")
killedCell <- sample(ncol(countMatrix),(ncol(countMatrix)*percent/100))
killedCellNames = sort(killedCell)
killedCellNames = colnames(countMatrix)[killedCellNames]
write.table(killedCellNames,killedCellFile,quote = FALSE,row.names = FALSE,col.names = FALSE)

spaGCNout = paste("/scratch/spaGCNout",suffix,sep="")
system(paste("mkdir ",spaGCNout,sep=""))
matrixFile = paste("/scratch/",matrixName,sep = "")
tissuePositionsFile = paste("/scratch/",tissuePositionsName,sep = "")
imageFile = paste("/scratch/",imageName,sep = "")
runSpaGNC = paste("python /home/mainScript.py ",spaGCNout," ", use_histology," yes ",
    matrixFile," ",tissuePositionsFile," ",imageFile," ",lResolution," ",pcaDimensions," ",p," ",killedCellFile,sep = "")
system(runSpaGNC)
mainVector = read.table(paste(spaGCNout,"/Clusters.txt",sep=""),header=TRUE,row.names=1)
mainVector <- mainVector + 1

write.table(mainVector,paste("./Permutation/clusterB_",index,".","txt",sep=""),sep="\t")
write.table(killedCell,paste("./Permutation/killC_",index,".","txt",sep=""),sep="\t")
rm(list=setdiff(ls(),"index"))
dir.create("./memory")
system(paste("cat /proc/meminfo >  ./memory/",index,".txt",sep=""))