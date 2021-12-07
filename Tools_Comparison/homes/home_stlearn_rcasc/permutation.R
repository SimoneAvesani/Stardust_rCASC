library("Seurat")
library("argparser")
library(dplyr)
library(Matrix)
library(stringi)

p <- arg_parser("permutation")
p <- add_argument(p, "percent", help="Percentage of cell removed for bootstrap algorithm ")
p <- add_argument(p, "matrixName", help="matrix count name")
p <- add_argument(p, "lResolution", help="double for the resolution of Louvain algorithm")
p <- add_argument(p, "pcaDimensions", help="PCA dimensions")
p <- add_argument(p, "index", help="Clulstering method: SIMLR tSne Griph")

argv <- parse_args(p)
cat(system("pwd"))
matrixName=argv$matrixName
lResolution=argv$lResolution
percent=as.numeric(argv$percent)
pcaDimensions=as.numeric(argv$pcaDimensions)
index=argv$index

source("./../../../home/functions.R")
countMatrix=read.table(paste("/scratch/",matrixName,sep=""),sep="\t",header=TRUE,row.name=1)
countMatrix <- countMatrix[,sort(colnames(countMatrix))]
suffix = stri_rand_strings(length=5,n=1)
killedCellFile = paste("/scratch/killedCell",suffix,".txt",sep="")
killedCell <- sample(ncol(countMatrix),(ncol(countMatrix)*percent/100))
killedCellNames = sort(killedCell)
killedCellNames = colnames(countMatrix)[killedCellNames]
killedCellNames = gsub("[.]","-",killedCellNames)
write.table(killedCellNames,killedCellFile,quote = FALSE,row.names = FALSE,col.names = FALSE)

stLearnout = paste("/scratch/stLearnout",suffix,sep="")
system(paste("mkdir ",stLearnout,sep=""))
inputDataset = "/scratch/filtered_feature_bc_matrix"
lResolution = lResolution
runstLearn = paste("python /home/mainScript.py ",inputDataset," ",stLearnout," ",pcaDimensions," ",lResolution," yes ",killedCellFile," ",sep = "")
system(runstLearn)
mainVector = read.table(paste(stLearnout,"/Clusters.txt",sep=""),header=TRUE,row.names=1)
mainVector = mainVector + 1

write.table(mainVector,paste("./Permutation/clusterB_",index,".","txt",sep=""),sep="\t")
write.table(killedCell,paste("./Permutation/killC_",index,".","txt",sep=""),sep="\t")
rm(list=setdiff(ls(),"index"))
dir.create("./memory")
system(paste("cat /proc/meminfo >  ./memory/",index,".txt",sep=""))