library("Seurat")
library("argparser")
library(dplyr)
library(Matrix)

p <- arg_parser("permutation")
p <- add_argument(p, "percent", help="matrix count name")
p <- add_argument(p, "matrixName", help="Percentage of cell removed for bootstrap algorithm ")
p <- add_argument(p, "tissuePositionFile", help="file with spot coordinates")
p <- add_argument(p, "profileDistance", help="parameter of minkowski distance bw transcriptional profiles")
p <- add_argument(p, "spotDistance", help="parameter of minkowski distance bw spot position")
p <- add_argument(p, "spotDistanceTransformationWeight", help="transformation of spot distance measure")
p <- add_argument(p, "res", help="resolution for louvain algorithm")
p <- add_argument(p, "format", help="matrix format like csv, txt...")
p <- add_argument(p, "separator", help="matrix separator ")
p <- add_argument(p, "logTen", help="1 or 0 if is matrix is already in log10 or if is not")
p <- add_argument(p, "pcaDimensions", help="PCA dimensions")
p <- add_argument(p, "sparse", help="PCA dimensions")
p <- add_argument(p, "index", help="Clulstering method: SIMLR tSne Griph")

argv <- parse_args(p)
cat(system("pwd"))
matrixName=argv$matrixName
tissuePositionFile=argv$tissuePositionFile
profileDistance=argv$profileDistance
spotDistance=argv$spotDistance
spotDistanceTransformationWeight=as.double(argv$spotDistanceTransformationWeight)
res=argv$res
percent=as.numeric(argv$percent)
format=argv$format
separator=argv$separator
logTen=as.numeric(argv$logTen)
pcaDimensions=as.numeric(argv$pcaDimensions)
index=argv$index
sparse=argv$sparse

source("./../../../home/functions.R")
if(separator=="tab"){separator="\t"}

if(sparse=="FALSE"){
    countMatrix <- as.matrix(read.table(paste("./../../",matrixName,".",format,sep=""),sep=separator,header=TRUE,row.names=1))
    if(logTen==1){countMatrix=10^(countMatrix)}
    countMatrix <- countMatrix[,sort(colnames(countMatrix))]
    killedCell <- sample(ncol(countMatrix),(ncol(countMatrix)*percent/100))
    countMatrix <- countMatrix[,-killedCell]
}else{
    countMatrix <- Read10X(data.dir = "./..")
    if(logTen==1){
        stop("Sparse Matrix in Seurat has to be raw count")
    }
    countMatrix <- countMatrix[,sort(colnames(countMatrix))]
    killedCell <- sample(countMatrix@Dim[2],(countMatrix@Dim[2]*percent/100))
    countMatrix <- as.matrix(countMatrix)
    countMatrix <- countMatrix[,-killedCell]
    countMatrix <- Matrix(countMatrix, sparse = TRUE) 
}
res=as.double(res)
tissuePosition <- as.matrix(read.table(paste("./../../",tissuePositionFile,sep=""),header=TRUE,sep="\t",row.names=1))
d <- dim(tissuePosition)[2]
tissuePosition <- tissuePosition[,(d-1):d]
tissuePosition <- tissuePosition[sort(rownames(tissuePosition)),]
tissuePosition <- tissuePosition[-killedCell,]
pbmc <- CreateSeuratObject(countMatrix)
pbmc <- SCTransform(pbmc, assay = "RNA", verbose = FALSE)
pbmc <- RunPCA(pbmc, assay = "SCT", verbose = FALSE)
m = pbmc@reductions[["pca"]]@cell.embeddings[,1:pcaDimensions]
distPCA = dist(m,method="minkowski",p=as.numeric(profileDistance))  
distCoord <- dist(tissuePosition,method="minkowski",p=as.numeric(spotDistance))
distCoord <- distCoord*((max(distPCA)*as.double(spotDistanceTransformationWeight))/(max(distCoord)))
finalDistance <- as.matrix(distPCA + distCoord)
neighbors <- FindNeighbors(finalDistance)
neighbors <- list(neighbors_nn=neighbors$nn,neighbors_snn=neighbors$snn)
pbmc@graphs <- neighbors
pbmc <- FindClusters(pbmc, dims.use = 1:pcaDimensions, resolution=res,verbose = FALSE, graph.name = "neighbors_snn")
mainVector <- pbmc@active.ident
clustering.output <- mainVector

write.table(mainVector,paste("./Permutation/clusterB_",index,".",format,sep=""),sep=separator)
write.table(killedCell,paste("./Permutation/killC_",index,".",format,sep=""),sep=separator)
rm(list=setdiff(ls(),"index"))
dir.create("./memory")
system(paste("cat /proc/meminfo >  ./memory/",index,".txt",sep=""))