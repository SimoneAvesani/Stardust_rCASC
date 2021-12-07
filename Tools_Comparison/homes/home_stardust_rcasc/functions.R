euc.dist = function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

silhouette=function(nCluster,clustering.output){
 
    dataPlot=cbind(as.numeric(clustering.output[,3]),as.numeric(clustering.output[,4])) 
    nCluster=length(unique(clustering.output[,2]))
    mainVector=as.numeric(clustering.output[,2])
    intraScore=c()   
    extraScore=c()
    neighbor=c()
    silhouetteValue=c()

    for(k in 1:(length(dataPlot)/2))
    {
        a=0
        count=0
        #per ogni altro elemento nel suo cluster
        for(j in 1:(length(dataPlot)/2)){
            if(mainVector[k]==mainVector[j])
                {   
                    if(k != j ){
                        a=a+euc.dist(dataPlot[k,],dataPlot[j,])
                        count=count+1
                    }
                }
        }
        intraScore[k]=a/count 
    }
    extraScoreTemp=c()
    extraCountTemp=c()
    for(k in 1:(length(dataPlot)/2))
    {   
        for(s in 1:nCluster){
            extraScoreTemp[s]=0
            extraCountTemp[s]=0
        }
        for(j in 1:(length(dataPlot)/2))
        {
            if(mainVector[k] != mainVector[j]){
                extraScoreTemp[mainVector[j]]=extraScoreTemp[mainVector[j]]+ euc.dist(dataPlot[k,],dataPlot[j,])
                extraCountTemp[mainVector[j]]=extraCountTemp[mainVector[j]]+1
            }
        
        }
        extraScoreTemp=extraScoreTemp[-mainVector[k]]
        extraCountTemp=extraCountTemp[-mainVector[k]]
        extraScore[k]=min(extraScoreTemp/extraCountTemp)
        minIndex=which.min(extraScoreTemp/extraCountTemp)
        if(minIndex>=mainVector[k]){neighbor[k]=minIndex+1}else{neighbor[k]=minIndex}
    }       
    for(u in 1:length(extraScore)){silhouetteValue[u]=(extraScore[u]-intraScore[u])/max(extraScore[u],intraScore[u])}
    silhouette=matrix(cbind(extraScore,intraScore,mainVector,neighbor,silhouetteValue),nrow=length(extraScore))
    colnames(silhouette) = c("extraScore","intraScore","ClusterBelong","Neighbor","SilhouetteValue") # the first row will be the header
    return(cbind(clustering.output,extraScore,intraScore,neighbor,silhouetteValue))
}

clustering=function(matrixName,tissuePositionFile,profileDistance,spotDistance,
    spotDistanceTransformationWeight,res,nPerm,permAtTime,percent,nCluster,logTen,format,
    separator,pcaDimensions){
    if(separator=="tab"){separator2="\t"}else{separator2=separator}
    if(sparse=="FALSE"){
        countMatrix=read.table(paste("./../",matrixName,".",format,sep=""),sep=separator2,header=TRUE,row.name=1)
        if(logTen==1){countMatrix=10^(countMatrix)}
    }else{
        countMatrix <- Read10X(data.dir = "./")
        if(logTen==1){
            stop("Sparse Matrix in Seurat has to be raw count")
        }
    }
    res=as.double(res)
    countMatrix <- countMatrix[,sort(colnames(countMatrix))]
    tissuePosition <- as.matrix(read.table(paste("./../",tissuePositionFile,sep=""),header=TRUE,sep="\t",row.names=1))
    d <- dim(tissuePosition)[2]
    tissuePosition <- tissuePosition[,(d-1):d]
    tissuePosition <- tissuePosition[sort(rownames(tissuePosition)),]

    pbmc <- CreateSeuratObject(countMatrix)
    pbmc <- SCTransform(pbmc, assay = "RNA", verbose = FALSE)
    pbmc <- RunPCA(pbmc, assay = "SCT", verbose = FALSE)
    if(pcaDimensions <= 2){
        pcaDimensions = 2
    }
    m <- pbmc@reductions[["pca"]]@cell.embeddings[,1:pcaDimensions]
    distPCA = dist(m,method="minkowski",p=as.numeric(profileDistance))  
    distCoord <- dist(tissuePosition,method="minkowski",p=as.numeric(spotDistance))
    distCoord <- distCoord*((max(distPCA)*as.double(spotDistanceTransformationWeight))/(max(distCoord)))
    finalDistance <- as.matrix(distPCA + distCoord)
    neighbors <- FindNeighbors(finalDistance)
    neighbors <- list(neighbors_nn=neighbors$nn,neighbors_snn=neighbors$snn)
    pbmc@graphs <- neighbors
    pbmc <- FindClusters(pbmc,dims.use = 1:pcaDimensions, resolution=res,verbose = FALSE, graph.name = "neighbors_snn")

    pbmc.new <- RunTSNE(object = distPCA + distCoord)
    mainVector <- as.numeric(pbmc@active.ident) 
    Coordinates <- pbmc.new@cell.embeddings

    nCluster <- max(mainVector)
    dir.create(paste("./",nCluster,sep=""))
    dir.create(paste("./",nCluster,"/Permutation",sep=""))
    setwd(paste("./",nCluster,sep=""))

    clustering.output <- cbind(rownames(Coordinates),mainVector,Coordinates[,1],Coordinates[,2])
    clustering.output <- silhouette(length(unique(mainVector)),clustering.output)
    colnames(clustering.output) <- c("cellName","Belonging_Cluster","xChoord","yChoord","extraScore","intraScore","neighbor","silhouetteValue")
    write.table(clustering.output,paste(matrixName,"_clustering.output.",format,sep=""),sep=separator2, row.names = F)
    cycles <- nPerm/permAtTime
    cat(getwd())
    for(i in 1:cycles){
            system(paste("for X in $(seq ",permAtTime,")
        do
        nohup Rscript ./../../../home/permutation.R ",percent," ",matrixName," ",tissuePositionFile," ",profileDistance," ",spotDistance," ",spotDistanceTransformationWeight," ",res," ",format," ",separator," ",logTen," ",pcaDimensions," ",sparse," $(($X +",(i-1)*permAtTime," )) & 

        done"))
        d=1
        while(length(list.files("./Permutation",pattern=paste("*.",format,sep="")))!=i*permAtTime*2){
            if(d==1){cat(paste("Cluster number ",nCluster," ",((permAtTime*i))/nPerm*100," % complete \n"))}
            d=2
        }
        system("echo 3 > /proc/sys/vm/drop_caches")
        system("sync")
        gc()
    }
    cluster_p <- sapply(list.files("./Permutation/",pattern="cluster*"),FUN=function(x){a=read.table(paste("./Permutation/",x,sep=""),header=TRUE,col.names=1,sep=separator2)[[1]]})
    killedC <- sapply(list.files("./Permutation/",pattern="killC*"),FUN=function(x){a=read.table(paste("./Permutation/",x,sep=""),header=TRUE,col.names=1,sep=separator2)[[1]]})

    write.table(as.matrix(cluster_p,col.names=1),paste(matrixName,"_",nCluster,"_clusterP.",format,sep=""),sep=separator2,row.names=FALSE, quote=FALSE)
    write.table(as.matrix(killedC,col.names=1),paste(matrixName,"_",nCluster,"_killedCell.",format,sep=""),sep=separator2,row.names=FALSE, quote=FALSE)

    pdf("hist.pdf")
    clusters <- apply(cluster_p,2,FUN=function(x){max(x)})
    hist(clusters,xlab="nCluster",breaks=length(unique(cluster_p)))
    dev.off()

    write.table(sort(unique(clusters)),paste("./../rangeVector.",format,sep=""),sep=separator2,row.names=FALSE,col.names=FALSE)
    system("rm -r Permutation")
    return(length(unique(mainVector)))
}

silhouettePlot=function(matrixName,rangeVector,format,separator){
    if(separator=="tab"){separator="\t"} #BUG CORRECTION TAB PROBLEM 
    count=1
    l=list()
    for(i in rangeVector){
        l[[count]]=read.table(paste("./",i,"/",matrixName,"_clustering.output.",format,sep=""),sep=separator,header=TRUE)[,8]
        count=count+1
    }
    pdf(paste(matrixName,"_vioplot.pdf",sep=""))
    do.call(vioplot,c(l,list(names=rangeVector)))
    dev.off()
}