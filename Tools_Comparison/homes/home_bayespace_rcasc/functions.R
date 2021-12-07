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

clustering=function(matrixName,n_clusters,pcaDimensions,nPerm,permAtTime,
    percent,nCluster){
    n_clusters = as.integer(n_clusters)
    pcaDimensions = as.integer(pcaDimensions)
    countMatrix=read.table(paste("/scratch/",matrixName,sep=""),sep="\t",header=TRUE,row.name=1)
    countMatrix <- countMatrix[,sort(colnames(countMatrix))]
    pbmc <- CreateSeuratObject(countMatrix)
    pbmc <- SCTransform(pbmc, assay = "RNA", verbose = FALSE)
    pbmc <- RunPCA(pbmc, assay = "SCT", verbose = FALSE)
    pbmc.new <- RunTSNE(object = pbmc) 
    Coordinates <- pbmc.new@reductions[["tsne"]]@cell.embeddings

    sce <- readVisium("/scratch/")
    sce = sce[,sort(colnames(sce))]
    sce = sce[rowSums(counts(sce)) > 10,]
    sce <- spatialPreprocess(sce, platform="ST", n.PCs=pcaDimensions, n.HVGs=2000, log.normalize=TRUE)
    sce <- spatialCluster(sce, q=n_clusters, platform="ST", d=pcaDimensions,
                           init.method="mclust", model="t", gamma=2,
                           nrep=10000, burn.in=100,
                           save.chain=TRUE)
    system("mkdir /scratch/bayespaceout")
    png('/scratch/bayespaceout/bayespace.png')
    clusterPlot(sce)
    dev.off()
    mainVector = colData(sce)$spatial.cluster
    mainVector = as.numeric(mainVector)
    jumping_clusters = sort(unique(mainVector))
    for(i in 1:length(jumping_clusters)){
        mainVector[mainVector==jumping_clusters[i]] = i
    }
    nCluster <- max(mainVector)
    dir.create(paste("./",nCluster,sep=""))
    dir.create(paste("./",nCluster,"/Permutation",sep=""))
    setwd(paste("./",nCluster,sep=""))

    clustering.output <- cbind(rownames(Coordinates),mainVector,Coordinates[,1],Coordinates[,2])
    clustering.output <- silhouette(length(unique(mainVector)),clustering.output)
    colnames(clustering.output) <- c("cellName","Belonging_Cluster","xChoord","yChoord","extraScore","intraScore","neighbor","silhouetteValue")
    matrixNameBis = strsplit(matrixName,".",fixed = TRUE)[[1]][1]
    write.table(clustering.output,paste(matrixNameBis,"_clustering.output.","txt",sep=""),sep="\t", row.names = F)
    cycles <- nPerm/permAtTime
    cat(getwd())
    for(i in 1:cycles){
            system(paste("for X in $(seq ",permAtTime,")
        do
        nohup Rscript ./../../../home/permutation.R ",percent," ",n_clusters," ",pcaDimensions," $(($X +",(i-1)*permAtTime," )) & 

        done"))
        d=1
        while(length(list.files("./Permutation",pattern=paste("*.","txt",sep="")))!=i*permAtTime*2){
            if(d==1){cat(paste("Cluster number ",nCluster," ",((permAtTime*i))/nPerm*100," % complete \n"))}
            d=2
        }
        system("echo 3 > /proc/sys/vm/drop_caches")
        system("sync")
        gc()
    }
    cluster_p <- sapply(list.files("./Permutation/",pattern="cluster*"),FUN=function(x){a=read.table(paste("./Permutation/",x,sep=""),header=TRUE,col.names=1,sep="\t")[[1]]})
    killedC <- sapply(list.files("./Permutation/",pattern="killC*"),FUN=function(x){a=read.table(paste("./Permutation/",x,sep=""),header=TRUE,col.names=1,sep="\t")[[1]]})

    write.table(as.matrix(cluster_p,col.names=1),paste(matrixNameBis,"_",nCluster,"_clusterP.","txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)
    write.table(as.matrix(killedC,col.names=1),paste(matrixNameBis,"_",nCluster,"_killedCell.","txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

    pdf("hist.pdf")
    clusters <- apply(cluster_p,2,FUN=function(x){max(x)})
    hist(clusters,xlab="nCluster",breaks=length(unique(cluster_p)))
    dev.off()

    write.table(sort(unique(clusters)),paste("./../rangeVector.","txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE)
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