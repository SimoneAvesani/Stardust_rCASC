# Stardust installation and usage
Stardust can be installed as a standalone R package or can be used through the dedicated docker image. We suggest installing the following tools on a UNIX-like OS (like MacOS or a Linux distribution). The main aim of the tool is, given as input an expression matrix, the positions of spots and a space weight configuration, to derive a vector of cluster identities for each spot in the input data. Stardust depends on Seurat (for clustering and data visualization) and on rCASC for the computation of the stability scores and generation of the distributions comparison. Stardust, Seurat and rCASC installation instructions and usage are reported in the following subsections. Stardust tuning is reported in the last subsection. 

## Standalone R package
Stardust depends on Seurat, so we first need to install this package through the remotes package (in this way we can install a fixed version of Seurat). After Seurat installation we can install Stardust from our GitHub repository.
```R
# R code
install.packages("devtools")
install.packages("remotes")

# install Seurat
remotes::install_version("Seurat", version = "3.2.2")
# install Stardust
devtools::install_github("InfOmics/stardust")

```
Once the packages are installed, you can execute the following sample workflow based on the Mouse Kidney dataset. Let us download the input data.
```bash
# Bash code
# create a working directory and enter in it
mkdir MouseKidney && cd MouseKidney
# using wget download the expression matrix and spot positions of the Mouse Kidney            
# dataset. Download also the full dataset for the creation of a Seurat object for 
# visualization purposes

wget https://github.com/SimoneAvesani/Stardust_rCASC/raw/master/Datasets/MouseKidney/filtered_expression_matrix.txt.zip

wget https://github.com/SimoneAvesani/Stardust_rCASC/raw/master/Datasets/MouseKidney/spot_coordinates.txt

wget https://github.com/SimoneAvesani/Stardust_rCASC/raw/master/Datasets/MouseKidney/FullDataset.zip

# unzip the archives and delete unused data
unzip filtered_expression_matrix.txt.zip
unzip FullDataset.zip
rm -rf __MACOSX
rm filtered_expression_matrix.txt.zip FullDataset.zip

# start R
R
```
Now compute the cluster identities for each spot.
```R
# R code
# load Seurat and Stardust
library("Seurat")
library("stardust")

# load the count matrix and spot coordinates for the Mouse Kidney dataset
countMatrix = read.table("./filtered_expression_matrix.txt",row.names=1,header = TRUE)

spotPositions = read.table("./spot_coordinates.txt",row.names=1,header = TRUE)

# execute stardust passing to the method the count matrix, spot position and 
# the weight of spatial information relative to the transcriptional 
# similarity (spaceWeight can be a real number between 0 and 1)
output <- StardustOnSeurat(countMatrix = countMatrix, spotPositions = spotPositions, spaceWeight = 0.75)

# get the vector of cluster identities for each spot
clusters_identities = output@active.ident
```
Finally, load the full dataset (i.e. with histological images) published on 10X website (10X, Datasets) and assign to it the cluster identities. This object can then be visualized with Seurat.
```R
# R code
# create a full Seurat object with the data already downloaded
MouseKidney = Load10X_Spatial("./FullDataset/")

# assign cluster identities to the Seurat object
MouseKidney@active.ident = clusters_identities

# visualize the clusters overlaid to the tissue image
Seurat::SpatialDimPlot(MouseKidney)
```
### Standalone R package with docker
If you want a straight forward usage of Stardust you can also pull the dedicated docker container and skip all the possible dependency problems you could encounter with the package installation:
```bash
# Bash code
# First, pull the docker image and run it
docker pull giovannics/stardust
docker run --rm -it giovannics/stardust /bin/bash

# create a working directory and enter in it
mkdir MouseKidney && cd MouseKidney

# Using wget, download the count matrix and spot coordinates 
# (plus the full dataset for visualization purposes)
wget https://github.com/SimoneAvesani/Stardust_rCASC/raw/master/Datasets/MouseKidney/FullDataset.zip

wget https://github.com/SimoneAvesani/Stardust_rCASC/raw/master/Datasets/MouseKidney/filtered_expression_matrix.txt.zip

wget https://github.com/SimoneAvesani/Stardust_rCASC/raw/master/Datasets/MouseKidney/spot_coordinates.txt


# unzip the archives and delete unused data
unzip filtered_expression_matrix.txt.zip
unzip FullDataset.zip
rm -rf __MACOSX
rm filtered_expression_matrix.txt.zip FullDataset.zip

# start R
R
```

```R
# R code
# load Seurat and Stardust
library("Seurat")
library("stardust")

# load the count matrix and spot coordinates for the Mouse Kidney dataset
countMatrix = read.table("./filtered_expression_matrix.txt",row.names=1,header = TRUE)

spotPositions = read.table("./spot_coordinates.txt",row.names=1,header = TRUE)

# execute stardust passing to the method the count matrix, spot position and 
# the weight of spatial information relative to the transcriptional 
# similarity (spaceWeight can be a real number between 0 and 1)
output <- StardustOnSeurat(countMatrix = countMatrix, spotPositions = spotPositions, spaceWeight = 0.75)

# get the vector of cluster identities for each spot
clusters_identities = output@active.ident

# you can save the cluster identities to export them outside if you need them
write.table(clusters_identities,file="clusters_identities.txt")

# load the full dataset as a Seurat object and overwrite the cluster identities
MouseKidney = Load10X_Spatial("./FullDataset/")
MouseKidney@active.ident = factor(clusters_identities)

# save the plot as a png to export outside the container
png("spatialClustersPlot.png", units="px", width=800, height=800, res=150)
# alternatively you can save the plot as a jpg
jpeg('spatialClustersPlot.jpg')

Seurat::SpatialDimPlot(MouseKidney)
dev.off()
```

```bash
# Bash code (on a new terminal)
# get the container id
docker ps

# extract the cluster identities and clusters plot from the container
docker cp container_id:/MouseKidney/clusters_identities.txt .
docker cp container_id:/MouseKidney/spatialClustersPlot.png .

# you can now inspect in your local machine the clusters id that Stardust assigned to each spot and the clusters plot.
```

## Stability scores computation with rCASC
rCASC stability scores computation allows users to evaluate which Stardust configuration performs better on your particular dataset. rCASC is designed with a container architecture so that it can provide computational reproducibility across different machines. The package can be installed from our fork on GitHub and the required docker images can be pulled from Docker Hub.

```bash
# Bash code
# pull the image to run the permutations of Stardust 
# on multiple permutated datasets
docker pull giovannics/spatial2020seuratpermutation

# pull the image to run the stability scores computation
docker pull repbioinfo/seuratanalysis

# prepare a dedicated directory and download the count matrix and spot positions for 
# the Mouse Kidney dataset (as an example)
mkdir -p MouseKidney/scratch && cd MouseKidney

wget https://github.com/SimoneAvesani/Stardust_rCASC/raw/master/Datasets/MouseKidney/filtered_expression_matrix.txt.zip

wget https://github.com/SimoneAvesani/Stardust_rCASC/raw/master/Datasets/MouseKidney/spot_coordinates.txt

# unzip the archive and delete unused data
unzip filtered_expression_matrix.txt.zip
rm -rf __MACOSX
rm filtered_expression_matrix.txt.zip

# start R
R
```

```R
# R code
# install our rCASC fork on GitHub through the R package devtools
library(devtools)
install_github("InfOmics/rCASC")

# install ggplot that is a dependency for the figure generation
install.packages("ggplot2")
library("ggplot2")
```
When your dependencies are installed you can generate the violin plot image. In this way you can explore which configuration works best for your dataset by varying the parameters. 

If you want to evaluate the stability of one Stardust configuration you can call the StardustPermutation and the permAnalysisSeurat methods. 

```R
# load rCASC
library(rCASC)

# set the variables that contain the paths for the temporary files folder of rCASC, 
# the count matrix and the spot positions file
scratch.folder <- paste(getwd(),"/scratch",sep="")
file <- paste(getwd(),"/filtered_expression_matrix.txt",sep="")
tissuePosition <- paste(getwd(),"/spot_coordinates.txt",sep="")

# call the rCASC method to perform the permutations of a particular space
# configuration. The parameters meaning are:
# group ??? to create the docker image without superuser privileges
# scratch.folder ??? path of the folder that rCASC use for storing temporary files
# file ??? path of the count matrix file
# tissuePosition ??? path of the spot coordinates file
# spaceWeight ??? real number between 0 and 1 that describe how much space weight if
#               compared to the transcriptional similarity
# res ??? clustering resolution 
# nPerm ??? number of permutations to be computed
# permAtTime ??? number of permutation to compute in parallel
# percent ??? percentage of the input dataset to remove for each permutation
# pcaDimensions ??? number of principal components 
# separator ??? character separator of values in the input files

StardustPermutation(group="docker", scratch.folder=scratch.folder, file=file, tissuePosition=tissuePosition, spaceWeight=0.75, 
res=0.8, nPerm=80, permAtTime=8, percent=10, pcaDimensions=10, separator="\t")

# extract the number of clusters obtained in order to configure the next method call
cluster.path <- paste(data.folder=dirname(file), "Results", strsplit(basename(file),"\\.")[[1]][1], sep="/")
cluster <- as.numeric(list.dirs(cluster.path, full.names = FALSE, recursive = FALSE))

# call permAnalysisSeurat in order to compute the stability scores based on the 
# previous permutations. The parameters meaning are:
# group ??? to create the docker image without superuser privileges
# scratch.folder ??? path of the folder that rCASC use for storing temporary files
# file ??? path of the count matrix file
# nCluster ??? number of cluster obtained before
# separator ??? character separator of values in the input files
# sp ??? minimum number of percentage of cells that has to be in common between two 
#      permutation to be the same cluster.

permAnalysisSeurat(group="docker", scratch.folder = scratch.folder, file=file, nCluster=cluster, separator="\t", sp=0.8)
``` 
In ???Results/filtered_expression_matrix/9/filtered_expression_matrix_clustering.output.txt??? file, you will find the assigned cluster identity of each spot, and in ???Results/filtered_expression_matrix/9/filtered_expression_matrix_scoreSum.txt??? file its stability score for the configuration used (spaceWeight=0.75).

The coefficient of variation value can be computed from the "filtered_expression_matrix_scoreSum.txt??? file as follows. 

```R
# R code
# Read the stability scores and compute the coefficient of variation
mat <- read.table("filtered_expression_matrix_scoreSum.txt")
cv <- sd(mat$V2)/mean(mat$V2)
cv
``` 
For each compared method, a dedicated container image can be pulled to run the permutations. Note that each method requires specific input data that must be prepared in advance, see https://github.com/SimoneAvesani/Stardust_rCASC/tree/master/Tools_Comparison/homes for more details. 

```bash
# Bash code
# pull the image for each method
docker pull giovannics/bayespacepermutation
docker pull giovannics/giottopermutation
docker pull giovannics/spagcnpermutation
docker pull giovannics/stlearn-rcasc 
# start R
R
```

```R
# R code
library(rCASC) 

# For BayesSpace use
bayeSpacePermutation(group="docker", scratch.folder=scratch.folder, file=file, filtered_feature_bc_matrix=filtered_feature_bc_matrix, 
n_clusters=n_clusters, spatial=spatial,nPerm=80, permAtTime=8)

# For Giotto use
GiottoPermutation(group="docker", scratch.folder=scratch.folder, file=file, h5matrix.name=h5matrix.name, 
spotpositions.name=spotpositions.name, n_clusters=n_clusters, pcaDimensions=pcaDimensions, nPerm=80, permAtTime=8, percent=10)

# For SpaGCN use
spaGCNPermutation(group="docker", scratch.folder=scratch.folder, h5matrix.name=h5matrix.name, spotpositions.name=spotpositions.name, 
image.name=image.name, use_histology=TRUE, lResolution=lResolution, pcaDimensions=pcaDimensions, nPerm=80, permAtTime=8)

# For stLearn use
STLearnPermutation(group="docker", scratch.folder=scratch.folder, file=file, filtered_feature_bc_matrix=filtered_feature_bc_matrix, 
lResolution=res, nPerm=80, permAtTime=8,percent=10, pcaDimensions=pcaDimensions)

# For each method, extract the number of clusters obtained in order to configure 
# the next method call as in the workflow above
cluster.path <- paste(data.folder=dirname(file), "Results", strsplit(basename(file),"\\.")[[1]][1], sep="/")
cluster <- as.numeric(list.dirs(cluster.path, full.names=FALSE, recursive=FALSE))

# call permAnalysisSeurat in order to compute the stability scores based on the 
# previous permutations. The parameters meaning are:
# group ??? to create the docker image without superuser privileges
# scratch.folder ??? path of the folder that rCASC use for storing temporary files
# file ??? path of the count matrix file
# nCluster ??? number of cluster obtained before
# separator ??? character separator of values in the input files
# sp ??? minimum number of percentage of cells that has to be in common between two 
#      permutation to be the same cluster.

permAnalysisSeurat(group="docker", scratch.folder=scratch.folder, file=file, nCluster=cluster, separator="\t", sp=0.8)
``` 

## Stardust tuning 
First, pull the docker image and run it. Then, navigate into the runExample directory and prepare a directory to download the desired datasets. 

```bash
# Bash code
docker pull repbioinfo/tunestardust
docker run -it --privileged=true repbioinfo/tunestardust /bin/bash

# download the count matrix and spot positions for
# the Mouse Kidney dataset (as an example)
mkdir -p Datasets/MK/scratch && cd Datasets/MK

wget https://github.com/SimoneAvesani/Stardust_rCASC/raw/master/Datasets/MouseKidney/filtered_expression_matrix.txt.zip

wget https://github.com/SimoneAvesani/Stardust_rCASC/raw/master/Datasets/MouseKidney/spot_coordinates.txt

# unzip the archive and delete unused data
unzip filtered_expression_matrix.txt.zip
rm -rf __MACOSX
rm filtered_expression_matrix.txt.zip
``` 
From the runExample directory execute the runTuning.R script passing as an argument the dataset's folder name.

```bash
# Bash code
Rscript runTuning.R MK
```
In "Datasets/MK/results.txt", you will find the best space weight and resolution parameters returned. 
