# Stardust installation and usage
Stardust can be installed as a standalone R package or can be used through the dedicated docker image. We suggest installing the following tools on a UNIX-like OS (like MacOS or a Linux distribution). The main aim of the tool is, given as input an expression matrix, the positions of spots and a space weight configuration, to derive a vector of cluster identities for each spot in the input data. Stardust depends on Seurat (for clustering and data visualization) and on rCASC for the computation of the stability scores and generation of the distributions comparison as in Figure S1(b). Stardust, Seurat and rCASC installation instructions are reported in the following subsections. 

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
# using wget download the expression matrix and spot positions of the Mouse Kidney            # dataset. Download also the full dataset for the creation of a Seurat object for 
# visualization purposes

wget https://github.com/InfOmics/stardust/raw/validation_data/stardustData/Datasets/MouseKidney/filtered_expression_matrix.txt.zip

wget https://raw.githubusercontent.com/InfOmics/stardust/validation_data/stardustData/Datasets/MouseKidney/spot_coordinates.txt

wget https://github.com/InfOmics/stardust/raw/validation_data/stardustData/Datasets/MouseKidney/FullDataset.zip

# unzip the archives and delete unused data
unzip filtered_expression_matrix.txt.zip
rm -rf __MACOSX
rm filtered_expression_matrix.txt.zip FullDataset.zip

unzip FullDataset
rm FullDataset.zip

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

# Using wget, download the count matrix and spot coordinates (plus the full dataset for visualization purposes)
wget https://github.com/InfOmics/stardust/raw/validation_data/stardustData/Datasets/MouseKidney/FullDataset.zip

wget https://github.com/InfOmics/stardust/raw/validation_data/stardustData/Datasets/MouseKidney/filtered_expression_matrix.txt.zip

wget https://raw.githubusercontent.com/InfOmics/stardust/validation_data/stardustData/Datasets/MouseKidney/spot_coordinates.txt


# unzip the archives and delete unused data
unzip filtered_expression_matrix.txt.zip
rm -rf __MACOSX
rm filtered_expression_matrix.txt.zip 

unzip FullDataset
rm FullDataset.zip

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
rCASC stability scores computation allows users to evaluate which Stardust configuration performs better on your particular dataset. rCASC is designed with a container architecture so that it can provide computational reproducibility across different machines. The package can be installed from our fork on GitHub and the required docker containers can be pulled from Docker Hub.

```bash
# Bash code
# pull the container that implement the permutations of Stardust on multiple
# permutated datasets
docker pull giovannics/spatial2020seuratpermutation

# pull the container that implement the stability scores computation
docker pull repbioinfo/seuratanalysis

# start R
R
```


```R
# R code
# install out rCASC fork on GitHub through the R package devtools
install.packages("devtools")
library(devtools)
install_github("InfOmics/rCASC")
```
When your dependencies are installed you can generate the violin plots comparison image as in Figure S1(b). In this way you can explore which configuration works best for your dataset.
```bash
# Bash code
# prepare a dedicated directory and download the count matrix and spot positions for 
# the Mouse Kidney dataset
mkdir -p MouseKidney/scratch && cd MouseKidney

wget https://github.com/GiovanniCS/StardustData/raw/main/Datasets/MouseKidney/filtered_expression_matrix.txt.zip

wget https://raw.githubusercontent.com/GiovanniCS/StardustData/main/Datasets/MouseKidney/spot_coordinates.txt

# unzip the archive and delete unused data
unzip filtered_expression_matrix.txt.zip
rm -rf __MACOSX
rm filtered_expression_matrix.txt.zip

# start R
R
```

```R
# R code
# install ggplot that is a dependency for the figure generation
install.packages("ggplot2")
library("ggplot2")

# set the variables that contain the paths for the temporary files folder of rCASC,
# the count matrix and the spot positions file
scratch.folder <- paste(getwd(),"/scratch",sep="")
file <- paste(getwd(),"/filtered_expression_matrix.txt",sep="")
tissuePosition <- paste(getwd(),"/spot_coordinates.txt",sep="")

# call the rCASC method for generation of the violin plot comparison figure.
# It repeats the clustering permutation for 5 space configurations of Stardust (0,0.25,0.5,0.75 and 1). The parameters meaning are:
# group → to create the docker image without superuser privileges
# scratch.folder → path of the folder that rCASC use for storing temporary files
# file → path of the count matrix file
# tissuePosition → path of the spot coordinates file
# nPerm → number of permutations to be computed
# permAtTime → number of permutation to compute in parallel
# percent → percentage of the input dataset to remove for each permutation
# separator → character separator of values in the input files

StartdustConfigurations(group="docker",scratch.folder=scratch.folder,
file=file, tissuePosition=tissuePosition, nPerm=80, permAtTime=8, percent=10, separator="\t")
```
Under the “MouseKidney” folder you will see the figure produced and all the data used to create it. If you want to evaluate the stability of only one configuration (that is less computationally expensive) you can switch the StartdustConfigurations method call with the following one.
```R
# R code
# you need to have rCASC already installed, the containers and data in your current location of the file system as in the workflow above

# load rCASC
library(rCASC)

# set the variables that contain the paths for the temporary files folder of rCASC, the count matrix and the spot positions file
scratch.folder <- paste(getwd(),"/scratch",sep="")
file <- paste(getwd(),"/filtered_expression_matrix.txt",sep="")
tissuePosition <- paste(getwd(),"/spot_coordinates.txt",sep="")


# call the rCASC method to perform the permutations of a particular space
# configuration. The parameters meaning are:
# group → to create the docker image without superuser privileges
# scratch.folder → path of the folder that rCASC use for storing temporary files
# file → path of the count matrix file
# tissuePosition → path of the spot coordinates file
# spaceWeight → real number between 0 and 1 that describe how much space weight if
#               compared to the transcriptional similarity
# nPerm → number of permutations to be computed
# permAtTime → number of permutation to compute in parallel
# percent → percentage of the input dataset to remove for each permutation
# separator → character separator of values in the input files

StardustPermutation(group="docker",scratch.folder = scratch.folder,
file=file, tissuePosition=tissuePosition, spaceWeight=0.75, nPerm=80, permAtTime=8, percent=10, separator="\t")

# extract the number of clusters obtained in order to configure the next method call
cluster.path <- paste(data.folder=dirname(file), "Results", strsplit(basename(file),"\\.")[[1]][1], sep="/")
cluster <- as.numeric(list.dirs(cluster.path, full.names = FALSE, recursive = FALSE))

# call permAnalysisSeurat in order to compute the stability scores based on the 
# previous permutations. The parameters meaning are:
# group → to create the docker image without superuser privileges
# scratch.folder → path of the folder that rCASC use for storing temporary files
# file → path of the count matrix file
# nCluster → number of cluster obtained before
# separator → character separator of values in the input files
# sp → minimum number of percentage of cells that has to be in common between two 
#      permutation to be the same cluster.
permAnalysisSeurat(group="docker",scratch.folder = scratch.folder,file=file, nCluster=cluster,separator="\t",sp=0.8)
``` 
In “Results/filtered_expression_matrix/10/filtered_expression_matrix_clustering.output.txt” file, you will find the assigned cluster identity of each spot, and in “Results/filtered_expression_matrix/10/filtered_expression_matrix_scoreSum.txt” file its stability score for the configuration used (spaceWeight=0.75).