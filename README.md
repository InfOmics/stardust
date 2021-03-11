# Stardust installation and usage
Stardust can be installed as a standalone R package in your local machine or can be used through the dedicated docker image. The only purpose of the tool is, given as input an expression matrix, the positions of spots and a space weight configuration, to return a vector of cluster identities for each spot in the input data. If you want to visualize the results, you need Seurat R package and a full Seurat object of the dataset (i.e. with tissue images) in which you can overwrite the clusters assignments. Moreover, if you want to compute the stability scores of a particular configuration or if you want to generate the distributions comparison as in Figure S1(b) of "Supplementary data of “Stardust: improving spatial tranScripTomics data analysis through space awARe modularity optimization baseD clUSTering.”, you need to install our rCASC (Alessandrì et al., 2019) fork from GitHub. 

## Standalone R package
Before installing the package, let’s download some data that we will use in the workflow.
```bash
# create a working directory and enter in it
mkdir MouseKidney && cd MouseKidney
# using wget download the expression matrix and spot positions of the Mouse Kidney dataset. Download also the full dataset for the creation of a Seurat object for visualization purposes

wget https://github.com/InfOmics/stardust/raw/validation_data/stardustData/Datasets/MouseKidney/filtered_expression_matrix.txt.zip

wget https://raw.githubusercontent.com/InfOmics/stardust/validation_data/stardustData/Datasets/MouseKidney/spot_coordinates.txt

wget https://github.com/InfOmics/stardust/raw/validation_data/stardustData/Datasets/MouseKidney/FullDataset.zip

# unzip the archives and delete unused data
unzip filtered_expression_matrix.txt.zip
unzip FullDataset
rm -rf __MACOSX
rm filtered_expression_matrix.txt.zip FullDataset.zip

# start R
R
```

Stardust is dependent from Seurat, so we first need to install this package thorough devtools (in this way we can install a fixed version of Seurat):
```R
install.packages("devtools")
library("devtools")
# install Seurat
devtools::install_version(package = 'Seurat', version = 
   package_version('3.2.1'))
# install Stardust
devtools::install_github("InfOmics/stardust")
```
Once the package is installed, you can execute this sample workflow based on the Mouse Kidney dataset:
```R
# load Seurat and Stardust
library("Seurat")
library("stardust")

# load the count matrix and spot coordinates for the Mouse Kidney dataset
countMatrix = read.table("./filtered_expression_matrix.txt",row.names=1,header = TRUE)

spotPositions = read.table("./spot_coordinates.txt",row.names=1,header = TRUE)

# execute stardust passing to the method the count matrix, spot position and the weight of spatial information relative to the transcriptional similarity (spaceWeight can be a real number between 0 and 1)
output <- StardustOnSeurat(countMatrix = MouseKidney, spotPositions = MouseKidneyCoord, spaceWeight = 0.75)

# get the vector of cluster identities for each spot
clusters_identities = output@active.ident
```
If you have access to a full Seurat object (i.e. with histological images) of the dataset, you can overwrite the cluster identities slot and visualize the result with Seurat. We downloaded the necessary data in the first step.

```R
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
# First download the full dataset for the creation of a Seurat object for visualization purposes

wget https://github.com/InfOmics/stardust/raw/validation_data/stardustData/Datasets/MouseKidney/FullDataset.zip
unzip FullDataset

# pull the docker image and run it
docker pull giovannics/stardust
docker run -it giovannics/stardust /bin/bash
# You can now download the count matrix and spot coordinates and run Stardust without the need to install it. 

# create a working directory and enter in it
mkdir MouseKidney && cd MouseKidney

# using wget download the expression matrix and spot positions of the Mouse Kidney dataset.
wget https://github.com/InfOmics/stardust/raw/validation_data/stardustData/Datasets/MouseKidney/filtered_expression_matrix.txt.zip

wget https://raw.githubusercontent.com/InfOmics/stardust/validation_data/stardustData/Datasets/MouseKidney/spot_coordinates.txt


# unzip the archives and delete unused data
unzip filtered_expression_matrix.txt.zip
rm -rf __MACOSX
rm filtered_expression_matrix.txt.zip 

# start R
R
```

```R
# load Seurat and Stardust
library("Seurat")
library("stardust")

# load the count matrix and spot coordinates for the Mouse Kidney dataset
countMatrix = read.table("./filtered_expression_matrix.txt",row.names=1,header = TRUE)

spotPositions = read.table("./spot_coordinates.txt",row.names=1,header = TRUE)

# execute stardust passing to the method the count matrix, spot position and the weight of spatial information relative to the transcriptional similarity (spaceWeight can be a real number between 0 and 1)
output <- StardustOnSeurat(countMatrix = MouseKidney, spotPositions = MouseKidneyCoord, spaceWeight = 0.75)

# get the vector of cluster identities for each spot
clusters_identities = output@active.ident

# you can save the cluster identities to export them outside 
write.table(clusters_identities,file="clusters_identities.txt")
```

If you want to visualize the result over the tissue image you need to create a full Seurat object of the dataset outside the container ( you need Seurat installed on your machine) and extract the cluster identities. 
```bash
# On a new terminal..
# get the container id
docker ps
# extract the data from the container
docker cp container_id:/MouseKidney/clusters_identities.txt .

# start R
R
```

```R
# load Seurat
library("Seurat")

# create the full Seurat object
MouseKidney = Load10X_Spatial("./FullDataset/")

#load the cluster identities
clusters_identities = read.table("./clusters_identities.txt",row.names=1,header = TRUE)

# overwrite cluster identities in the Seurat object
MouseKidney@active.ident = factor(clusters_identities$x)

# visualize the clusters overlaid to the tissue image
Seurat::SpatialDimPlot(MouseKidney)
```

## Stability scores computation with rCASC
rCASC stability scores computation allows exploring which Stardust configuration performs better on your particular dataset. You first need to install our rCASC package fork from GitHub and then to download the required docker containers. rCASC is designed with a containers architecture so that it can provide computational reproducibility across different machines:

```bash
# pull the container that implement the permutations of Stardust on multiple permutated datasets
docker pull giovannics/spatial2020seuratpermutation

# pull the container that implement the stability scores computation
docker pull repbioinfo/seuratanalysis

# start R
R
```


```R
# install out rCASC fork on GitHub through the R package devtools
install.packages("devtools")
library(devtools)
install_github("InfOmics/rCASC")
```

When your dependencies are installed you can generate the violin plots comparison figure as in Figure S1(b). In this way you can explore which configuration works best for your dataset.
```bash
# prepare a dedicated directory and download the count matrix and spot positions for the Mouse Kidney dataset
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
# set the variables that contain the paths for the temporary files folder of rCASC, the count matrix and the spot positions file
scratch.folder <- paste(getwd(),"/scratch",sep="")
file <- paste(getwd(),"/filtered_expression_matrix.txt",sep="")
tissuePosition <- paste(getwd(),"/spot_coordinates.txt",sep="")

# call the rCASC method for generation of the violin plot comparison figure.
# It repeats the clustering permutation for 5 space configurations of Stardust (0,0.25,0.5,0.75 and 1). The parameters meaning are:
# scratch.folder → path of the folder that rCASC use for storing temporary files
# file → path of the count matrix file
# tissuePosition → path of the spot coordinates file
# nPerm → number of permutations to be computed
# permAtTime → number of permutation to compute in parallel
# percent → percentage of the input dataset to remove for each permutation
# separator → character separator of values in the input files

StartdustConfigurations(scratch.folder=scratch.folder, file=file, 
tissuePosition=tissuePosition, nPerm=80, permAtTime=8, percent=10, separator"\t")
```
Under the “MouseKidney” folder you will see the figure produced and all the data used to create it.
If you want to evaluate the stability of only one configuration (that is less computationally expensive) you can switch the StartdustConfigurations method call with the following one.
```R
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
