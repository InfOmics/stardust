# Stardust
Stardust is an R package that allows to integrate positional information with transcription similarity when clustering Spatial Transcriptomics data.  
It is clustering algorithm agnostic as it can work with any method that accepts as input a distance matrix. Currently the only method supported is the one provided by Seurat (Louvain algorithm).  

In order to setup the environment you need to install Seurat and devtools
```R
install.packages("devtools")
library("devtools")
devtools::install_version(package = 'Seurat', version = package_version('3.2.1'))
```
Then, using devtools, install stardust package from its github reporitory
```R
install_github("InfOmics/stardust")
```  
Here is an example of workflow:
```R
library("stardust")
data("MouseKidney")
data("MouseKidneyCoord")
seurat.object <- StardustOnSeurat(countMatrix = MouseKidney, spotPositions = MouseKidneyCoord, spaceWeight = 0.75)
plot(MouseKidneyCoord, col = seurat.object@active.ident)
``` 