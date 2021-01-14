
#' @title Stardust on Seurat
#' @description This function calls Seurat clustering feature, passing a custom distance matrix that integrate both
#' transcriptional similarity and spatial information
#' @param countMatrix, dataframe with spots id as columns and genes id as rows.
#' @param spotPositions, dataframe with x and y spot coordinates as columns and spots id as rows.
#' @param spaceWeight, double in [0,1]. Weight for the linear transformation of spot distance.
#'  1 means spot distance weight as much as profile distance, 0 means spot distance doesn't contribute at all in the overall 
#'  distance measure.
#' @param pcaDimensions, number of PCA dimension to include in the distance computation.
#' 
#' @author Giovanni Motterle, giovanni [dot] motterle [at] univr [dot] it, University of Verona
#'
#' @return Seurat object with cluster identifications for each spot under active.ident slot.
#' @export

StardustOnSeurat <- function(countMatrix, spotPositions, spaceWeight=1, pcaDimensions=5){
  
  #check valid value for spaceWeight
  spaceWeight = as.double(spaceWeight)
  if(is.na(spaceWeight))
    stop("Param spaceWeight is not a double number.")
  if(spaceWeight > 1 | spaceWeight < 0)
    stop("spaceWeight is not in [0,1]")
  
  #check valid value for countMatrix and spotPosition
  if(!is.data.frame(countMatrix))
    stop("countMatrix is not a dataframe")
  if(!is.data.frame(spotPositions))
    stop("spotPositions is not a dataframe")
  
  countMatrix <- countMatrix[,sort(colnames(countMatrix))]
  d <- dim(spotPositions)[2]
  spotPositions <- spotPositions[,(d-1):d]
  spotPositions <- spotPositions[sort(rownames(spotPositions)),]

  seuratObject <- Seurat::CreateSeuratObject(countMatrix)
  seuratObject <- suppressWarnings(Seurat::SCTransform(seuratObject, assay = "RNA", verbose = FALSE))
  seuratObject <- Seurat::RunPCA(seuratObject, assay = "SCT", verbose = FALSE)
  if(pcaDimensions <= 2){
    pcaDimensions = 2
  }
  m <- seuratObject@reductions[["pca"]]@cell.embeddings[,1:pcaDimensions]
  distPCA = dist(m,method="minkowski",p=2)  
  distCoord <- dist(spotPositions,method="minkowski",p=2)
  distCoord <- distCoord*((max(distPCA)*as.double(spaceWeight))/(max(distCoord)))
  finalDistance <- as.matrix(distPCA + distCoord)
  neighbors <- suppressWarnings(Seurat::FindNeighbors(finalDistance))
  neighbors <- list(neighbors_nn=neighbors$nn,neighbors_snn=neighbors$snn)
  seuratObject@graphs <- neighbors
  seuratObject <- suppressWarnings(Seurat::FindClusters(seuratObject, verbose = FALSE, graph.name = "neighbors_snn"))
  seuratObject
}