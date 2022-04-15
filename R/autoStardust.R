#' @title Automatic Stardust 
#' @description This function calls Seurat clustering feature, passing a custom distance matrix that integrates both
#' transcriptional similarity and spatial information.
#' @param countMatrix, dataframe with spots id as columns and genes id as rows.
#' @param spotPositions, dataframe with x and y spot coordinates as columns and spots id as rows.
#' @param pcaDimensions, number of PCA dimension to include in the distance computation.
#' @param res, double resolution for Louvain algorithm.
#' 
#' @author Eva Viesi, eva [dot] viesi [at] univr [dot] it, University of Verona
#' @author Simone Avesani, simone [dot] avesani [at] univr [dot] it, University of Verona
#' @author Giovanni Motterle, giovanni [dot] motterle [at] univr [dot] it, University of Verona
#'
#' @return Seurat object with cluster identifications for each spot under active.ident slot.
#' @export
autoStardust <- function(countMatrix, 
                         spotPositions, 
                         pcaDimensions=10, 
                         res=0.8){
  
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
  distCoord <- distCoord*(max(distPCA)/max(distCoord))
  expr_norm <- (distPCA - min(distPCA)) / (max(distPCA) - min(distPCA)) 
  distCoord <- (distCoord)*(as.double(as.vector(expr_norm)))
  finalDistance <- as.matrix(distPCA + distCoord)
  neighbors <- suppressWarnings(Seurat::FindNeighbors(finalDistance))
  neighbors <- list(neighbors_nn=neighbors$nn,neighbors_snn=neighbors$snn)
  seuratObject@graphs <- neighbors
  seuratObject <- suppressWarnings(Seurat::FindClusters(seuratObject, verbose = FALSE, graph.name = "neighbors_snn"))
  seuratObject
}