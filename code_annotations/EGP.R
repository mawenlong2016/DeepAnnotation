
EGP <- function(Matlist = NULL){
  # EGP
  if(sum(lengths(Matlist)) > 0){
    library(RobustRankAggreg)
    EGPMat <- aggregateRanks(glist = Matlist)
    EGPMat <- as.matrix(EGPMat)
    EGPMat <- cbind(EGPMat[,2], c(1:nrow(EGPMat)), EGPMat[,1])
    EGPMat <- as.data.frame(EGPMat, stringsAsFactors = F)
    EGPMat[,1] <- as.numeric(EGPMat[,1])
    EGPMat[,2] <- as.numeric(EGPMat[,2])
    EGPMat[,1] <- (-log10(EGPMat[,1]) - min(-log10(EGPMat[,1])))/(max(-log10(EGPMat[,1])) - min(-log10(EGPMat[,1])))
    rownames(EGPMat) <- EGPMat[,3]
  }else{
    EGPMat = NULL
  }
 
  return(EGPMat)
}
