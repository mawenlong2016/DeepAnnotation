
CompFunc <- function(FuncTerms=NULL, cpu=1, AmplitudePCA=NULL
                     , AmplitudeICA=NULL, AmplitudeNMF=NULL, rangeID=NULL){
  
  source("Ttest_MAP.R")
  source("EGP.R")
  AmplitudeMF <- cbind(AmplitudePCA,AmplitudeICA,AmplitudeNMF)
  geneName=rownames(AmplitudeMF)
  print("Gene functions complement begains!")
  positiveGene <- unique(FuncTerms[,1])
  geneSetTerm <- unique(FuncTerms[,2])
  
  PCAScore <- matrix(, length(geneName), length(rangeID))
  rownames(PCAScore) <- geneName
  colnames(PCAScore) <- geneSetTerm[rangeID]
  ICAScore <- NMFScore <- MFScore <- ENScore <- RandomScore <- PCAScore
  
  
  PCARank <- matrix(, length(geneName), length(rangeID))
  rownames(PCARank) <- geneName
  colnames(PCARank) <- geneSetTerm[rangeID]
  ICARank <- NMFRank <- MFRank <- ENRank <- RandomRank <- PCARank
  
  for(geneSetTermIdx in rangeID){#length(geneSetTerm)
    t1 <- Sys.time()
    overlapTarget <- unique(FuncTerms[which(FuncTerms[,2] == geneSetTerm[geneSetTermIdx]),1])
    # if(length(overlapTarget) > 2){
      print(paste0(geneSetTermIdx, " of ", length(geneSetTerm), " is ", geneSetTerm[geneSetTermIdx]," has genes ", length(overlapTarget)))
      getTrain_nega_all <- setdiff(positiveGene, overlapTarget)    
      targetTrain_nega_all <- setdiff(geneName, overlapTarget)      
 
      # LOOCV
      options(stringsAsFactors=F)
      options(scipen=999)
      suppressPackageStartupMessages(library(doParallel))
      suppressPackageStartupMessages(library(foreach))
      cl <-  makeCluster(cpu)
      registerDoParallel(cl)
      #do parallel computation
      tempList <- foreach(i = 1 : length(overlapTarget)) %dopar%{
        source("Ttest_MAP.R")
        source("EGP.R")
        targetTest <- overlapTarget[i]
        targetTrain <- overlapTarget[-i]
        targetTrain_nega <- setdiff(geneName, targetTrain)        
        
        # prediction based on student t-test and matrix factorization
        # PCA
        PCAMat <- Ttest_MAP(Amplitude = AmplitudePCA, targetTrain = targetTrain
                            , targetTrain_nega_all = targetTrain_nega_all
                            , targetTrain_nega = targetTrain_nega, targetTest=targetTest)
        if(length(which(PCAMat[,1] == targetTest)) > 0){
          PCAtemp <- cbind(targetTest, geneSetTerm[geneSetTermIdx], PCAMat[which(PCAMat[,1] == targetTest),2], which(PCAMat[,1] == targetTest), "label", "PCA")
          colnames(PCAtemp) <- c("geneName", "GO", "score", "rank", "annotation", "MF")
        }else{
          PCAtemp <- cbind(targetTest, geneSetTerm[geneSetTermIdx], NA, NA, "label", "PCA")
          colnames(PCAtemp) <- c("geneName", "GO", "score", "rank", "annotation", "MF")
        }
        # ICA
        ICAMat <- Ttest_MAP(Amplitude = AmplitudeICA, targetTrain = targetTrain
                            , targetTrain_nega_all = targetTrain_nega_all
                            , targetTrain_nega = targetTrain_nega, targetTest=targetTest)
        if(length(which(ICAMat[,1] == targetTest)) > 0){
          ICAtemp <- cbind(targetTest, geneSetTerm[geneSetTermIdx], ICAMat[which(ICAMat[,1] == targetTest),2], which(ICAMat[,1] == targetTest), "label", "ICA")
          colnames(ICAtemp) <- c("geneName", "GO", "score", "rank", "annotation", "MF")
        }else{
          ICAtemp <- cbind(targetTest, geneSetTerm[geneSetTermIdx], NA, NA, "label", "ICA")
          colnames(ICAtemp) <- c("geneName", "GO", "score", "rank", "annotation", "MF")
        }
        # NMF
        NMFMat <- Ttest_MAP(Amplitude = AmplitudeNMF, targetTrain = targetTrain
                            , targetTrain_nega_all = targetTrain_nega_all
                            , targetTrain_nega = targetTrain_nega, targetTest=targetTest)
        if(length(which(NMFMat[,1] == targetTest)) > 0){
          NMFtemp <- cbind(targetTest, geneSetTerm[geneSetTermIdx], NMFMat[which(NMFMat[,1] == targetTest),2], which(NMFMat[,1] == targetTest), "label", "NMF")
          colnames(NMFtemp) <- c("geneName", "GO", "score", "rank", "annotation", "MF")
        }else{
          NMFtemp <- cbind(targetTest, geneSetTerm[geneSetTermIdx], NA, NA, "label", "NMF")
          colnames(NMFtemp) <- c("geneName", "GO", "score", "rank", "annotation", "MF")
        }
        # MF
        MFMat <- Ttest_MAP(Amplitude = AmplitudeMF, targetTrain = targetTrain
                           , targetTrain_nega_all = targetTrain_nega_all
                           , targetTrain_nega = targetTrain_nega, targetTest=targetTest)
        if(length(which(MFMat[,1] == targetTest)) > 0){
          MFtemp <- cbind(targetTest, geneSetTerm[geneSetTermIdx], MFMat[which(MFMat[,1] == targetTest),2], which(MFMat[,1] == targetTest), "label", "MF")
          colnames(MFtemp) <- c("geneName", "GO", "score", "rank", "annotation", "MF")
        }else{
          MFtemp <- cbind(targetTest, geneSetTerm[geneSetTermIdx], NA, NA, "label", "MF")
          colnames(MFtemp) <- c("geneName", "GO", "score", "rank", "annotation", "MF")
        }
        # randomforest
        # time is too long
        
        # random selection
        # background
        Randomtemp <- cbind(targetTest, geneSetTerm[geneSetTermIdx], runif(1,0,1), sample(1:(length(geneName) - length(targetTrain)),1), "label", "Random")
        colnames(Randomtemp) <- c("geneName", "GO", "score", "rank", "annotation", "MF")        
        
        # ensemble
        ENMat <- EGP(Matlist = list(PCAMat[,1], ICAMat[,1], NMFMat[,1], MFMat[,1]))
        if(length(which(ENMat[,3] == targetTest)) > 0){
          ENtemp <- cbind(targetTest, geneSetTerm[geneSetTermIdx], ENMat[which(ENMat[,3] == targetTest),1], which(ENMat[,3] == targetTest), "label", "ensemble")
          colnames(ENtemp) <- c("geneName", "GO", "score", "rank", "annotation", "MF")
        }else{
          ENtemp <- cbind(targetTest, geneSetTerm[geneSetTermIdx], NA, NA, "label", "ensemble")
          colnames(ENtemp) <- c("geneName", "GO", "score", "rank", "annotation", "MF")
        }       
        return(list(PCAtemp = PCAtemp, ICAtemp = ICAtemp, NMFtemp = NMFtemp, 
                    MFtemp = MFtemp, ENtemp = ENtemp, 
                    Randomtemp = Randomtemp                    
        ))
      }
      names(tempList) <- c(1 : length(overlapTarget))
      stopCluster(cl)  
      
      for(pp in 1:length(tempList)){ 
        PCAScore[cbind(tempList[[pp]]$PCAtemp[,1], tempList[[pp]]$PCAtemp[,2])] <- as.numeric(tempList[[pp]]$PCAtemp[,3])
        ICAScore[cbind(tempList[[pp]]$ICAtemp[,1], tempList[[pp]]$ICAtemp[,2])] <- as.numeric(tempList[[pp]]$ICAtemp[,3])
        NMFScore[cbind(tempList[[pp]]$NMFtemp[,1], tempList[[pp]]$NMFtemp[,2])] <- as.numeric(tempList[[pp]]$NMFtemp[,3])
        MFScore[cbind(tempList[[pp]]$MFtemp[,1], tempList[[pp]]$MFtemp[,2])] <- as.numeric(tempList[[pp]]$MFtemp[,3])
        ENScore[cbind(tempList[[pp]]$ENtemp[,1], tempList[[pp]]$ENtemp[,2])] <- as.numeric(tempList[[pp]]$ENtemp[,3])
        RandomScore[cbind(tempList[[pp]]$Randomtemp[,1], tempList[[pp]]$Randomtemp[,2])] <- as.numeric(tempList[[pp]]$Randomtemp[,3])
        
        PCARank[cbind(tempList[[pp]]$PCAtemp[,1], tempList[[pp]]$PCAtemp[,2])] <- as.numeric(tempList[[pp]]$PCAtemp[,4])
        ICARank[cbind(tempList[[pp]]$ICAtemp[,1], tempList[[pp]]$ICAtemp[,2])] <- as.numeric(tempList[[pp]]$ICAtemp[,4])
        NMFRank[cbind(tempList[[pp]]$NMFtemp[,1], tempList[[pp]]$NMFtemp[,2])] <- as.numeric(tempList[[pp]]$NMFtemp[,4])
        MFRank[cbind(tempList[[pp]]$MFtemp[,1], tempList[[pp]]$MFtemp[,2])] <- as.numeric(tempList[[pp]]$MFtemp[,4])
        ENRank[cbind(tempList[[pp]]$ENtemp[,1], tempList[[pp]]$ENtemp[,2])] <- as.numeric(tempList[[pp]]$ENtemp[,4])
        RandomRank[cbind(tempList[[pp]]$Randomtemp[,1], tempList[[pp]]$Randomtemp[,2])] <- as.numeric(tempList[[pp]]$Randomtemp[,4])
      }
      
      #########
      # predict
      targetTrain <- overlapTarget
      targetTrain_nega <- setdiff(geneName, targetTrain)
      
      # PCA
      PCAMat <- Ttest_MAP(Amplitude = AmplitudePCA, targetTrain = targetTrain
                          , targetTrain_nega_all = targetTrain_nega_all
                          , targetTrain_nega = targetTrain_nega, targetTest=targetTest)
      PCAtemp <- cbind(rownames(PCAMat), geneSetTerm[geneSetTermIdx], PCAMat[,2], 1:nrow(PCAMat), "Unlabel", "PCA")
      colnames(PCAtemp) <- c("geneName", "GO", "score", "rank", "annotation", "MF")
      
      # ICA
      ICAMat <- Ttest_MAP(Amplitude = AmplitudeICA, targetTrain = targetTrain
                          , targetTrain_nega_all = targetTrain_nega_all
                          , targetTrain_nega = targetTrain_nega, targetTest=targetTest)
      ICAtemp <- cbind(rownames(ICAMat), geneSetTerm[geneSetTermIdx], ICAMat[,2], 1:nrow(ICAMat), "Unlabel", "ICA")
      colnames(ICAtemp) <- c("geneName", "GO", "score", "rank", "annotation", "MF")
      
      # NMF
      NMFMat <- Ttest_MAP(Amplitude = AmplitudeNMF, targetTrain = targetTrain
                          , targetTrain_nega_all = targetTrain_nega_all
                          , targetTrain_nega = targetTrain_nega, targetTest=targetTest)
      NMFtemp <- cbind(rownames(NMFMat), geneSetTerm[geneSetTermIdx], NMFMat[,2], 1:nrow(NMFMat), "Unlabel", "NMF")
      colnames(NMFtemp) <- c("geneName", "GO", "score", "rank", "annotation", "MF")
      
      # MF
      MFMat <- Ttest_MAP(Amplitude = AmplitudeMF, targetTrain = targetTrain
                         , targetTrain_nega_all = targetTrain_nega_all
                         , targetTrain_nega = targetTrain_nega, targetTest=targetTest)
      MFtemp <- cbind(rownames(MFMat), geneSetTerm[geneSetTermIdx], MFMat[,2], 1:nrow(MFMat), "Unlabel", "MF")
      colnames(MFtemp) <- c("geneName", "GO", "score", "rank", "annotation", "MF")
      
      # randomforest
      # time is too long
      
      # random selection
      # background     
      Randomtemp <- cbind(setdiff(geneName, overlapTarget), geneSetTerm[geneSetTermIdx], runif(length(geneName)-length(overlapTarget),0,1), sample(1:(length(geneName)-length(overlapTarget)),length(geneName)-length(overlapTarget)), "Unlabel", "Random")
      colnames(Randomtemp) <- c("geneName", "GO", "score", "rank", "annotation", "MF")
      
      # ensemble
      ENMat <- EGP(Matlist = list(PCAMat[,1], ICAMat[,1], NMFMat[,1], MFMat[,1]))
      ENtemp <- cbind(rownames(ENMat), geneSetTerm[geneSetTermIdx], ENMat[,1], 1:nrow(ENMat), "Unlabel", "ensemble")
      colnames(ENtemp) <- c("geneName", "GO", "score", "rank", "annotation", "MF")
      
      
      PCAScore[cbind(PCAtemp[,1], PCAtemp[,2])] <- as.numeric(PCAtemp[,3])
      ICAScore[cbind(ICAtemp[,1], ICAtemp[,2])] <- as.numeric(ICAtemp[,3])
      NMFScore[cbind(NMFtemp[,1], NMFtemp[,2])] <- as.numeric(NMFtemp[,3])
      MFScore[cbind(MFtemp[,1], MFtemp[,2])] <- as.numeric(MFtemp[,3])
      ENScore[cbind(ENtemp[,1], ENtemp[,2])] <- as.numeric(ENtemp[,3])
      RandomScore[cbind(Randomtemp[,1], Randomtemp[,2])] <- as.numeric(Randomtemp[,3])
      
      
      PCARank[cbind(PCAtemp[,1], PCAtemp[,2])] <- as.numeric(PCAtemp[,4])
      ICARank[cbind(ICAtemp[,1], ICAtemp[,2])] <- as.numeric(ICAtemp[,4])
      NMFRank[cbind(NMFtemp[,1], NMFtemp[,2])] <- as.numeric(NMFtemp[,4])
      MFRank[cbind(MFtemp[,1], MFtemp[,2])] <- as.numeric(MFtemp[,4])
      ENRank[cbind(ENtemp[,1], ENtemp[,2])] <- as.numeric(ENtemp[,4])
      RandomRank[cbind(Randomtemp[,1], Randomtemp[,2])] <- as.numeric(Randomtemp[,4])
      t2 <- Sys.time()
      print(t2 - t1)

    
  }
  return(list(PCAScore=PCAScore, ICAScore=ICAScore, NMFScore=NMFScore
              , MFScore=MFScore, ENScore=ENScore, RandomScore=RandomScore
              , PCARank=PCARank,ICARank=ICARank,NMFRank=NMFRank
              , MFRank=MFRank, ENRank=ENRank, RandomRank=RandomRank))
}