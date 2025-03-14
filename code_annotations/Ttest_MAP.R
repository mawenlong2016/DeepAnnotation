Ttest_MAP <- function(Amplitude = NULL, targetTrain = NULL
                      , targetTrain_nega_all = NULL
                      , targetTrain_nega = NULL, targetTest=NULL){
  # Ttest
  targetTrain <- intersect(targetTrain, rownames(Amplitude))
  targetTrain_nega_all <- intersect(targetTrain_nega_all, rownames(Amplitude))
  targetTrain_nega <- intersect(targetTrain_nega, rownames(Amplitude))
  if(length(targetTrain) >= 2){
    targetPoolZscore <- matrix(NA ,ncol(Amplitude), 1)
    for (j in 1:ncol(Amplitude)) {
      dat1 <- Amplitude[targetTrain,j]
      dat2 <- Amplitude[targetTrain_nega_all,j]
      temp <- t.test(dat1, dat2)
      if(temp$p.value == 0){temp$p.value=1.0e-100}
      # targetPoolPCA[j,1] <- t.test(dat1, dat2)$statistic#T-test
      targetPoolZscore[j,1] <- qnorm(temp$p.value/2) * if(temp$statistic > 0) -1 else 1#T-test
    }
    calCor <- apply(Amplitude[targetTrain_nega, ], 1, function(ll){
      temp <- cor.test(targetPoolZscore[,1], ll)
      if(is.na(temp$estimate)){
        return(0)#return(NA)
      }else{
        if(temp$p.value == 0){temp$p.value=1.0e-100}
        return(qnorm(temp$p.value/2) * if(temp$estimate > 0) -1 else 1)
      }
    })
    #target test rank based on T-test
    tempT <- na.exclude(calCor)
    tempTOrder <- tempT[order(tempT, decreasing = TRUE)]
    TMat <- as.data.frame(cbind(names(tempTOrder), tempTOrder), stringsAsFactors = F)
    TMat[,1] <- as.character(TMat[,1])
    TMat[,2] <- as.numeric(TMat[,2])
    TMat[,2] <- (TMat[,2] - min(TMat[,2]))/(max(TMat[,2]) - min(TMat[,2]))
    rownames(TMat) <- TMat[,1]
  }else if(length(targetTrain) == 1){
    targetPoolZscore <- Amplitude[targetTrain,]
    targetPoolZscore <- as.matrix(targetPoolZscore)
    calCor <- apply(Amplitude[targetTrain_nega, ], 1, function(ll){
      temp <- cor.test(targetPoolZscore[,1], ll)
      if(is.na(temp$estimate)){
        return(0)#return(NA)
      }else{
        if(temp$p.value == 0){temp$p.value=1.0e-100}
        return(qnorm(temp$p.value/2) * if(temp$estimate > 0) -1 else 1)
      }
    })
    #target test rank based on T-test
    tempT <- na.exclude(calCor)
    tempTOrder <- tempT[order(tempT, decreasing = TRUE)]
    TMat <- as.data.frame(cbind(names(tempTOrder), tempTOrder), stringsAsFactors = F)
    TMat[,1] <- as.character(TMat[,1])
    TMat[,2] <- as.numeric(TMat[,2])
    TMat[,2] <- (TMat[,2] - min(TMat[,2]))/(max(TMat[,2]) - min(TMat[,2]))
    rownames(TMat) <- TMat[,1]
  }else{
    targetTrain=targetTest
    targetPoolZscore <- Amplitude[targetTrain,]
    targetPoolZscore <- as.matrix(targetPoolZscore)
    calCor <- apply(Amplitude[targetTrain_nega, ], 1, function(ll){
      temp <- cor.test(targetPoolZscore[,1], ll)
      if(is.na(temp$estimate)){
        return(0)#return(NA)
        }else{
          if(temp$p.value == 0){temp$p.value=1.0e-100}
          return(qnorm(temp$p.value/2) * if(temp$estimate > 0) -1 else 1)
        }
        })
    #target test rank based on T-test
    tempT <- na.exclude(calCor)
    tempTOrder <- tempT[order(tempT, decreasing = TRUE)]
    TMat <- as.data.frame(cbind(names(tempTOrder), tempTOrder), stringsAsFactors = F)
    TMat[,1] <- as.character(TMat[,1])
    TMat[,2] <- as.numeric(TMat[,2])
    TMat[,2] <- (TMat[,2] - min(TMat[,2]))/(max(TMat[,2]) - min(TMat[,2]))
    rownames(TMat) <- TMat[,1]
  }
  # else{TMat <- NULL}
  
  return(TMat)
}
