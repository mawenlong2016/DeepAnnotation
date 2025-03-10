# load the prediction results
# the original Figure2.RData could be downloaded via: https://doi.org/10.5281/zenodo.14586911
# setwd('/data3/mwl/DeepAnnotation_review/gigadb')
# topPCC : denotes the top-ranked samples' PCC
# meanNDCG: denotes the top-ranked samples' meanN DCG
# meanRE: denotes the top-ranked samples' mean RE
load('./figure_data/Figure2.RData')
if(!file.exists('./figure_plot')){
		system(paste0("mkdir ./figure_plot"), intern =TRUE)
}
pred_phenotype_cal <- cbind(pred_phenotype, pred_phenotype_lightgbm[,2],pred_phenotype_KAML[,2],pred_phenotype_BLUP[,2],pred_phenotype_BayesR[,2],pred_phenotype_MBLUP[,2],pred_phenotype_BayesRC[,2], Genotype[,2], SNPAnnotation[,2], GeneFunction[,2], Network[,2])
colnames(pred_phenotype_cal) <- c("real", "rrBLUP", "LightGBM", "KAML", "BLUP", "BayesR", "MBLUP", "BayesRC","Genotype", "SNPAnnotation", "GeneFunction", "Network")
# top-ranked PCC calculation
cor_spatial_result <- NULL
for(ratio in 1:100){
	a1 = order(pred_phenotype_cal[,1], decreasing=T)[1:max(ceiling(nrow(pred_phenotype_cal)*ratio/100),5)]
	cor_spatial_result <- rbind(cor_spatial_result, cor(pred_phenotype_cal[a1,])[1,])
}
rownames(cor_spatial_result) <- paste0("top ",1:100,"%")
topPCC <- cor_spatial_result[1:10,]
# NDCG calculation
##from plos one, 2015, A Ranking Approach to Genomic Selection
NDCG <- function( realScores, predScores, topK = 10){  
  if( length(realScores) != length(predScores)){
    stop("Error: different length between realScores and predScores")
  }
  if( length(realScores) < topK ) {
    stop("Error: too large topK")
  }  
  scoreMat <- cbind(realScores,predScores)
  scoreMatSortbyPred <- scoreMat[order(scoreMat[,2],decreasing = TRUE),]
  scoreMatSortByReal <- scoreMat[order(scoreMat[,1],decreasing = TRUE),]  
  DCG <- rep(0, topK)
  IDCG <- rep(0, topK)
  for(idx in 1:topK){
    DCG[idx] <-  scoreMatSortbyPred[idx,1]/log(idx+1,2)
    IDCG[idx] <- scoreMatSortByReal[idx,1]/log(idx+1,2)
  }  
  NDCG <- sum(DCG)/sum(IDCG) 
  names(NDCG) <- paste0("NDCG","_top",topK)
  return(NDCG)
}
NDCGEvaluation <- function(method = NULL,topAlpha,realScores,predScores, allNDCG = T){  
  if(allNDCG == T){
    ndcg <- 1 : round(length(realScores)*(max(topAlpha)/100))
  }else{
    ndcg <- round(length(realScores)*(topAlpha/100))
  }  
  if(method == "NDCG"){
    result <- sapply(ndcg,function(ii){
      NDCG(realScores = realScores,predScores = predScores,topK = ii)      
    })
    if(allNDCG == F){
      names(result) <- paste0("NDCG_top",topAlpha)
    }
  }else if(method == "meanNDCG"){
    ndcg <- 1 : round(length(realScores)*(max(topAlpha)/100))
    NDCGEval <- sapply(ndcg,function(ii){
      NDCG(realScores = realScores,predScores = predScores,topK = ii)
    })
    result <- sapply(topAlpha,function(ii){
      iii <- round(ii*length(realScores)/100)
      sum(NDCGEval[1:iii])/iii
    })
    names(result) <- paste0("meanNDCG_top",topAlpha)
  }
  result
}
meanNDCG = apply(pred_phenotype_cal[,2:ncol(pred_phenotype_cal)],2,function(x){NDCGEvaluation(method = "meanNDCG",topAlpha = 1:10,realScores = pred_phenotype_cal[,1],predScores = x,allNDCG = F)})
# relative efficiency (RE) calculation
##from heredity, 2014, Genomic-enabled prediction with classification algorithms
##from plant genome, 2018, Applications of Machine Learning Methods to Genomic Selection in Breeding Wheat for Rust Resistance
RE <- function(realScores, predScores, topK = 10){  
  if( length(realScores) != length(predScores)){
    stop("Error: different length between realScores and predScores")
  }
  if( length(realScores) < topK ) {
    stop("Error: too large topK")
  }  
  scoreMat <- cbind(realScores,predScores)
  scoreMatSortbyPred <- scoreMat[order(scoreMat[,2],decreasing = TRUE),]
  scoreMatSortByReal <- scoreMat[order(scoreMat[,1],decreasing = TRUE),]  
  RE <- (mean(scoreMatSortbyPred[1:topK,1]) - mean(scoreMatSortByReal[,1]))/(mean(scoreMatSortByReal[1:topK,1]) - mean(scoreMatSortByReal[,1]))
  names(RE) <- paste0("RE","_top",topK)
  return(RE)
}
REevaluation <- function(method = NULL, topAlpha,realScores,predScores, allRE = T){  
  if(allRE == T){
    re <- 1 : round(length(realScores)*(max(topAlpha)/100))
  }else{
    re <- round(length(realScores)*(topAlpha/100))
  }
  if(method == "RE"){
    result <- sapply(re,function(ii){
      RE(realScores = realScores,predScores = predScores,topK = ii)      
    })
    if(allRE == F){
      names(result) <- paste0("RE_top",topAlpha)
    } 
  }else if(method == "meanRE"){
    re <- 1 : round(length(realScores)*(max(topAlpha)/100))
    REEval <- sapply(re,function(ii){
      RE(realScores = realScores,predScores = predScores,topK = ii)
    })
    result <- sapply(topAlpha,function(ii){
      iii <- round(ii*length(realScores)/100)
      sum(REEval[1:iii])/iii
    })
    names(result) <- paste0("meanRE_top",topAlpha)
  }   
  result
}
meanRE = apply(pred_phenotype_cal[,2:ncol(pred_phenotype_cal)],2,function(x){REevaluation(method = 'meanRE', topAlpha = 1:10,realScores = pred_phenotype_cal[,1],predScores = x,allRE = F)})
# pheatmap plot
pdf("./figure_plot/LMP_Figure2C_top_rank_evaluation.pdf",height=16, width=12)
library(pheatmap)
pheatmap(topPCC[1:10,c(2,3,4,5,6,9,7,8,10,11,12)], cluster_rows = F, cluster_cols = F, scale='none', border_color = NA)
pheatmap(meanNDCG[1:10,c(1,2,3,4,5,8,6,7,9,10,11)], cluster_rows = F, cluster_cols = F, scale='none', border_color = NA)
pheatmap(meanRE[1:10,c(1,2,3,4,5,8,6,7,9,10,11)], cluster_rows = F, cluster_cols = F, scale='none', border_color = NA)
dev.off()
cor_spatial_pvalue_all <- NULL
for(l in 9:12){
	cor_spatial_pvalue <- NULL
	for(k in 2:12){
		cor_spatial_pvalue <- c(cor_spatial_pvalue, t.test(cor_spatial_result[1:10,l], cor_spatial_result[1:10,k], alternative = 'greater', paired = T)$p.value)
	}
	cor_spatial_pvalue_all <- rbind(cor_spatial_pvalue_all, cor_spatial_pvalue)
}
colnames(cor_spatial_pvalue_all) <- c("rrBLUP", "LightGBM", "KAML", "BLUP", "BayesR", "MBLUP", "BayesRC","Genotype", "SNPAnnotation", "GeneFunction", "Network")
meanNDCG_pvalue_all <- NULL
for(l in 8:11){
	meanNDCG_pvalue <- NULL
	for(k in 1:11){
		meanNDCG_pvalue <- c(meanNDCG_pvalue, t.test(meanNDCG[1:10,l], meanNDCG[1:10,k], alternative = 'greater', paired = T)$p.value)
	}
	meanNDCG_pvalue_all <- rbind(meanNDCG_pvalue_all, meanNDCG_pvalue)
}
colnames(meanNDCG_pvalue_all) <- c("rrBLUP", "LightGBM", "KAML", "BLUP", "BayesR", "MBLUP", "BayesRC","Genotype", "SNPAnnotation", "GeneFunction", "Network")
meanRE_pvalue_all <- NULL
for(l in 8:11){
	meanRE_pvalue <- NULL
	for(k in 1:11){
		meanRE_pvalue <- c(meanRE_pvalue, t.test(meanRE[1:10,l], meanRE[1:10,k], alternative = 'greater', paired = T)$p.value)
	}
	meanRE_pvalue_all <- rbind(meanRE_pvalue_all, meanRE_pvalue)
}
colnames(meanRE_pvalue_all) <- c("rrBLUP", "LightGBM", "KAML", "BLUP", "BayesR", "MBLUP", "BayesRC","Genotype", "SNPAnnotation", "GeneFunction", "Network")