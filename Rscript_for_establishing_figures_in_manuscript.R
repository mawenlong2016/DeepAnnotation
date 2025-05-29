
# Rscript for establishing figures in manuscript.

# ******************** Main Figures******************************** #
########################----Figure 2----########################
# Figure 2. Prediction performance of DeepAnnotation through cross-validation compared with rrBLUP, LightGBM, KAML, BLUP, BayesR, MultiBLUP, BayesRC, and DeepAnnotation on LMP trait.
#--------------------------------------------------------------
# setwd('/data3/mwl/DeepAnnotation_review/gigadb')
load('./figure_data/Figure2.RData')
if(!file.exists('./figure_plot')){
		system(paste0("mkdir ./figure_plot"), intern =TRUE)
}
#--------------------------------------------------------------
# Figure 2A, Pearson correlation coefficient scores of different hyperparameters combinations based on a 2-fold cross-validation experiment from DeepAnnotation using different types of functional annotation data. 
options(scipen=100)
learning_rate <- c(0.1, 0.01, 0.001, 0.0001, 0.00001)
unit <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130)
regularizer <- c("l1", "l2", "l1_l2","none")
regularizer_rate <- c(0.1, 0.01, 0.001, 0.0001)
momentum <- c(0.99, 0.98, 0.95, 0.9, 0.8)
dropout <- c(0.3, 0.2, 0.1)
unit_ratio <- c("1,1,5,5,1,1", "1,2,3,4,5,6", "1,2,4,5,6,7", "1,3,2,2,3,1", 
	            "1,2,2,5,4,3", "2,2,1,1,2,2")
para_flag <- cbind(c(1,125),c(126,250),c(251,375),c(375,500))
vivo_list_plot <- list()
omic_name <- c("genome","snp","function","network")
for(kk in 1:ncol(para_flag)){
	vivo_list <- list()
	value_result_mean <- NULL
	value_result_median <- NULL
	opt_result <- opt_result_back[para_flag[1,kk]:para_flag[2,kk],]
	for(i in 1:length(learning_rate)){
		vivo_list[[paste0("LR_",learning_rate[i])]] <- as.numeric(opt_result[which(opt_result[,14] == as.vector(learning_rate[i])),3])
		value_result_mean <- c(value_result_mean, mean(vivo_list[[paste0("LR_",learning_rate[i])]]))
		value_result_median <- c(value_result_median, median(vivo_list[[paste0("LR_",learning_rate[i])]]))
	}
	for(i in 1:length(unit)){
		vivo_list[[paste0("Unit_",unit[i])]] <- as.numeric(opt_result[which(opt_result[,13] == as.vector(unit[i])),3])
		value_result_mean <- c(value_result_mean, mean(vivo_list[[paste0("Unit_",unit[i])]]))
		value_result_median <- c(value_result_median, median(vivo_list[[paste0("Unit_",unit[i])]]))
	}
	for(i in 1:length(regularizer)){
		vivo_list[[paste0("Regu_",regularizer[i])]] <- as.numeric(opt_result[which(opt_result[,12] == as.vector(regularizer[i])),3])
		value_result_mean <- c(value_result_mean, mean(vivo_list[[paste0("Regu_",regularizer[i])]]))
		value_result_median <- c(value_result_median, median(vivo_list[[paste0("Regu_",regularizer[i])]]))
	}
	for(i in 1:length(regularizer_rate)){
		vivo_list[[paste0("ReRa_",regularizer_rate[i])]] <- as.numeric(opt_result[which(opt_result[,11] == as.vector(regularizer_rate[i])),3])
		value_result_mean <- c(value_result_mean, mean(vivo_list[[paste0("ReRa_",regularizer_rate[i])]]))
		value_result_median <- c(value_result_median, median(vivo_list[[paste0("ReRa_",regularizer_rate[i])]]))
	}
	for(i in 1:length(momentum)){
		vivo_list[[paste0("MOME_",momentum[i])]] <- as.numeric(opt_result[which(opt_result[,10] == as.vector(momentum[i])),3])
		value_result_mean <- c(value_result_mean, mean(vivo_list[[paste0("MOME_",momentum[i])]]))
		value_result_median <- c(value_result_median, median(vivo_list[[paste0("MOME_",momentum[i])]]))
	}
	for(i in 1:length(dropout)){
		vivo_list[[paste0("Drop_",dropout[i])]] <- as.numeric(opt_result[which(opt_result[,9] == as.vector(dropout[i])),3])
		value_result_mean <- c(value_result_mean, mean(vivo_list[[paste0("Drop_",dropout[i])]]))
		value_result_median <- c(value_result_median, median(vivo_list[[paste0("Drop_",dropout[i])]]))
	}
	for(i in 1:length(unit_ratio)){
		vivo_list[[paste0("UR_",unit_ratio[i])]] <- as.numeric(opt_result[which(opt_result[,8] == as.vector(unit_ratio[i])),3])
		value_result_mean <- c(value_result_mean, mean(vivo_list[[paste0("UR_",unit_ratio[i])]]))
		value_result_median <- c(value_result_median, median(vivo_list[[paste0("UR_",unit_ratio[i])]]))
	}
	vivo_list_plot[[paste0(omic_name[kk])]] <- list(vivo_list = vivo_list,
													value_result_mean = value_result_mean,
													value_result_median = value_result_median)
}
write.table(cbind(vivo_list_plot[[1]]$value_result_median, 
	vivo_list_plot[[2]]$value_result_median,
	vivo_list_plot[[3]]$value_result_median,
	vivo_list_plot[[4]]$value_result_median,
	vivo_list_plot[[1]]$value_result_mean, 
	vivo_list_plot[[2]]$value_result_mean, 
	vivo_list_plot[[3]]$value_result_mean, 
	vivo_list_plot[[4]]$value_result_mean), file = "./figure_plot/LMP_Figure2A_opt_genome.txt"
	, quote = F, row.names=F, col.names=F, sep="\t")
library(RColorBrewer)
# install.packages('vioplot')
library(vioplot)
col_break = brewer.pal(9, "Set3")
pdf("./figure_plot/LMP_Figure2A_opt_genome.pdf",height=16, width=12)
par(mfcol = c(4,1))
for(kk in 1:length(vivo_list_plot)){
	vivo_list <- vivo_list_plot[[kk]]$vivo_list
	value_result_mean <- vivo_list_plot[[kk]]$value_result_mean
	value_result_median <- vivo_list_plot[[kk]]$value_result_median
	vioplot(vivo_list, pch = 20, 
		border = col_break[9], col = col_break[9], 
		cex = 0.5,ylab = "PCC", main = names(vivo_list_plot)[kk],
		xaxt = "n", yaxt = "n")
	points(value_result_mean, pch = 18, col = 'black', lwd = 6)
	axis(side = 2)
	axis(side = 1, labels = names(vivo_list)
		, las = 2
		, at = 1:length(vivo_list))
}
dev.off()
#--------------------------------------------------------------
# Figure 2B, PCC evaluation of different models based on 5-fold cross-validation (CV) experiment. In the "CV together" approach, each testing fold from the 5-fold cross-validation were consolidated before metric computation. On the other hand, in the "CV separate" approach, metrics were computed separately for each testing fraction of the 5-fold cross-validation. The paired t-test P-values of DeepAnnotation compared with other models were displayed on the black boxes with coral represents significance with P-value < 0.05.
pred_phenotype_cal <- cbind(pred_phenotype, pred_phenotype_lightgbm[,2],pred_phenotype_KAML[,2],pred_phenotype_BLUP[,2],pred_phenotype_BayesR[,2],pred_phenotype_MBLUP[,2],pred_phenotype_BayesRC[,2], Genotype[,2], SNPAnnotation[,2], GeneFunction[,2], Network[,2])
colnames(pred_phenotype_cal) <- c("real", "rrBLUP", "LightGBM", "KAML", "BLUP", "BayesR", "MBLUP", "BayesRC","Genotype", "SNPAnnotation", "GeneFunction", "Network")
cor_cv_result <- NULL
for(cv in 1:nrow(Test)){
	cor_cv_result <- rbind(cor_cv_result, cor(pred_phenotype_cal[cv_flag[[paste0("cv",cv)]][["test"]],])[1,])
}
cor_cv_pvalue_all <- NULL
for(l in 9:12){
	cor_cv_pvalue <- NULL
	for(k in 2:12){
		cor_cv_pvalue <- c(cor_cv_pvalue, t.test(cor_cv_result[,l], cor_cv_result[,k], paired = T, alternative = 'greater')$p.value)
	}
	cor_cv_pvalue_all <- rbind(cor_cv_pvalue_all, cor_cv_pvalue)
}
colnames(cor_cv_pvalue_all) <- c("rrBLUP", "LightGBM", "KAML", "BLUP", "BayesR", "MBLUP", "BayesRC","Genotype", "SNPAnnotation", "GeneFunction", "Network")
cor_cv_pvalue_all
cor(pred_phenotype_cal)[,1]
model_name <- c("rrBLUP", "LightGBM", "KAML", "BLUP", "BayesR", "MBLUP", "BayesRC","Genotype", "SNPAnnotation", "GeneFunction", "Network")
vivo_list_sparial_pcc <- list()
vivo_list_sparial_pcc_mean <- NULL
for(k in 1:length(model_name)){
	vivo_list_sparial_pcc[[model_name[k]]] <- as.numeric(cor_cv_result[,k+1])
	vivo_list_sparial_pcc_mean <- c(vivo_list_sparial_pcc_mean, mean(cor_cv_result[,k+1]))
}
library(RColorBrewer)
col_break = brewer.pal(12, "Set3")
col_break = col_break[c(9, 5, 6, 11, 12, 1, 2, 7, 8, 10, 4)]
library(vioplot)
pdf("./figure_plot/LMP_Figure2B_PCC_cv.pdf",height=16, width=12)
vioplot(vivo_list_sparial_pcc, pch = 20, 
	border = col_break, col = col_break, 
	cex = 0.5,ylab = "PCC", 
	main = "PCC", 
	xaxt = "n", yaxt = "n")
axis(side = 2)
axis(side = 1
	, las = 2
	, at = 1:length(vivo_list_sparial_pcc))
points(vivo_list_sparial_pcc_mean, pch = 18, col = 'black', lwd = 6)
dev.off()
pdf("./figure_plot/LMP_Figure2B_PCC.pdf",height=16, width=12)
barplot(cor(pred_phenotype_cal)[,1][2:12])
dev.off()
#--------------------------------------------------------------
# Figure 2C, Relative efficiency (RE) values between the predicted and observed phenotypic values of top-ranked samples from top 1% to top 10%. The paired t-test P-values were displayed on the black boxes with coral represents significance with P-value < 0.05. 
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
      names(result) <- paste0("top",topAlpha,"%")
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
    names(result) <- paste0("top",topAlpha,"%")
  }   
  result
}
RE_calculate = apply(pred_phenotype_cal[,2:ncol(pred_phenotype_cal)],2,function(x){REevaluation(method = 'RE', topAlpha = 1:10,realScores = pred_phenotype_cal[,1],predScores = x,allRE = F)})
pdf("./figure_plot/LMP_Figure2C_PCC_pheat.pdf",height=16, width=12)
library(pheatmap)
pheatmap(RE_calculate[1:10,c(1,2,3,4,5,8,6,7,9,10,11)], cluster_rows = F, cluster_cols = F, scale='none', border_color = NA)
dev.off()
RE_pvalue_all <- NULL
for(l in 8:11){
	RE_pvalue <- NULL
	for(k in 1:11){
		RE_pvalue <- c(RE_pvalue, t.test(RE_calculate[1:10,l], RE_calculate[1:10,k], alternative = 'greater', paired = T)$p.value)
	}
	RE_pvalue_all <- rbind(RE_pvalue_all, RE_pvalue)
}
colnames(RE_pvalue_all) <- c("rrBLUP", "LightGBM", "KAML", "BLUP", "BayesR", "MBLUP", "BayesRC","Genotype", "SNPAnnotation", "GeneFunction", "Network")
RE_pvalue_all
#--------------------------------------------------------------
# Figure 2D, Elapsed training time of different models for 5 CVs. Each column of vertical points represents one CV. For DeepAnnotation, the training times are not including the hyperparameter optimization. 
library(RColorBrewer)
col_break = brewer.pal(12, "Set3")
col_break = col_break[c(9, 5, 6, 11, 12, 1, 2, 7, 8, 10, 4)]
# model_running_time <- read.table("./figure_data/model_running_time.txt", header = T)
pdf("./figure_plot/LMP_Figure2D_time.pdf",height=6, width=6)
plot(model_running_time[1:5,1], pch = 20, type='b', col=col_break[1], ylim = c(0,3.2))
for(i in 2:ncol(model_running_time))
lines(model_running_time[1:5,i], pch = 20, type = 'b', col=col_break[i])
dev.off()
time_pvalue_all <- NULL
for(k in 8:11){
	time_pvalue <- NULL
	for(l in 1:11){
		time_pvalue <- c(time_pvalue, t.test(model_running_time[1:5,k], model_running_time[1:5,l], alternative = 'less', paired = T)$p.value)
	}
	time_pvalue_all <- rbind(time_pvalue_all, time_pvalue)
}
colnames(time_pvalue_all) <- c("rrBLUP", "LightGBM", "KAML", "BLUP", "BayesR", "MBLUP", "BayesRC","Genotype", "SNPAnnotation", "GeneFunction", "Network")
time_pvalue_all
# save(list = ls(), file = "./figure_data/Figure2.RData")
########################----Figure 2----########################
########################----Figure 3----########################
# Figure 3. Prediction performance of DeepAnnotation for independent test dataset of LMP.
#--------------------------------------------------------------
# setwd('/data3/mwl/DeepAnnotation_review/gigadb')
load('./figure_data/Figure3.RData')
if(!file.exists('./figure_plot')){
		system(paste0("mkdir ./figure_plot"), intern =TRUE)
}
#--------------------------------------------------------------
# Figure 3A,Training and validation status of DeepAnnotation for all functional annotations of each epoch based on the 5-fold cross-validation. 
# overall training steps
options(scipen=999)
Network_step_training <- do.call(rbind, strsplit(Network_step[grep('Predict training Step', Network_step[,1]),1], "( )|:"))
Network_step_star <- which(as.numeric(Network_step_training[,4]) == 1)
Network_step_end <- Network_step_star-1
Network_step_end <- Network_step_end[-1]
Network_step_end <- c(Network_step_end, nrow(Network_step_training))
Network_step_flag <- rbind(Network_step_star, Network_step_end)
Network_step_range <- Network_step_flag[2,] - Network_step_flag[1,] + 1
Network_step_range
# ************************************ #
Network_train_loss <- NULL
for(i in 1:5){
	Network_train_loss <- cbind(Network_train_loss,as.numeric(Network_step_training[Network_step_star[i]:(Network_step_star[i]+min(Network_step_range)-1),7]))
}
Network_train_loss <- cbind(Network_train_loss, apply(Network_train_loss, 1, mean))
# ************************************ #
Network_step_validing <- do.call(rbind, strsplit(Network_step[grep('Predict validation Step', Network_step[,1]),1], "( )|:"))
Network_step_star <- which(as.numeric(Network_step_validing[,4]) == 1)
Network_step_end <- Network_step_star-1
Network_step_end <- Network_step_end[-1]
Network_step_end <- c(Network_step_end, nrow(Network_step_validing))
Network_step_flag <- rbind(Network_step_star, Network_step_end)
Network_step_range <- Network_step_flag[2,] - Network_step_flag[1,] + 1
Network_step_range
# ************************************ #
Network_valid_loss <- NULL
for(i in 1:5){
	Network_valid_loss <- cbind(Network_valid_loss,as.numeric(Network_step_validing[Network_step_star[i]:(Network_step_star[i]+min(Network_step_range)-1),7]))
}
Network_valid_loss <- cbind(Network_valid_loss, apply(Network_valid_loss, 1, mean))
# ************************************ #
train_min_flag <- NULL
train_loss_bk <- 99999999
temp_num <- 0
for(i in 1:nrow(Network_train_loss)){
	temp <- NULL
	if(Network_train_loss[i,6] < train_loss_bk){
		train_loss_bk <- Network_train_loss[i,6]
		train_min_flag <- rbind(train_min_flag, cbind(i, Network_train_loss[i,6], temp_num))
		temp_num <- 0	
	}else{
		temp_num <- temp_num + 1
	}
}
which.max(train_min_flag[,3])
train_min_flag[(which.max(train_min_flag[,3])-1):(which.max(train_min_flag[,3])+1),]
valid_min_flag <- NULL
valid_loss_bk <- 99999999
temp_num <- 0
for(i in 1:nrow(Network_valid_loss)){
	temp <- NULL
	if(Network_valid_loss[i,6] < valid_loss_bk){
		valid_loss_bk <- Network_valid_loss[i,6]
		valid_min_flag <- rbind(valid_min_flag, cbind(i, Network_valid_loss[i,6], temp_num))
		temp_num <- 0	
	}else{
		temp_num <- temp_num + 1
	}
}
which.max(valid_min_flag[,3])
valid_min_flag[(which.max(valid_min_flag[,3])-1):(which.max(valid_min_flag[,3])+1),]
pdf("./figure_plot/LMP_Figure3A_Network_train_barplot.pdf"
	,width = 12, height = 4)
a1 <- train_min_flag[which(train_min_flag[,3] > 0),3]
barplot(a1, las = 2)
dev.off()
pdf("./figure_plot/LMP_Figure3A_Network_valid_barplot.pdf"
	,width = 12, height = 4)
a2 <- valid_min_flag[which(valid_min_flag[,3] > 0),3]
barplot(a2, las = 2)
dev.off()
library(RColorBrewer)
col_loss <- brewer.pal(8, "Set2")[c(8,7)]
pdf("./figure_plot/LMP_Figure3A_Network_loss_1_100.pdf"
	)
plot(Network_train_loss[1:100,6], type = "l", lwd = 2, col = col_loss[1],ylim=c(0,5000))
lines(Network_valid_loss[1:100,6], lwd = 2, col = col_loss[2])
abline(v = which.min(Network_train_loss[1:100,6]), col = "black")
abline(v = which.min(Network_valid_loss[1:100,6]), col = "black")
dev.off()

pdf("./figure_plot/LMP_Figure3_Network_loss_100_1008.pdf"
	)
plot(Network_train_loss[100:nrow(Network_train_loss),6], type = "l", lwd = 2, col = col_loss[1],ylim=c(5,20))
lines(Network_valid_loss[100:nrow(Network_train_loss),6], lwd = 2, col = col_loss[2])
abline(v = which.min(Network_train_loss[100:nrow(Network_train_loss),6]), col = "black")
abline(v = which.min(Network_valid_loss[100:nrow(Network_train_loss),6]), col = "black")
dev.off()
which.min(Network_train_loss[,6]) # 994
which.min(Network_valid_loss[,6]) # 973
#--------------------------------------------------------------
# Figure 3B,The distribution of RE scores of top 1~20 ranked samples under the distribution of 5 CV results.
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
      names(result) <- paste0("top",topAlpha,"%")
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
    names(result) <- paste0("top",topAlpha,"%")
  }   
  result
}

RE_calculate_cv <- NULL
for(cv in 1:5){
	pred_phenotype_cal <- cbind(pred_phenotype_test[,1],pred_phenotype_test[,cv+1],pred_phenotype_test_lightgbm[,cv+1],pred_phenotype_test_KAML[,cv+1],pred_phenotype_test_BLUP[,cv+1],pred_phenotype_test_BayesR[,cv+1],pred_phenotype_test_MBLUP[,cv+1],pred_phenotype_test_BayesRC[,cv+1], pred_phenotype_Network[,cv])
	colnames(pred_phenotype_cal) <- c("real", "rrBLUP", "LightGBM", "KAML", "BLUP", "BayesR", "MBLUP", "BayesRC","Network")
	RE_calculate = apply(pred_phenotype_cal[,2:ncol(pred_phenotype_cal)],2,function(x){REevaluation(method = 'RE', topAlpha = 1:100,realScores = pred_phenotype_cal[,1],predScores = x,allRE = T)})
	RE_calculate_cv[[cv]] <- RE_calculate
}
RE_calculate <- (RE_calculate_cv[[1]]+RE_calculate_cv[[2]]+RE_calculate_cv[[3]]+RE_calculate_cv[[4]]+RE_calculate_cv[[5]])/5

RE_calculate_pvalue <- NULL
for(k in 1:7){
  RE_calculate_pvalue <- c(RE_calculate_pvalue, t.test(RE_calculate[1:20,"Network"], RE_calculate[1:20,k], alternative = 'greater', paired=T)$p.value)
}
RE_calculate_pvalue

pdf("./figure_plot/LMP_Figure3B_RE.pdf",height=16, width=12)
vivo_list <- list()
for(i in 1:ncol(RE_calculate)){
	vivo_list[[i]] <- as.numeric(RE_calculate[1:20,i])
}
names(vivo_list) <- c("rrBLUP", "LightGBM", "KAML", "BLUP", "BayesR", "MBLUP", "BayesRC","Network")
value_result_mean <- apply(RE_calculate[1:20,1:8],2,mean)
library(vioplot)
library(RColorBrewer)
col_break = brewer.pal(12, "Set3")
col_break = col_break[c(9, 5, 6, 11, 12, 1, 2, 4)]
vioplot(vivo_list, pch = 20, 
	border = col_break, col = col_break, 
	cex = 0.5,ylab = "PCC",
	xaxt = "n", yaxt = "n")
points(value_result_mean, pch = 18, col = 'black', lwd = 6)
axis(side = 2)
axis(side = 1, labels = names(vivo_list)
	, las = 2
	, at = 1:length(vivo_list))
dev.off()
# save(list = ls(), file = "./figure_data/Figure3.RData")
########################----Figure 3----########################
########################----Figure 4----########################
# Figure 4. Interpretability of DeepAnnotation allows for the fine-mapping of potential causal SNPs. 
#--------------------------------------------------------------
# setwd('/data3/mwl/DeepAnnotation_review/gigadb')
load('./figure_data/Figure4.RData')
if(!file.exists('./figure_plot')){
		system(paste0("mkdir ./figure_plot"), intern =TRUE)
}
options(scipen=999)
#--------------------------------------------------------------
library(GO.db)
# sig_term_annotate <- cbind(sig_term[,3], sig_term[,1], Term(sig_term[,3]), Ontology(sig_term[,3]))
# annotate_sig_gene <- gold_gene[FuncTerms[which(FuncTerms[,2] == sig_term_annotate[which(as.numeric(sig_term_annotate[,2]) < 0.05),1]),1],]
annotate_sig_gene[which(is.na(annotate_sig_gene[,1]) == F),]
#--------------------------------------------------------------
# Figure 4B, Metaterm 8 acts as a co-regulatory module that regulates skeletal muscle development through epigenetic regulation.
sig_term_annotate[which(as.numeric(sig_term_annotate[,2]) < 0.05),]
library(pheatmap)
pdf(file="./figure_plot/Figure4B_LMP_metagene8_pheatmap.pdf")
pheatmap(-log10(metagene[,1]), cluster_rows=F, cluster_cols=F, cutree_rows=F, cutree_cols=F, show_rownames=T, show_colnames=F)
dev.off()
#--------------------------------------------------------------
# Figure 4C,  Potential causal genes in metaterm 8 and their functional properties compared with others. 
calCor_AM <- calCor
candidate_HDAC1 <- names(which(calCor[name2] > 1.28)) # name2 # names(which(calCor[name2] > 1.28))
length(candidate_HDAC1)
pdf(file="./figure_plot/Figure4C_LMP_PCAScore_HDAC1_vioplot.pdf")
plot(density(calCor[candidate_HDAC1]),col="red", xlim=c(-4.2,7.2),ylim=c(0,0.45))
lines(density(calCor[setdiff(name3, c('ENSSSCG00000003613',candidate_HDAC1))]),col="black")
dev.off()
t.test(calCor[candidate_HDAC1], calCor[setdiff(name3, c('ENSSSCG00000003613',candidate_HDAC1))], alternative='greater')$p.value
HDAC1_candidate <- gold_gene[intersect(candidate_HDAC1, gold_gene[,3]),]
length(which(HDAC1_candidate[,1] < 0.01))
HDAC1_candidate_sig <- HDAC1_candidate[which(HDAC1_candidate[,1] < 0.01),3]
HDAC1_candidate_sig_ensemble_1Mb[which(as.numeric(HDAC1_candidate_sig_ensemble_1Mb[,2]) < 0),2] <- 1
write.table(HDAC1_candidate_sig_ensemble_1Mb, file = "./figure_plot/HDAC1_candidate_sig_ensemble_1Mb.bed",quote=F,row.names=F,col.names=F,sep='\t')
write.table(HDAC1_candidate_sig_ensemble, file = "./figure_plot/HDAC1_candidate_sig_ensemble.bed",quote=F,row.names=F,col.names=F,sep='\t')
HDAC1_candidate_sig_ensemble_1Mb_noncoding_SNP_infor_sig <- HDAC1_candidate_sig_ensemble_1Mb_noncoding_SNP_infor[which(as.numeric(HDAC1_candidate_sig_ensemble_1Mb_noncoding_SNP_infor[,"element-p-value"]) < 0.05),]
HDAC1_candidate_sig_ensemble_1Mb_noncoding_SNP_infor_sig <- HDAC1_candidate_sig_ensemble_1Mb_noncoding_SNP_infor_sig[which(as.numeric(HDAC1_candidate_sig_ensemble_1Mb_noncoding_SNP_infor_sig[,"SNP-p-value"]) < 0.05),]
HDAC1_candidate_sig_ensemble_1Mb_noncoding_SNP_infor_sig
write.table(HDAC1_candidate_sig_ensemble_1Mb_noncoding_SNP_infor_sig, file = "./figure_plot/HDAC1_candidate_sig_ensemble_1Mb_noncoding_SNP_infor_sig.txt", row.names=F, col.names=T, sep="\t",quote=F)
# save(list = ls(), file = "./figure_data/Figure4.RData")
#--------------------------------------------------------------
########################----Figure 4----########################
########################----Figure 5----########################
# Figure 5, Prediction performance of DeepAnnotation in terms of robustness. 
#--------------------------------------------------------------
# setwd('/data3/mwl/DeepAnnotation_review/gigadb')
load('./figure_data/Figure5.RData')
if(!file.exists('./figure_plot')){
		system(paste0("mkdir ./figure_plot"), intern =TRUE)
}
#--------------------------------------------------------------
# Figure 5A, Optimized hyperparameters of DeepAnnotation on LMP, LMD, BF traits. 
#--------------------------------------------------------------
options(scipen=100)
learning_rate <- c(0.1, 0.01, 0.001, 0.0001, 0.00001)
unit <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130)
regularizer <- c("l1", "l2", "l1_l2","none")
regularizer_rate <- c(0.1, 0.01, 0.001, 0.0001)
momentum <- c(0.99, 0.98, 0.95, 0.9, 0.8)
dropout <- c(0.3, 0.2, 0.1)
unit_ratio <- c("1,1,5,5,1,1", "1,2,3,4,5,6", "1,2,4,5,6,7", "1,3,2,2,3,1", 
	            "1,2,2,5,4,3", "2,2,1,1,2,2")
para_flag <- cbind(c(1,125),c(126,250),c(251,375),c(375,500))
vivo_list_plot <- list()
omic_name <- c("genome","snp","function","network")
for(kk in 1:ncol(para_flag)){
	vivo_list <- list()
	value_result_mean <- NULL
	value_result_median <- NULL
	opt_result <- opt_result_LMD[para_flag[1,kk]:para_flag[2,kk],]
	for(i in 1:length(learning_rate)){
		vivo_list[[paste0("LR_",learning_rate[i])]] <- as.numeric(opt_result[which(opt_result[,14] == as.vector(learning_rate[i])),3])
		value_result_mean <- c(value_result_mean, mean(vivo_list[[paste0("LR_",learning_rate[i])]]))
		value_result_median <- c(value_result_median, median(vivo_list[[paste0("LR_",learning_rate[i])]]))
	}
	for(i in 1:length(unit)){
		vivo_list[[paste0("Unit_",unit[i])]] <- as.numeric(opt_result[which(opt_result[,13] == as.vector(unit[i])),3])
		value_result_mean <- c(value_result_mean, mean(vivo_list[[paste0("Unit_",unit[i])]]))
		value_result_median <- c(value_result_median, median(vivo_list[[paste0("Unit_",unit[i])]]))
	}
	for(i in 1:length(regularizer)){
		vivo_list[[paste0("Regu_",regularizer[i])]] <- as.numeric(opt_result[which(opt_result[,12] == as.vector(regularizer[i])),3])
		value_result_mean <- c(value_result_mean, mean(vivo_list[[paste0("Regu_",regularizer[i])]]))
		value_result_median <- c(value_result_median, median(vivo_list[[paste0("Regu_",regularizer[i])]]))
	}
	for(i in 1:length(regularizer_rate)){
		vivo_list[[paste0("ReRa_",regularizer_rate[i])]] <- as.numeric(opt_result[which(opt_result[,11] == as.vector(regularizer_rate[i])),3])
		value_result_mean <- c(value_result_mean, mean(vivo_list[[paste0("ReRa_",regularizer_rate[i])]]))
		value_result_median <- c(value_result_median, median(vivo_list[[paste0("ReRa_",regularizer_rate[i])]]))
	}
	for(i in 1:length(momentum)){
		vivo_list[[paste0("MOME_",momentum[i])]] <- as.numeric(opt_result[which(opt_result[,10] == as.vector(momentum[i])),3])
		value_result_mean <- c(value_result_mean, mean(vivo_list[[paste0("MOME_",momentum[i])]]))
		value_result_median <- c(value_result_median, median(vivo_list[[paste0("MOME_",momentum[i])]]))
	}
	for(i in 1:length(dropout)){
		vivo_list[[paste0("Drop_",dropout[i])]] <- as.numeric(opt_result[which(opt_result[,9] == as.vector(dropout[i])),3])
		value_result_mean <- c(value_result_mean, mean(vivo_list[[paste0("Drop_",dropout[i])]]))
		value_result_median <- c(value_result_median, median(vivo_list[[paste0("Drop_",dropout[i])]]))
	}
	for(i in 1:length(unit_ratio)){
		vivo_list[[paste0("UR_",unit_ratio[i])]] <- as.numeric(opt_result[which(opt_result[,8] == as.vector(unit_ratio[i])),3])
		value_result_mean <- c(value_result_mean, mean(vivo_list[[paste0("UR_",unit_ratio[i])]]))
		value_result_median <- c(value_result_median, median(vivo_list[[paste0("UR_",unit_ratio[i])]]))
	}
	vivo_list_plot[[paste0(omic_name[kk])]] <- list(vivo_list = vivo_list,
													value_result_mean = value_result_mean,
													value_result_median = value_result_median)
}
write.table(cbind(vivo_list_plot[[1]]$value_result_median, 
	vivo_list_plot[[2]]$value_result_median,
	vivo_list_plot[[3]]$value_result_median,
	vivo_list_plot[[4]]$value_result_median,
	vivo_list_plot[[1]]$value_result_mean, 
	vivo_list_plot[[2]]$value_result_mean, 
	vivo_list_plot[[3]]$value_result_mean, 
	vivo_list_plot[[4]]$value_result_mean), file = "./figure_plot/LMD_Figure5_opt_genome.txt"
	, quote = F, row.names=F, col.names=F, sep="\t")
for(kk in 1:ncol(para_flag)){
	vivo_list <- list()
	value_result_mean <- NULL
	value_result_median <- NULL
	opt_result <- opt_result_BF[para_flag[1,kk]:para_flag[2,kk],]
	for(i in 1:length(learning_rate)){
		vivo_list[[paste0("LR_",learning_rate[i])]] <- as.numeric(opt_result[which(opt_result[,14] == as.vector(learning_rate[i])),3])
		value_result_mean <- c(value_result_mean, mean(vivo_list[[paste0("LR_",learning_rate[i])]]))
		value_result_median <- c(value_result_median, median(vivo_list[[paste0("LR_",learning_rate[i])]]))
	}
	for(i in 1:length(unit)){
		vivo_list[[paste0("Unit_",unit[i])]] <- as.numeric(opt_result[which(opt_result[,13] == as.vector(unit[i])),3])
		value_result_mean <- c(value_result_mean, mean(vivo_list[[paste0("Unit_",unit[i])]]))
		value_result_median <- c(value_result_median, median(vivo_list[[paste0("Unit_",unit[i])]]))
	}
	for(i in 1:length(regularizer)){
		vivo_list[[paste0("Regu_",regularizer[i])]] <- as.numeric(opt_result[which(opt_result[,12] == as.vector(regularizer[i])),3])
		value_result_mean <- c(value_result_mean, mean(vivo_list[[paste0("Regu_",regularizer[i])]]))
		value_result_median <- c(value_result_median, median(vivo_list[[paste0("Regu_",regularizer[i])]]))
	}
	for(i in 1:length(regularizer_rate)){
		vivo_list[[paste0("ReRa_",regularizer_rate[i])]] <- as.numeric(opt_result[which(opt_result[,11] == as.vector(regularizer_rate[i])),3])
		value_result_mean <- c(value_result_mean, mean(vivo_list[[paste0("ReRa_",regularizer_rate[i])]]))
		value_result_median <- c(value_result_median, median(vivo_list[[paste0("ReRa_",regularizer_rate[i])]]))
	}
	for(i in 1:length(momentum)){
		vivo_list[[paste0("MOME_",momentum[i])]] <- as.numeric(opt_result[which(opt_result[,10] == as.vector(momentum[i])),3])
		value_result_mean <- c(value_result_mean, mean(vivo_list[[paste0("MOME_",momentum[i])]]))
		value_result_median <- c(value_result_median, median(vivo_list[[paste0("MOME_",momentum[i])]]))
	}
	for(i in 1:length(dropout)){
		vivo_list[[paste0("Drop_",dropout[i])]] <- as.numeric(opt_result[which(opt_result[,9] == as.vector(dropout[i])),3])
		value_result_mean <- c(value_result_mean, mean(vivo_list[[paste0("Drop_",dropout[i])]]))
		value_result_median <- c(value_result_median, median(vivo_list[[paste0("Drop_",dropout[i])]]))
	}
	for(i in 1:length(unit_ratio)){
		vivo_list[[paste0("UR_",unit_ratio[i])]] <- as.numeric(opt_result[which(opt_result[,8] == as.vector(unit_ratio[i])),3])
		value_result_mean <- c(value_result_mean, mean(vivo_list[[paste0("UR_",unit_ratio[i])]]))
		value_result_median <- c(value_result_median, median(vivo_list[[paste0("UR_",unit_ratio[i])]]))
	}
	vivo_list_plot[[paste0(omic_name[kk])]] <- list(vivo_list = vivo_list,
													value_result_mean = value_result_mean,
													value_result_median = value_result_median)
}
write.table(cbind(vivo_list_plot[[1]]$value_result_median, 
	vivo_list_plot[[2]]$value_result_median,
	vivo_list_plot[[3]]$value_result_median,
	vivo_list_plot[[4]]$value_result_median,
	vivo_list_plot[[1]]$value_result_mean, 
	vivo_list_plot[[2]]$value_result_mean, 
	vivo_list_plot[[3]]$value_result_mean, 
	vivo_list_plot[[4]]$value_result_mean), file = "./figure_plot/BF_Figure5_opt_genome.txt"
	, quote = F, row.names=F, col.names=F, sep="\t")
#--------------------------------------------------------------
# Figure 5B,PCC evaluation of different models based on 5-fold cross-validation experiment for LMD.
cor_cv_result_LMD <- NULL
for(cv in 1:length(cv_flag_LMD)){
	cor_cv_result_LMD <- rbind(cor_cv_result_LMD, cor(pred_phenotype_cal_LMD[cv_flag_LMD[[paste0("cv",cv)]][["test"]],])[1,])
}
cor_cv_pvalue_all_LMD <- NULL
for(l in 9:12){
	cor_cv_pvalue <- NULL
	for(k in 2:12){
		cor_cv_pvalue <- c(cor_cv_pvalue, t.test(cor_cv_result_LMD[,l], cor_cv_result_LMD[,k], paired = T, alternative = 'greater')$p.value)
	}
	cor_cv_pvalue_all_LMD <- rbind(cor_cv_pvalue_all_LMD, cor_cv_pvalue)
}
colnames(cor_cv_pvalue_all_LMD) <- c("rrBLUP", "LightGBM", "KAML", "BLUP", "BayesR", "MBLUP", "BayesRC","Genotype", "SNPAnnotation", "GeneFunction", "Network")
cor_cv_pvalue_all_LMD
model_name <- c("rrBLUP", "LightGBM", "KAML", "BLUP", "BayesR", "MBLUP", "BayesRC","Genotype", "SNPAnnotation", "GeneFunction", "Network")
vivo_list_sparial_pcc_LMD <- list()
vivo_list_sparial_pcc_mean_LMD <- NULL
for(k in 1:length(model_name)){
	vivo_list_sparial_pcc_LMD[[model_name[k]]] <- as.numeric(cor_cv_result_LMD[,k+1])
	vivo_list_sparial_pcc_mean_LMD <- c(vivo_list_sparial_pcc_mean_LMD, mean(cor_cv_result_LMD[,k+1]))
}
library(RColorBrewer)
col_break = brewer.pal(12, "Set3")
col_break = col_break[c(9, 5, 6, 11, 12, 1, 2, 7, 8, 10, 4)]
library(vioplot)
pdf("./figure_plot/LMD_Figure5B_PCC_cv.pdf",height=16, width=12)
vioplot(vivo_list_sparial_pcc_LMD, pch = 20, 
	border = col_break, col = col_break, 
	cex = 0.5,ylab = "PCC", 
	main = "PCC", 
	xaxt = "n", yaxt = "n")
axis(side = 2)
axis(side = 1
	, las = 2
	, at = 1:length(vivo_list_sparial_pcc_LMD))
points(vivo_list_sparial_pcc_mean_LMD, pch = 18, col = 'black', lwd = 6)
dev.off()
pdf("./figure_plot/LMD_Figure5B_PCC.pdf",height=16, width=12)
barplot(cor(pred_phenotype_cal_LMD)[,1][2:12])
dev.off()
#--------------------------------------------------------------
# Figure 5C,PCC evaluation of different models based on 5-fold cross-validation experiment for BF.
cor_cv_result_BF <- NULL
for(cv in 1:length(cv_flag_BF)){
	cor_cv_result_BF <- rbind(cor_cv_result_BF, cor(pred_phenotype_cal_BF[cv_flag_BF[[paste0("cv",cv)]][["test"]],])[1,])
}
cor_cv_pvalue_all_BF <- NULL
for(l in 9:12){
	cor_cv_pvalue <- NULL
	for(k in 2:12){
		cor_cv_pvalue <- c(cor_cv_pvalue, t.test(cor_cv_result_BF[,l], cor_cv_result_BF[,k], paired = T, alternative = 'greater')$p.value)
	}
	cor_cv_pvalue_all_BF <- rbind(cor_cv_pvalue_all_BF, cor_cv_pvalue)
}
colnames(cor_cv_pvalue_all_BF) <- c("rrBLUP", "LightGBM", "KAML", "BLUP", "BayesR", "MBLUP", "BayesRC","Genotype", "SNPAnnotation", "GeneFunction", "Network")
cor_cv_pvalue_all_BF
model_name <- c("rrBLUP", "LightGBM", "KAML", "BLUP", "BayesR", "MBLUP", "BayesRC","Genotype", "SNPAnnotation", "GeneFunction", "Network")
vivo_list_sparial_pcc_BF <- list()
vivo_list_sparial_pcc_mean_BF <- NULL
for(k in 1:length(model_name)){
	vivo_list_sparial_pcc_BF[[model_name[k]]] <- as.numeric(cor_cv_result_BF[,k+1])
	vivo_list_sparial_pcc_mean_BF <- c(vivo_list_sparial_pcc_mean_BF, mean(cor_cv_result_BF[,k+1]))
}
library(RColorBrewer)
col_break = brewer.pal(12, "Set3")
col_break = col_break[c(9, 5, 6, 11, 12, 1, 2, 7, 8, 10, 4)]
library(vioplot)
pdf("./figure_plot/BF_Figure5C_PCC_cv.pdf",height=16, width=12)
vioplot(vivo_list_sparial_pcc_BF, pch = 20, 
	border = col_break, col = col_break, 
	cex = 0.5,ylab = "PCC", 
	main = "PCC", 
	xaxt = "n", yaxt = "n")
axis(side = 2)
axis(side = 1
	, las = 2
	, at = 1:length(vivo_list_sparial_pcc_BF))
points(vivo_list_sparial_pcc_mean_BF, pch = 18, col = 'black', lwd = 6)
dev.off()
pdf("./figure_plot/BF_Figure5C_PCC.pdf",height=16, width=12)
barplot(cor(pred_phenotype_cal_BF)[,1][2:12])
dev.off()
# save(list = ls(), file = "./figure_data/Figure5.RData")
########################----Figure 5----########################
########################----Figure 6----########################
# Figure 6, Evaluation of statistically significant variants.
#--------------------------------------------------------------
# setwd('/data3/mwl/DeepAnnotation_review/gigadb')
load('./figure_data/Figure6.RData')
if(!file.exists('./figure_plot')){
		system(paste0("mkdir ./figure_plot"), intern =TRUE)
}
options(scipen=999)
#--------------------------------------------------------------
# Figure 6A,Training and validation performance of the Genotype model across epochs, based on 5-fold cross-validation. 
Genotype_step_training <- do.call(rbind, strsplit(Genotype_step[grep('Predict training Step', Genotype_step[,1]),1], "( )|:"))
Genotype_step_star <- which(as.numeric(Genotype_step_training[,4]) == 1)
Genotype_step_end <- Genotype_step_star-1
Genotype_step_end <- Genotype_step_end[-1]
Genotype_step_end <- c(Genotype_step_end, nrow(Genotype_step_training))
Genotype_step_flag <- rbind(Genotype_step_star, Genotype_step_end)
Genotype_step_range <- Genotype_step_flag[2,] - Genotype_step_flag[1,] + 1
Genotype_train_loss <- NULL
for(i in 1:5){
	Genotype_train_loss <- cbind(Genotype_train_loss,as.numeric(Genotype_step_training[Genotype_step_star[i]:(Genotype_step_star[i]+min(Genotype_step_range)-1),7]))
}
Genotype_train_loss <- cbind(Genotype_train_loss, apply(Genotype_train_loss, 1, mean))
Genotype_step_validing <- do.call(rbind, strsplit(Genotype_step[grep('Predict validation Step', Genotype_step[,1]),1], "( )|:"))
Genotype_step_star <- which(as.numeric(Genotype_step_validing[,4]) == 1)
Genotype_step_end <- Genotype_step_star-1
Genotype_step_end <- Genotype_step_end[-1]
Genotype_step_end <- c(Genotype_step_end, nrow(Genotype_step_validing))
Genotype_step_flag <- rbind(Genotype_step_star, Genotype_step_end)
Genotype_step_range <- Genotype_step_flag[2,] - Genotype_step_flag[1,] + 1
Genotype_step_range
Genotype_valid_loss <- NULL
for(i in 1:5){
	Genotype_valid_loss <- cbind(Genotype_valid_loss,as.numeric(Genotype_step_validing[Genotype_step_star[i]:(Genotype_step_star[i]+min(Genotype_step_range)-1),7]))
}
Genotype_valid_loss <- cbind(Genotype_valid_loss, apply(Genotype_valid_loss, 1, mean))
train_min_flag <- NULL
train_loss_bk <- 99999999
temp_num <- 0
for(i in 1:nrow(Genotype_train_loss)){
	temp <- NULL
	if(Genotype_train_loss[i,6] < train_loss_bk){
		train_loss_bk <- Genotype_train_loss[i,6]
		train_min_flag <- rbind(train_min_flag, cbind(i, Genotype_train_loss[i,6], temp_num))
		temp_num <- 0	
	}else{
		temp_num <- temp_num + 1
	}
}
which.max(train_min_flag[,3])
train_min_flag[(which.max(train_min_flag[,3])-1):(which.max(train_min_flag[,3])+1),]
valid_min_flag <- NULL
valid_loss_bk <- 99999999
temp_num <- 0
for(i in 1:nrow(Genotype_valid_loss)){
	temp <- NULL
	if(Genotype_valid_loss[i,6] < valid_loss_bk){
		valid_loss_bk <- Genotype_valid_loss[i,6]
		valid_min_flag <- rbind(valid_min_flag, cbind(i, Genotype_valid_loss[i,6], temp_num))
		temp_num <- 0	
	}else{
		temp_num <- temp_num + 1
	}
}
which.max(valid_min_flag[,3])
valid_min_flag[(which.max(valid_min_flag[,3])-1):(which.max(valid_min_flag[,3])+1),]
pdf("./figure_plot/LMP_Figure6A_Genotype_train_barplot.pdf"
	,width = 12, height = 4)
a1 <- train_min_flag[which(train_min_flag[,3] > 0),3]
barplot(a1, las = 2)
dev.off()
pdf("./figure_plot/LMP_Figure6A_Genotype_valid_barplot.pdf"
	,width = 12, height = 4)
a2 <- valid_min_flag[which(valid_min_flag[,3] > 0),3]
barplot(a2, las = 2)
dev.off()
library(RColorBrewer)
col_loss <- brewer.pal(8, "Set2")[c(8,7)]
pdf("./figure_plot/LMP_Figure6A_Genotype_loss_1_100.pdf"
	)
plot(Genotype_train_loss[1:100,6], type = "l", lwd = 2, col = col_loss[1], ylim=c(0,5000))
lines(Genotype_valid_loss[1:100,6], lwd = 2, col = col_loss[2])
abline(v = which.min(Genotype_train_loss[1:100,6]), col = "black")
abline(v = which.min(Genotype_valid_loss[1:100,6]), col = "black")
dev.off()
pdf("./figure_plot/LMP_Figure6A_Genotype_loss_100_924.pdf"
	)
plot(Genotype_train_loss[100:nrow(Genotype_train_loss),6], type = "l", lwd = 2, col = col_loss[1],ylim=c(4,20))
lines(Genotype_valid_loss[100:nrow(Genotype_train_loss),6], lwd = 2, col = col_loss[2])
abline(v = which.min(Genotype_train_loss[100:nrow(Genotype_train_loss),6]), col = "black")
abline(v = which.min(Genotype_valid_loss[100:nrow(Genotype_train_loss),6]), col = "black")
dev.off()
which.min(Genotype_train_loss[,6]) # 907
which.min(Genotype_valid_loss[,6]) # 564
#--------------------------------------------------------------
# Figure 6B,Estimated heritability of significant SNPs (P-value < 0.05) identified by different models for the LMP trait.
library(RColorBrewer)
col_break = brewer.pal(12, "Set3")
col_break = col_break[c(9, 5, 6, 11, 12, 1, 2)]
heri_real_result_plot <- heri_real_result[c(5,3,4,2,1,7,6)]
names(heri_real_result_plot) <- c("rrBLUP", "BLUP", "BayesR", "Genotype", "Network","metaG","metaN")
pdf("./figure_plot/LMP_Figure6B_heri_real.pdf",height=16, width=12)
barplot(heri_real_result_plot, col = col_break)
dev.off()
# ----------------------------------------------- 
# Figure 6C, Genomic prediction accuracy (measure by Pearson correlation coefficient, PCC) calculated using rrBLUP model for different SNP sets: I) 11,084 potential causal variants identified by DeepAnnotation (including 10,041 non-coding or cis-regulatory variants and 1,043 coding or RNA secondary structure-related variants with P-value < 0.01 from the backtracking strategy). II) 644,964 SNPs located within 1Mb of the 11,084 DeepAnnotation predicted causal variants. III) 4,850 SNPs with GWAS P-value < 1.0e-6. IV) 11,083 SNPs with GWAS P-value < 1.0e-5. V) 24,833 SNPs with GWAS P-value < 1.0e-4. VI) Ten sets of 11,084 SNPs randomly selected across the genome (mean PCC reported).
PCC_result <- t(t(cor(pred_phenotype[test_sample,])[1,]))
PCC_result_plot <- c("SNP_random"=mean(PCC_result[2:11,1]),(PCC_result[12:16,1]))
pdf("./figure_plot/LMP_Figure6C_PCC_plot.pdf",height=16, width=12)
barplot(PCC_result_plot, col = col_break)
dev.off()
########################----Figure 6----########################
# ******************** Main Figures******************************** #

# *********************Supplementary Figures******************************** #
########################----Supplementary Figure S2----########################
# Supplementary Figure S2, Average PCC from different sets of SNPs with rrBLUP model.
#--------------------------------------------------------------
# setwd('/data3/mwl/DeepAnnotation_review/gigadb')
load('./figure_data/FigureS2.RData')
if(!file.exists('./figure_plot')){
		system(paste0("mkdir ./figure_plot"), intern =TRUE)
}
pdf(file="./figure_plot/Supp_FigureS2_pcc_barplot.pdf")
barplot(round(apply(rr_all[,2:5], 2, mean),3), col = 1:4, ylim=c(0.4,0.42))
dev.off()
pdf(file="./figure_plot/Supp_FigureS2_log2num_barplot.pdf")
barplot(log2(c(11633164,32451,344093,754033)), col=c(2,3,6,5))
dev.off()
# save(list = ls(), file = "./figure_data/FigureS2.RData")
########################----Supplementary Figure S2----########################
########################----Supplementary Figure S3----########################
# Supplementary Figure S3,  Receiver operating characteristic (ROC) curves of ATAC data learned from DeepSEA model for different tissues of different species.
# setwd('/data3/mwl/DeepAnnotation_review/gigadb')
load('./figure_data/FigureS3.RData')
if(!file.exists('./figure_plot')){
		system(paste0("mkdir ./figure_plot"), intern =TRUE)
}
library(RColorBrewer)
library(ROCR)
library(pROC)
tissues <- colnames(predict)[5:ncol(predict)]
col_ROC = brewer.pal(9, "Set1")
species <- c("Duroc","ES","LW","MS")
pdf("./figure_plot/Supp_FigureS3_ATAC_ROC.pdf"
	,width=8, height=8)
par(mfcol = c(2,2))
AUC1 <- NULL
i = 1
tempP <- predict[,grep(species[i],colnames(predict))]
tempB <- b1[,grep(species[i],colnames(b1))]
pred <- prediction(as.numeric(tempP[,1]), tempB[,1])	
perf <- performance(pred,'tpr','fpr')
plot(perf, col = col_ROC[1], xlim = c(0,1), ylim = c(0,1))
legend(x = 0.5, y = 1, border = NA
     , bty = "n", col=col_ROC[1:ncol(tempP)], lty=c(rep(1,length(tissues)))
	, legend=colnames(tempP))
AUC1 <- c(AUC1, unlist(slot(performance(pred,'auc'),"y.values")))
for(ii in 2:ncol(tempP)){
	pred <- prediction(as.numeric(tempP[,ii]), tempB[,ii])	
	perf <- performance(pred,'tpr','fpr')
	plot(perf, col = col_ROC[ii], add = T)
	AUC1 <- c(AUC1, unlist(slot(performance(pred,'auc'),"y.values")))

}
for(i in 2:length(species)){
	tempP <- predict[,grep(species[i],colnames(predict))]
	tempB <- b1[,grep(species[i],colnames(b1))]	
	pred <- prediction(as.numeric(tempP[,1]), tempB[,1])	
	perf <- performance(pred,'tpr','fpr')
	plot(perf, col = col_ROC[1], xlim = c(0,1), ylim = c(0,1))
	legend(x = 0.5, y = 1, border = NA
	     , bty = "n", col=col_ROC[1:ncol(tempP)], lty=c(rep(1,length(tissues)))
		, legend=colnames(tempP))
	AUC1 <- c(AUC1, unlist(slot(performance(pred,'auc'),"y.values")))
	for(ii in 2:ncol(tempP)){
		pred <- prediction(as.numeric(tempP[,ii]), tempB[,ii])	
		perf <- performance(pred,'tpr','fpr')
		plot(perf, col = col_ROC[ii], add = T)
		AUC1 <- c(AUC1, unlist(slot(performance(pred,'auc'),"y.values")))
	}
}
dev.off()
names(AUC1) <- tissues
# save(list = ls(), file = "./figure_data/FigureS3.RData")
########################----Supplementary Figure S3----########################
# *********************Supplementary Figures******************************** #

