load("CompFunc.RData")
#############################
# load 1-6300
preflag <- 1
nowflag <- 6300
PCAScore=NULL
ICAScore=NULL
NMFScore=NULL
MFScore=NULL
ENScore=NULL
RandomScore=NULL
PCARank=NULL
ICARank=NULL
NMFRank=NULL
MFRank=NULL
ENRank=NULL
RandomRank=NULL

positiveGene <- unique(FuncTerms[,1])
geneSetTerm <- unique(FuncTerms[,2])
for(geneSetTermIdx in preflag:nowflag){
	temp <- read.table(paste0("./",geneSetTerm[geneSetTermIdx]
		,".txt"))
	rownames(temp) <- temp[,1]
	temp <- temp[,-1]
	PCAScore <- cbind(PCAScore, temp[,1])
	ICAScore <- cbind(ICAScore, temp[,2])
	NMFScore <- cbind(NMFScore, temp[,3])
	MFScore <- cbind(MFScore, temp[,4])
	ENScore <- cbind(ENScore, temp[,5])
	RandomScore <- cbind(RandomScore, temp[,6])
	PCARank <- cbind(PCARank, temp[,7])
	ICARank <- cbind(ICARank, temp[,8])
	NMFRank <- cbind(NMFRank, temp[,9])
	MFRank <- cbind(MFRank, temp[,10])
	ENRank <- cbind(ENRank, temp[,11])
	RandomRank <- cbind(RandomRank, temp[,12])
}
colnames(PCAScore) <- colnames(ICAScore) <- colnames(NMFScore) <- colnames(MFScore)<- colnames(ENScore) <- colnames(RandomScore) <- colnames(PCARank) <- colnames(ICARank)<- colnames(NMFRank) <- colnames(MFRank) <- colnames(ENRank) <- colnames(RandomRank)<- geneSetTerm[preflag:nowflag]
CompFunc_result <- list(PCAScore=PCAScore, ICAScore=ICAScore, NMFScore=NMFScore
              , MFScore=MFScore, ENScore=ENScore, RandomScore=RandomScore
              , PCARank=PCARank,ICARank=ICARank,NMFRank=NMFRank
              , MFRank=MFRank, ENRank=ENRank, RandomRank=RandomRank)
save(CompFunc_result, file = paste0("CompFunc_result_"
	,preflag, "_",nowflag,".RData"), version = 2)
#############################
# combined all result togather
PCAScore=NULL
ICAScore=NULL
NMFScore=NULL
MFScore=NULL
ENScore=NULL
RandomScore=NULL
PCARank=NULL
ICARank=NULL
NMFRank=NULL
MFRank=NULL
ENRank=NULL
RandomRank=NULL
IDpool <- c(seq(1,10171,900), 10140)
for(flag in 2:length(IDpool)){ # 2:length(IDpool) # original
	t1 <- Sys.time()
	preflag <- IDpool[flag-1]
	nowflag <- IDpool[flag]-1	
	load(file = paste0("CompFunc_result_"
										,preflag, "_",nowflag,".RData")
		)
	PCAScore <- cbind(PCAScore, CompFunc_result$PCAScore)
	ICAScore <- cbind(ICAScore, CompFunc_result$ICAScore)
	NMFScore <- cbind(NMFScore, CompFunc_result$NMFScore)
	MFScore <- cbind(MFScore, CompFunc_result$MFScore)
	ENScore <- cbind(ENScore, CompFunc_result$ENScore)
	RandomScore <- cbind(RandomScore, CompFunc_result$RandomScore)
	PCARank <- cbind(PCARank, CompFunc_result$PCARank)
	ICARank <- cbind(ICARank, CompFunc_result$ICARank)
	NMFRank <- cbind(NMFRank, CompFunc_result$NMFRank)
	MFRank <- cbind(MFRank, CompFunc_result$MFRank)
	ENRank <- cbind(ENRank, CompFunc_result$ENRank)
	RandomRank <- cbind(RandomRank, CompFunc_result$RandomRank)
	t2 <- Sys.time()
	print(t2 - t1)

}
save(PCAScore, file = paste0("PCAScore.RData"))
save(ICAScore, file = paste0("ICAScore.RData"))
save(NMFScore, file = paste0("NMFScore.RData"))
save(MFScore, file = paste0("MFScore.RData"))
save(ENScore, file = paste0("ENScore.RData"))
save(RandomScore, file = paste0("RandomScore.RData"))
save(PCARank, file = paste0("PCARank.RData"))
save(ICARank, file = paste0("ICARank.RData"))
save(NMFRank, file = paste0("NMFRank.RData"))
save(MFRank, file = paste0("MFRank.RData"))
save(ENRank, file = paste0("ENRank.RData"))
save(RandomRank, file = paste0("RandomRank.RData"))