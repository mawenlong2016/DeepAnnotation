
load("CompFunc.RData")
load("CompFunc.R")
source("Ttest_MAP.R")
source("EGP.R")
IDpool <- list(1:5000,5001:10000,10001:15016)
for(flag in 1:length(IDpool)){
	preflag <- range(IDpool[[flag]])[1]
	nowflag <- range(IDpool[[flag]])[2]
	t1 <- Sys.time()
	print(paste0("Function complement of ",preflag, "--", nowflag, " begins !"))
	CompFunc_result <- CompFunc(FuncTerms=FuncTerms, cpu=50, AmplitudePCA=AmplitudePCA
		, AmplitudeICA=AmplitudeICA, AmplitudeNMF=AmplitudeNMF, rangeID=c(preflag:nowflag))
	print(paste0("Function complement of ",preflag, "--", nowflag, " ends !"))
	t2 <- Sys.time()
	print(t2 - t1)
	save(CompFunc_result, file = paste0("CompFunc_result_"
										,preflag, "_",nowflag,".RData")
		, version = 2)

}