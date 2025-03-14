# Pig genes (Sscrofa11.1)
# http://www.ensembl.org/biomart/martview, download at 20220407
# 
work_dir = "./"
GO_annotation <- read.delim(paste0(work_dir,"gene_annotation.txt-20220407")
	, sep = "\t", quote = "", header = T, stringsAsFactors = F)
GenePosition <- GO_annotation[,c(1,10,11,12,13)]
GOAnnotate <- GO_annotation[,c(1,19,21,22,20)]
GenePosition <- GenePosition[!duplicated(GenePosition),]
GOAnnotate <- GOAnnotate[which(GOAnnotate[,2]!=""),]
GOAnnotate <- GOAnnotate[!duplicated(GOAnnotate),]
colnames(GOAnnotate) <- c("Gene", "GO", "Code", "Domain", "Name")
colnames(GenePosition) <- c("Gene", "Chr", "Type", "Start", "End")
GOAnnotate <- GOAnnotate[which(GOAnnotate[,3]!=""),]
GOAnnotate[which(GOAnnotate[,2] == "GO:0140672"),4] <- "cellular_component"
GOAnnotate[which(GOAnnotate[,2] == "GO:0140672"),5] <- "ATAC complex"
GOAnnotate[which(GOAnnotate[,2] == "GO:0140684"),4] <- "molecular_function"
GOAnnotate[which(GOAnnotate[,2] == "GO:0140684"),5] <- "histone H3-tri/dimethyl-lysine-9 demethylase activity"
GOAnnotate[which(GOAnnotate[,2] == "GO:0140662"),4] <- "molecular_function"
GOAnnotate[which(GOAnnotate[,2] == "GO:0140662"),5] <- "ATP-dependent protein folding chaperone"
GOAnnotate[which(GOAnnotate[,2] == "GO:0006288"),4] <- "biological_process"
GOAnnotate[which(GOAnnotate[,2] == "GO:0006288"),5] <- "base-excision repair, DNA ligation"
GOAnnotate[which(GOAnnotate[,2] == "GO:0140658"),4] <- "molecular_function"
GOAnnotate[which(GOAnnotate[,2] == "GO:0140658"),5] <- "ATP-dependent chromatin remodeler activity"
GOAnnotate[which(GOAnnotate[,2] == "GO:0140682"),4] <- "molecular_function"
GOAnnotate[which(GOAnnotate[,2] == "GO:0140682"),5] <- "histone H3-di/monomethyl-lysine-4 FAD-dependent demethylase activity"
GOAnnotate[which(GOAnnotate[,2] == "GO:0008934"),4] <- "molecular_function"
GOAnnotate[which(GOAnnotate[,2] == "GO:0008934"),5] <- "inositol monophosphate 1-phosphatase activity"
GOAnnotate[which(GOAnnotate[,2] == "GO:0004828"),4] <- "molecular_function"
GOAnnotate[which(GOAnnotate[,2] == "GO:0004828"),5] <- "serine-tRNA ligase activity"
GOAnnotate[which(GOAnnotate[,2] == "GO:0120319"),4] <- "molecular_function"
GOAnnotate[which(GOAnnotate[,2] == "GO:0120319"),5] <- "long-chain fatty acid omega-1 hydroxylase activity"
GOAnnotate[which(GOAnnotate[,2] == "GO:0120325"),4] <- "molecular_function"
GOAnnotate[which(GOAnnotate[,2] == "GO:0120325"),5] <- "NuRD complex binding"
GOAnnotate[which(GOAnnotate[,2] == "GO:0106386"),4] <- "molecular_function"
GOAnnotate[which(GOAnnotate[,2] == "GO:0106386"),5] <- "(3R)-hydroxyacyl-CoA dehydrogenase (NAD) activity"
GOAnnotate[which(GOAnnotate[,2] == "GO:0106383"),4] <- "biological_process"
GOAnnotate[which(GOAnnotate[,2] == "GO:0106383"),5] <- "dAMP salvag"
verifyCode = c("IDA","IMP","IPI", "IGI", "IEP", "TAS")
GTOGOVerify_evaluation <- GOAnnotate[which(is.na(match(GOAnnotate[,3], verifyCode)) == F),]
GTOGOVerify_evaluation <- GTOGOVerify_evaluation[!duplicated(GTOGOVerify_evaluation),]
save(GTOGOVerify_evaluation,file="GTOGOVerify_evaluation.RData")
GTOGOVerify <- GOAnnotate[!duplicated(GOAnnotate),]
# kegg pathway annotation from KEGG
# ssc00001.keg
# https://www.genome.jp/kegg-bin/download_htext?htext=ssc00001.keg&format=htext&filedir=
# download at 20220407
# extract kegg pathways

# awk '$1=="C" &&$NF~"PATH:" || $1=="D"' ssc00001.keg-20220407 | grep -P "PATH|\tK" | sed 's#^C[[:space:]]*##; s#^D[[:space:]]*##; s# \[#\t\[#; s# #\t#' | awk 'BEGIN{FS=OFS="\t"}{if($NF~"PATH:") a=$3"\t"$2;else print $1,a}' | awk 'BEGIN{FS=OFS="\t"}{a[$1]=a[$1]$2",";b[$1]=b[$1]"|"$3; next}END{for (i in a) print i,a[i],b[i]}' | sed 's#,\t|#\t#; s#\[PATH:#path:#g; s#\]##g' > KEGG_ssc.txt
# perl -alne '{if(/^C/){/PATH:ssc(\d+)/;$kegg=$1}else{print "ssc$kegg\t$F[1]" if /^D/ and $kegg;}}' ssc00001.keg-20220407 > kegg2gene_ssc.txt
# awk '$1=="C" &&$NF~"PATH:"' ssc00001.keg-20220407 | awk '{print $1,"ssc"$2,$3}' > KEGG_names.txt
kegg2gene <- read.table(paste0(work_dir,"kegg2gene_ssc.txt"))
write.table(kegg2gene[,2], file = paste0(work_dir, "entrezID.txt"), )
KEGG_names <- read.table(paste0(work_dir,"KEGG_names.txt"))
rownames(KEGG_names) <- KEGG_names[,2]
# https://david.ncifcrf.gov/
# upload entrezID at 20220407 
# Conversion Summary
# ID Count	In DAVID DB	Conversion
# 8205	Yes	Successful
# 630	Yes	None
# 24	No	None
# 630	Ambiguous	Pending
entrezID_gene <- read.table(paste0(work_dir,"entrezID_geneID.txt"), sep = "\t", header = T, stringsAsFactors=F)
a = entrezID_gene
a1 = a[!duplicated(a[,1]),]
rownames(a1) <- a1[,1]
a2 = a[duplicated(a[,1]),]
a = a2
colnames(kegg2gene) <- c("KEGG","entrezID")
kegg2gene <- cbind(kegg2gene, "geneID"=a1[kegg2gene[,2],2])
temp3 <- NULL
while(1){
	print(dim(a2))
	a1 = a[!duplicated(a[,1]),]
	rownames(a1) <- a1[,1]
	a2 = a[duplicated(a[,1]),]
	temp1 <- kegg2gene[which(is.na(match(kegg2gene[,2], a2[,1])) == F),1:2]
	temp1 <- cbind(temp1, "geneID"=a2[temp1[,2],2])
	temp3 <- rbind(temp3, temp1)	
	a = a2
	if(nrow(a2) == 0)
		break

}
kegg2gene <- rbind(kegg2gene, temp3)
kegg2gene <- as.data.frame(kegg2gene, stringsAsFactors=F)
kegg2gene <- kegg2gene[!duplicated(kegg2gene),]
# gene functions prediction
work_dir = ""
geneExp <- read.table('all_TPM', header=T, stringsAsFactors=F, sep = "\t")
geneExp_back = geneExp
geneExp <- geneExp[,c(7:ncol(geneExp))]
geneExp <- log2(geneExp + 1)
sampleNum <- ncol(geneExp)
nonExpressedNum = ceiling(sampleNum * 0.1)
nonzero_num <- apply(geneExp, 1, function(ll){
	return(length(which(ll > 0)))
	})
geneExp <- geneExp[which(nonzero_num >= 1),]
FuncTerms <- rbind(cbind(geneID=GTOGOVerifyHireWhole[,1],term=GTOGOVerifyHireWhole[,2], Source="GO"),
	cbind(geneID=kegg2gene[,3],term=kegg2gene[,1],Source="KEGG"))
###########################
library(GO.db)
GOID <- unique(FuncTerms[which(FuncTerms[,3]=="GO"),2])
GOID <- cbind(GOID=GOID, Def=Term(GOID))
KEGG <- unique(FuncTerms[which(FuncTerms[,3]=="KEGG"),2])
KEGG_names <- read.table(paste0(work_dir,"KEGG_names.txt"), stringsAsFactors=F)
rownames(KEGG_names) <- KEGG_names[,2]
KEGG <- cbind(KEGG=KEGG, Def=KEGG_names[KEGG,3])
TermDef <- rbind(GOID, KEGG)
rownames(TermDef) <- TermDef[,1]
save(TermDef, file=paste0(work_dir,"TermDef.RData"))
############################################
# decomposition gene expression
source(paste0(work_dir,"MFdecomposition.R"))
rownames(geneExp) <- geneExp_back[,1]
t1 <- Sys.time()
genePCA <- MFdecomposition(geneExp = geneExp, cpu = 30, methods = "PCA")
t2 <- Sys.time()
print(t2 - t1)
t1 <- Sys.time()
geneICA <- MFdecomposition(geneExp = geneExp, cpu = 30, methods = "ICA", factors = dim(genePCA$Amplitude)[2])
t2 <- Sys.time()
print(t2 - t1)
t1 <- Sys.time()
geneNMF <- MFdecomposition(geneExp = geneExp, cpu = 30, methods = "NMF", factors = dim(genePCA$Amplitude)[2])
t2 <- Sys.time()
print(t2 - t1)
save(list = ls(), file= paste0(work_dir,"2_easyMF.RData-20220411"))

AmplitudePCA <- genePCA$Amplitude
AmplitudeICA <- geneICA$Amplitude
AmplitudeNMF_back <- geneNMF$Amplitude
AmplitudeNMF <- geneNMF$Amplitude
rownames(AmplitudePCA) <- rownames(AmplitudeICA) <- rownames(AmplitudeNMF)<- rownames(AmplitudeNMF_back) <- rownames(geneExp)
a1 = apply(AmplitudeNMF, 2, range)
AmplitudeNMF <- AmplitudeNMF[,which(a1[2,] > 0)]
FuncTerms <- FuncTerms[!duplicated(FuncTerms),]
save(FuncTerms, file="FuncTerms_all.RData")
FuncTerms <- FuncTerms[which(is.na(match(FuncTerms[,1], rownames(geneExp))) == F),]
save(AmplitudePCA, AmplitudeICA, AmplitudeNMF, FuncTerms
	,  file = paste0(work_dir,"CompFunc.RData"), version = 2)
save(FuncTerms, file = paste0(work_dir,"FuncTerms.RData"))
save(GTOGOVerifyHireWhole, file = paste0(work_dir,"GTOGOVerifyHireWhole.RData"))

Rscript run_CompFunc.R
Rscript transform_function.R