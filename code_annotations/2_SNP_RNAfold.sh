
# install ViennaRNA RNAfold
# http://rpm.pbone.net/info_idpl_54834908_distro_centos6_com_ViennaRNA-1.8.5-13.7.x86_64.rpm.html
# sudo rpm -ivh ViennaRNA-1.8.5-13.7.x86_64.rpm
cd /data/mawenlong/DeepAnnotation/workflow/data/variant
SNP_Duroc_Gene_vep_coding.vcf
Rscript ../../code_deepsea/6_extractSeq.R-pig SNP_Duroc_Gene_vep_coding.vcf 
python ../../code/2_extractMuFasta.py SNP_Duroc_Gene_vep_coding.vcf.wt1100.fasta 1100

RNAfold -d2 -noLP < SNP_Duroc_Gene_vep_coding.vcf.wt1100.fasta.ref.fasta > vep_coding_ref.out
RNAfold -d2 -noLP < SNP_Duroc_Gene_vep_coding.vcf.mut1100.fasta.ref.fasta > vep_coding_alt.out

con <- file("vep_coding_ref.out", "r")
vep_coding_result <- NULL
vep_coding_ref <- NULL
oneline = readLines(con, n = 1)	
#print(oneline)	
t1 <- Sys.time()
while(1){
	if(length(oneline) == 0){
        break
    }

	if(length(grep(">",oneline))){		
		count <- 0	
		num <- 0
		temp1 <- unlist(strsplit(oneline, split = ">|_"))
		while(1){
			oneline = readLines(con, n = 1)
			if(length(oneline) == 0){
		        break
		    }
#			print(oneline)
			if(length(grep(">",oneline))){
				break
			}else if(length(grep(")",oneline))){
				temp <- unlist(strsplit(oneline, split = "\\(|)"))
				count <- count + as.numeric(temp[length(temp)])	
				num <- num + 1			
			}else{

			}
		}
		vep_coding_ref <- rbind(vep_coding_ref,
									cbind(temp1[4],as.numeric(temp1[6])-550,".",temp1[2],temp1[3],count,num))

	}	
}
close(con)
t2 <- Sys.time()
print(t2 - t1)

con <- file("vep_coding_alt.out", "r")
vep_coding_alt <- NULL
oneline = readLines(con, n = 1)	
#print(oneline)	
t1 <- Sys.time()
while(1){
	if(length(oneline) == 0){
        break
    }

	if(length(grep(">",oneline))){		
		count <- 0	
		num <- 0
		temp1 <- unlist(strsplit(oneline, split = ">|_"))
		while(1){
			oneline = readLines(con, n = 1)
			if(length(oneline) == 0){
		        break
		    }
#			print(oneline)
			if(length(grep(">",oneline))){
				break
			}else if(length(grep(")",oneline))){
				temp <- unlist(strsplit(oneline, split = "\\(|)"))
				count <- count + as.numeric(temp[length(temp)])	
				num <- num + 1			
			}else{

			}
		}
		vep_coding_alt <- rbind(vep_coding_alt,
									cbind(temp1[4],as.numeric(temp1[6])-550,".",temp1[2],temp1[3],count,num))

	}	
}
close(con)
t2 <- Sys.time()
print(t2 - t1)

save(vep_coding_ref, vep_coding_alt, file = "vep_coding_result.RData")
