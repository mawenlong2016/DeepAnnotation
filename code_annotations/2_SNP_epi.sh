awk '{print $1"\t"$2"\t"$3"\t""Duroc_Fat"}' merge_Duroc_Fat.bed > 2_TFBind.bed
awk '{print $1"\t"$2"\t"$3"\t""Duroc_Heart"}' merge_Duroc_Heart.bed >> /2_TFBind.bed
awk '{print $1"\t"$2"\t"$3"\t""Duroc_Liver"}' merge_Duroc_Liver.bed >> 2_TFBind.bed
awk '{print $1"\t"$2"\t"$3"\t""Duroc_Muscle"}' merge_Duroc_Muscle.bed >> 2_TFBind.bed
awk '{print $1"\t"$2"\t"$3"\t""Duroc_Spleen"}' GSM4256382_Duroc_2W_1_Spleen_ATAC.narrowPeak >> 2_TFBind.bed
awk '{print $1"\t"$2"\t"$3"\t""ES_Fat"}' merge_ES_Fat.bed >> 2_TFBind.bed
awk '{print $1"\t"$2"\t"$3"\t""ES_Heart"}' GSM4256385_ES_2W_1_Heart_ATAC.narrowPeak >> 2_TFBind.bed
awk '{print $1"\t"$2"\t"$3"\t""ES_Liver"}' GSM4256387_ES_2W_2_Liver_ATAC.narrowPeak >> 2_TFBind.bed
awk '{print $1"\t"$2"\t"$3"\t""ES_Muscle"}' merge_ES_Muscle.bed >> 2_TFBind.bed
awk '{print $1"\t"$2"\t"$3"\t""ES_Spleen"}' GSM4256390_ES_2W_1_Spleen_ATAC.narrowPeak >> 2_TFBind.bed
awk '{print $1"\t"$2"\t"$3"\t""LW_Fat"}' merge_LW_Fat.bed >> 2_TFBind.bed
awk '{print $1"\t"$2"\t"$3"\t""LW_Liver"}' GSM4256393_LW_2W_1_Liver_ATAC.narrowPeak >> 2_TFBind.bed
awk '{print $1"\t"$2"\t"$3"\t""LW_Muscle"}' merge_LW_Muscle.bed >> 2_TFBind.bed
awk '{print $1"\t"$2"\t"$3"\t""MS_Fat"}' merge_MS_Fat.bed >> 2_TFBind.bed
awk '{print $1"\t"$2"\t"$3"\t""MS_Liver"}' GSM4256400_MS_2W_2_Liver_ATAC.narrowPeak >> 2_TFBind.bed
awk '{print $1"\t"$2"\t"$3"\t""MS_Muscle"}' GSM4256401_MS_2W_2_Muscle_ATAC.narrowPeak >> 2_TFBind.bed

cd variant
conda activate DeepSEA
python ../../code/2_DeepSEA-pre.py -c 1_ChrPos.txt -t 2_TFBind.bed -f 3_Feature.bed
docker run --runtime=nvidia -it --net=host --name DeepAnnotation-DeepSEA -v /data/mawenlong/DeepAnnotation/workflow:/data deepseat:20211229 /bin/bash
cd /data/data/variant/
export CUDA_VISIBLE_DEVICES=0
sh ../../code_deepsea/8_trainDeepSEA.sh
# return back to conda DeepSEAT
python ../../code_deepsea/9_predictDeepSEA-pre.py 4_testBin.bed ./predict_results test 
# return back to docker
python ../../code_deepsea/9_predictDeepSEA.py 4_testBin.bed ./predict_results test 
# return back to conda DeepSEAT
python ../../code_deepsea/9_predictDeepSEA-end.py 4_testBin.bed ./predict_results test 
Rscript ../../code_deepsea/10_evaluation.R

# predict SNP epigenome function from deepsea model
awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' duroc.vcf.bed > variant.vcf
python ../../code_deepsea/9_predictDeepSEA-pre.py variant.vcf
#Rscript  ../../code_deepsea/6_extractSeq.R-pig variant.vcf 1100
python ../../code_deepsea/9_predictDeepSEA.py variant.vcf
python ../../code_deepsea/9_predictDeepSEA-end.py variant.vcf
python ../../code_deepsea/deepsea_predict/3_h5ToOutput.py variant.vcf variant.vcf.wt1100.fasta.ref.h5.pred.h5 variant.vcf.mut1100.fasta.ref.h5.pred.h5