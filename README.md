# ___DeepAnnotation: A novel interpretable deep learning-based genomic selection model that integrates comprehensive functional annotations___ <br>

The python package 'DeepAnnotation' can be used to perform genomic selection (GS), which is a promising
breeding strategy for agricultural breeding. DeepAnnotation predicts phenotypes from comprehensive multi-omics functional annotations with interpretable deep learning framework. The effectiveness
of DeepAnnotation has been demonstrated in predicting three pork production traits (lean meat percentage at 100 kg [LMP], loin muscle depth at 100 kg [LMD], back fat thickness at 100 kg [BF]) on a population
of 1940 Duroc boars with 11633164 SNPs (_Sus scrofa_) from [GigaDB](https://ftp.cngb.org/pub/gigadb/pub/10.5524/100001_101000/100894). The comprehensive functional annotation data and data used for training DeepAnnotation are available at [**zenodo**](https://zenodo.org/records/10215565).
<br>
## Version and download <br>
* [Version 1.0](https://github.com/mawenlong2016/DeepAnnotation/archive/refs/heads/main.zip) -First version released on Oct, third, 2023<br>
## DeepAnnotation dependencies
DeepAnnotation denpends on the following software environments：<br>
1. [Python](https://www.python.org) - The python (version >=3.6) is needed. <br>
2. We suggest [conda](https://docs.conda.io/projects/miniconda/en/latest/) to install the required packages. <br>
## DeepAnnotation Quick Installation
### Installed by personal computer
1. Build DeepAnnotation environment by conda
```bash
$ conda create -n DeepAnnotation python=3.6 numpy==1.19.2 tensorflow-gpu==2.2.0
$ conda activate DeepAnnotation
$ conda install matplotlib
$ conda install scikit-learn
```
2. Download the DeepAnnotation source files
```bash
# direct from github
$ git clone https://github.com/mawenlong2016/DeepAnnotation.git     
# direct from website URL
$ wget https://github.com/mawenlong2016/DeepAnnotation/archive/refs/heads/main.zip  
$ unzip main.zip
$ mv DeepAnnotation-main DeepAnnotation
```
3. Prepare your files with the same format that we have prepared in the examples of Duroc (./DeepAnnotation/species/Duroc)
```bash
# Due to the size limitation of github webline
# The Genotype_train_example.h5 was split into a couple of files with less than 25MB
# Therefore, we should combine these slices before we run DeepAnnotation example
$ cat ./DeepAnnotation/species/Duroc/Genotype_train_example_slice* > ./DeepAnnotation/species/Duroc/Genotype_train_example.h5
```
4. activate DeepAnnotation conda environment and enter the DeepAnnotation directory
```bash
$ conda activate DeepAnnotation
$ cd ./DeepAnnotation
```
### Installed by docker
1. **Prerequisite**: Before proceeding with the installation, ensure that Docker and Nvidia-Docker is installed and running. If you haven't installed Docker yet, you can follow [the official Docker tutorial](https://docs.docker.com/get-docker/) for installation instructions. If you haven't installed Nvidia-Docker yet, you can follow [the official Nvidia Docker Container tutorial](https://docs.nvidia.com/deeplearning/frameworks/user-guide/index.html) for installation instructions.
2. Pull the docker image
```bash
$ docker pull wenlong2023/deepannotation:20231129
```
3. To launch DeepAnnotation, please run the given command in your terminal after performing the previous steps. 
```bash
# launch by docker
$ docker run --runtime=nvidia -it --net=host wenlong2023/deepannotation:20231129 /bin/bash
# launch by nvidia-docker
$ docker run -it --net=host wenlong2023/deepannotation:20231129 /bin/bash
# start conda environment
$ source ~/.bashrc
```
4. activate DeepAnnotation conda environment and enter the DeepAnnotation directory
```bash
$ conda activate DeepAnnotation
$ cd /DeepAnnotation
```
## DeepAnnotation exemplified usage <br>
```bash
# for only genotype data
# train DeepAnnotation with cross-validation strategy
$ python ./code/DeepAnnotation_cv.py \
--training_steps 300 --early_stop 30 --kfold 2 --learning_rate 0.01 \
--genotype ./species/Duroc/Genotype_train_example.h5 \
--phenotype ./species/Duroc/phenotype_train_example.txt \
--codingSNP ./species/Duroc/codingSNP_flag_example.txt  \
--work_path ./species/Duroc/Network_example \
--seed 1 --verbose True --save_prcocess True --save_earlymodel True --phe_flag 5 \
--batch_ratio 0.1 --max_gpu 0.9 --decay_steps 100 --decay_rate 0.9 --seed 1  \
--architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn \
--unit 90 --regularizer l2 --regularizer_rate 0.1 --momentum 0.95 --dropout 0.1 \
--bias_architecture none,none,none \
--unit_ratio 1,2,4,5,6,7 \
--calculate -1
# train DeepAnnotation with specific hyperparameters and make predictions for candidate samples
$ python ./code/DeepAnnotation_train.py \
--training_steps 300  --learning_rate 0.01 \
--genotype ./species/Duroc/Genotype_train_example.h5 \
--phenotype ./species/Duroc/phenotype_train_example.txt \
--codingSNP ./species/Duroc/codingSNP_flag_example.txt  \
--work_path ./species/Duroc/Train \
--architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn \
--unit 90 --regularizer l2 --regularizer_rate 0.1 --momentum 0.95 --dropout 0.1 \
--bias_architecture none,none,none \
--unit_ratio 1,2,4,5,6,7 \
--calculate -1 --phe_flag 5 --verbose True \
--predict True --batch_ratio 0.1 --max_gpu 0.9 --decay_steps 100 --decay_rate 0.9 \
--genotype_test ./species/Duroc/Genotype_test_example.h5
# make predictions for candidate samples with trained DeepAnnotation model
$ python ./code/DeepAnnotation_predict.py \
--training_steps 300  --learning_rate 0.01 \
--genotype ./species/Duroc/Genotype_train_example.h5 \
--phenotype ./species/Duroc/phenotype_train_example.txt \
--codingSNP ./species/Duroc/codingSNP_flag_example.txt \
--work_path ./species/Duroc/Train \
--architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn \
--unit 90 --regularizer l2 --regularizer_rate 0.1 --momentum 0.95 --dropout 0.1 \
--bias_architecture none,none,none \
--unit_ratio 1,2,4,5,6,7 \
--calculate -1 --phe_flag 5 --batch_ratio 0.1 --max_gpu 0.9 \
--decay_steps 100 --decay_rate 0.9 --seed 1 \
--genotype_test ./species/Duroc/Genotype_test_example.h5 \
--model_dir ./species/Duroc/Train/final_model

# for both genotype and SNP functional impact data
# train DeepAnnotation with cross-validation strategy
$ python ./code/DeepAnnotation_cv.py \
--training_steps 300 --early_stop 30 --kfold 2 --learning_rate 0.01 \
--genotype ./species/Duroc/Genotype_train_example.h5 \
--phenotype ./species/Duroc/phenotype_train_example.txt \
--codingSNP ./species/Duroc/codingSNP_flag_example.txt \
--work_path ./species/Duroc/Network_example --seed 1 --verbose True \
--save_prcocess True --save_earlymodel True --phe_flag 5 \
--batch_ratio 0.1 --max_gpu 0.9 --decay_steps 100 --decay_rate 0.9 --seed 1 \
--architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn \
--unit 90 --regularizer l2 --regularizer_rate 0.1 --momentum 0.95 --dropout 0.1 \
--bias_architecture fnn,none,none \
--unit_ratio 1,2,4,5,6,7 \
--calculate -1
# train DeepAnnotation with specific hyperparameters and make predictions for candidate samples
$ python ./code/DeepAnnotation_train.py \
--training_steps 300  --learning_rate 0.01 \
--genotype ./species/Duroc/Genotype_train_example.h5 \
--phenotype ./species/Duroc/phenotype_train_example.txt  \
--codingSNP ./species/Duroc/codingSNP_flag_example.txt  \
--work_path ./species/Duroc/Train \
--architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn \
--unit 90 --regularizer l2 --regularizer_rate 0.1 --momentum 0.95 --dropout 0.1 \
--bias_architecture fnn,none,none \
--unit_ratio 1,2,4,5,6,7 \
--calculate -1 --phe_flag 5 --verbose True --predict True \
--batch_ratio 0.1 --max_gpu 0.9 --decay_steps 100 --decay_rate 0.9 \
--genotype_test ./species/Duroc/Genotype_test_example.h5
# make predictions for candidate samples with trained DeepAnnotation model
$ python ./code/DeepAnnotation_predict.py \
--training_steps 300  --learning_rate 0.01 \
--genotype ./species/Duroc/Genotype_train_example.h5 \
--phenotype ./species/Duroc/phenotype_train_example.txt \
--codingSNP ./species/Duroc/codingSNP_flag_example.txt \
--work_path ./species/Duroc/Train \
--architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn \
--unit 90 --regularizer l2 --regularizer_rate 0.1 --momentum 0.95 --dropout 0.1 \
--bias_architecture fnn,none,none \
--unit_ratio 1,2,4,5,6,7 \
--calculate -1 --phe_flag 5 --batch_ratio 0.1 --max_gpu 0.9 \
--decay_steps 100 --decay_rate 0.9 --seed 1 \
--genotype_test ./species/Duroc/Genotype_test_example.h5 \
--model_dir ./species/Duroc/Train/final_model


# for genotype, SNP functional impact, and gene annotation data
# train DeepAnnotation with cross-validation strategy
$ python ./code/DeepAnnotation_cv.py \
--training_steps 300 --early_stop 30 --kfold 2 --learning_rate 0.01 \
--genotype ./species/Duroc/Genotype_train_example.h5 \
--phenotype ./species/Duroc/phenotype_train_example.txt \
--function ./species/Duroc/GeneFunc_example.h5  \
--codingSNP ./species/Duroc/codingSNP_flag_example.txt  \
--work_path ./species/Duroc/Network_example \
--seed 1 --verbose True --save_prcocess True --save_earlymodel True \
--phe_flag 5 --batch_ratio 0.1 --max_gpu 0.9 --decay_steps 100 --decay_rate 0.9 --seed 1  \
--architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn \
--unit 90 --regularizer l2 --regularizer_rate 0.1 --momentum 0.95 --dropout 0.1 \
--bias_architecture fnn,fnn,none \
--unit_ratio 1,2,4,5,6,7 \
--calculate -1
# train DeepAnnotation with specific hyperparameters and make predictions for candidate samples
$ python ./code/DeepAnnotation_train.py \
--training_steps 300  --learning_rate 0.01 \
--genotype ./species/Duroc/Genotype_train_example.h5 \
--phenotype ./species/Duroc/phenotype_train_example.txt \
--function ./species/Duroc/GeneFunc_example.h5 \
--codingSNP ./species/Duroc/codingSNP_flag_example.txt  \
--work_path ./species/Duroc/Train \
--architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn \
--unit 90 --regularizer l2 --regularizer_rate 0.1 --momentum 0.95 --dropout 0.1 \
--bias_architecture fnn,fnn,none \
--unit_ratio 1,2,4,5,6,7 \
--calculate -1 --phe_flag 5 --verbose True --predict True \
--batch_ratio 0.1 --max_gpu 0.9 --decay_steps 100 --decay_rate 0.9 \
--genotype_test ./species/Duroc/Genotype_test_example.h5
# make predictions for candidate samples with trained DeepAnnotation model
$ python ./code/DeepAnnotation_predict.py \
--training_steps 300  --learning_rate 0.01 \
--genotype ./species/Duroc/Genotype_train_example.h5 \
--phenotype ./species/Duroc/phenotype_train_example.txt \
--function ./species/Duroc/GeneFunc_example.h5 \
--codingSNP ./species/Duroc/codingSNP_flag_example.txt \
--work_path ./species/Duroc/Train \
--architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn \
--unit 90 --regularizer l2 --regularizer_rate 0.1 --momentum 0.95 --dropout 0.1 \
--bias_architecture fnn,fnn,none \
--unit_ratio 1,2,4,5,6,7 \
--calculate -1 --phe_flag 5 --batch_ratio 0.1 --max_gpu 0.9 \
--decay_steps 100 --decay_rate 0.9 --seed 1 \
--genotype_test ./species/Duroc/Genotype_test_example.h5 \
--model_dir ./species/Duroc/Train/final_model


# for all genotype, SNP functional impact, gene annotation, and metaterm regulatory module data
# train DeepAnnotation with cross-validation strategy
$ python ./code/DeepAnnotation_cv.py \
--training_steps 300 --early_stop 30 --kfold 2 --learning_rate 0.01 \
--genotype ./species/Duroc/Genotype_train_example.h5 \
--phenotype ./species/Duroc/phenotype_train_example.txt \
--function ./species/Duroc/GeneFunc_example.h5  \
--metagene ./species/Duroc/metagene_example.h5 \
--codingSNP ./species/Duroc/codingSNP_flag_example.txt \
--work_path ./species/Duroc/Network_example \
--seed 1 --verbose True --save_prcocess True --save_earlymodel True \
--phe_flag 5 --batch_ratio 0.1 --max_gpu 0.9 --decay_steps 100 --decay_rate 0.9 --seed 1  \
--architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn \
--unit 90 --regularizer l2 --regularizer_rate 0.1 --momentum 0.95 --dropout 0.1 \
--bias_architecture fnn,fnn,fnn \
--unit_ratio 1,2,4,5,6,7 \
--calculate -1
# train DeepAnnotation with specific hyperparameters and make predictions for candidate samples
$ python ./code/DeepAnnotation_train.py \
--training_steps 300  --learning_rate 0.01 \
--genotype ./species/Duroc/Genotype_train_example.h5 \
--phenotype ./species/Duroc/phenotype_train_example.txt \
--function ./species/Duroc/GeneFunc_example.h5 \
--metagene ./species/Duroc/metagene_example.h5 \
--codingSNP ./species/Duroc/codingSNP_flag_example.txt  \
--work_path ./species/Duroc/Train \
--architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn \
--unit 90 --regularizer l2 --regularizer_rate 0.1 --momentum 0.95 --dropout 0.1 \
--bias_architecture fnn,fnn,fnn \
--unit_ratio 1,2,4,5,6,7 \
--calculate -1 --phe_flag 5 --verbose True --predict True --batch_ratio 0.1 \
--max_gpu 0.9 --decay_steps 100 --decay_rate 0.9 \
--genotype_test ./species/Duroc/Genotype_test_example.h5
# make predictions for candidate samples with trained DeepAnnotation model
$ python ./code/DeepAnnotation_predict.py \
--training_steps 300  --learning_rate 0.01 \
--genotype ./species/Duroc/Genotype_train_example.h5 \
--phenotype ./species/Duroc/phenotype_train_example.txt \
--function ./species/Duroc/GeneFunc_example.h5 \
--metagene ./species/Duroc/metagene_example.h5 \
--codingSNP ./species/Duroc/codingSNP_flag_example.txt \
--work_path ./species/Duroc/Train \
--architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn \
--unit 90 --regularizer l2 --regularizer_rate 0.1 --momentum 0.95 --dropout 0.1 \
--bias_architecture fnn,fnn,fnn \
--unit_ratio 1,2,4,5,6,7 \
--calculate -1 --phe_flag 5 --batch_ratio 0.1 --max_gpu 0.9 \
--decay_steps 100 --decay_rate 0.9 --seed 1 \
--genotype_test ./species/Duroc/Genotype_test_example.h5 \
--model_dir ./species/Duroc/Train/final_model
```
## Ask questions
Please use [DeepAnnotation/issues](https://github.com/mawenlong2016/DeepAnnotation/issues) for how to use DeepAnnotation and reporting bugs.
