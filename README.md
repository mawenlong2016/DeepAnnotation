# ___DeepAnnotation: A novel interpretable deep learning-based genomic selection model that integrates comprehensive functional annotations___ <br>

The python package 'DeepAnnotation' can be used to perform genomic selection (GS), which is a promising
breeding strategy for agricultural breeding. DeepAnnotation predicts phenotypes from comprehensive multi-omics functional annotations with interpretable deep learning framework. The effectiveness
of DeepAnnotation has been demonstrated in predicting three pork production traits (lean meat percentage at 100 kg [LMP], loin muscle depth at 100 kg [LMD], back fat thickness at 100 kg [BF]) on a population
of 1940 Duroc boars with 11633164 SNPs (_Sus scrofa_) from [**GigaDB**](https://doi.org/10.5524/100894)). The comprehensive functional annotation data and data used for training DeepAnnotation are available at [**Zenodo**](https://zenodo.org/records/14586911).
<br>
## Version and download <br>
* [Version 1.0](https://github.com/mawenlong2016/DeepAnnotation/archive/refs/heads/main.zip) -First version released on Jan, second, 2024<br>
## DeepAnnotation dependencies
DeepAnnotation denpends on the following software environments：<br>
1. [Python](https://www.python.org) - The python (version >=3.6) is needed. <br>
2. We suggest [conda](https://docs.conda.io/projects/miniconda/en/latest/) or [docker](https://www.docker.com/) to install the required packages. <br>
3. We suggest [git](https://docs.github.com/) or [wget](https://www.gnu.org/software/wget/) to access the DeepAnnotation source files. <br>
## Data accession
The full datasets for training and testing DeepAnnotation are available at [**Zenodo**](https://zenodo.org/records/14586911). Specifically, we provide an example dataset (compressed file [**DeepAnnotation-main.zip**](https://zenodo.org/records/14586911/files/DeepAnnotation-main.zip?download=1) at **Zenodo repository**) to learn how to use DeepAnnotation, which is exemplified in the following **DeepAnnotation exemplified usage** section. Moreover, we also provide the example dataset and package the environment at [**DockerHub**](https://hub.docker.com/r/wenlong2023/deepannotation), which could be accessed following the instruction at **Installed by docker** section.
## DeepAnnotation Quick Installation
### Installed by personal computer
1. Build DeepAnnotation environment by conda
```bash
$ conda create -n DeepAnnotation python=3.6 numpy==1.19.2 tensorflow-gpu==2.2.0
$ conda activate DeepAnnotation
$ conda install matplotlib
$ conda install scikit-learn
$ pip install framework-reproducibility
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
3. Prepare your files with the same format that we have prepared in the examples of Duroc (./DeepAnnotation/species/Duroc), which is available at **Data accession** section.
```bash
# The example dataset could be accessed via:
$ wget -c https://zenodo.org/records/14586911/files/DeepAnnotation-main.zip?download=1
# unzip the compressed file
$ unzip DeepAnnotation-main.zip
# The Genotype_train_example.h5 was split into a couple of files with less than 25MB
# Therefore, we should combine these slices before we run DeepAnnotation example
$ cat ./DeepAnnotation-main/species/Duroc/Genotype_train_example_slice* > ./DeepAnnotation-main/species/Duroc/Genotype_train_example.h5
# Create a new folder under the DeepAnnotation directory and copy the example dataset
$ mkdir -p ./DeepAnnotation/species/Duroc
$ cp -r ./DeepAnnotation-main/species/Duroc/  ./DeepAnnotation/species/Duroc
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
$ docker pull wenlong2023/deepannotation:20240101
```
3. To launch DeepAnnotation, please run the given command in your terminal after performing the previous steps. 
```bash
# launch by docker
$ docker run --runtime=nvidia -it --net=host wenlong2023/deepannotation:20240101 /bin/bash
# launch by nvidia-docker
$ docker run -it --net=host wenlong2023/deepannotation:20240101 /bin/bash
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
--training_steps 382 --early_stop 30 --kfold 2 --learning_rate 0.1 \
--genotype ./species/Duroc/Genotype_train_example.h5 \
--phenotype ./species/Duroc/phenotype_train_example.txt \
--work_path ./species/Duroc/Network_example \
--verbose True --save_prcocess True --save_earlymodel True \
--phe_flag 5 --batch_ratio 0.1 --max_gpu 0.9 --decay_steps 100 --decay_rate 0.9 \
--architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn \
--unit 110 --regularizer l2 --regularizer_rate 0.01 --momentum 0.9 --dropout 0.1 \
--bias_architecture none,none,none \
--unit_ratio 1,2,4,5,6,7 \
--calculate -1
# train DeepAnnotation with specific hyperparameters and make predictions for candidate samples
$ python ./code/DeepAnnotation_train.py \
--training_steps 382  --learning_rate 0.1 \
--genotype ./species/Duroc/Genotype_train_example.h5 \
--phenotype ./species/Duroc/phenotype_train_example.txt \
--work_path ./species/Duroc/Train \
--architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn \
--unit 110 --regularizer l2 --regularizer_rate 0.01 --momentum 0.9 --dropout 0.1 \
--bias_architecture none,none,none \
--unit_ratio 1,2,4,5,6,7 \
--calculate -1 --phe_flag 5 --verbose True \
--predict True --batch_ratio 0.1 --max_gpu 0.9 --decay_steps 100 --decay_rate 0.9 \
--genotype_test ./species/Duroc/Genotype_test_example.h5
# make predictions for candidate samples with trained DeepAnnotation model
$ python ./code/DeepAnnotation_predict.py \
--training_steps 382  --learning_rate 0.1 \
--genotype ./species/Duroc/Genotype_train_example.h5 \
--phenotype ./species/Duroc/phenotype_train_example.txt \
--work_path ./species/Duroc/Train \
--architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn \
--unit 110 --regularizer l2 --regularizer_rate 0.01 --momentum 0.9 --dropout 0.1 \
--bias_architecture none,none,none \
--unit_ratio 1,2,4,5,6,7 \
--calculate -1 --phe_flag 5 --batch_ratio 0.1 --max_gpu 0.9 \
--decay_steps 100 --decay_rate 0.9 \
--genotype_test ./species/Duroc/Genotype_test_example.h5 \
--model_dir ./species/Duroc/Train/final_model

# for both genotype and SNP functional impact data
# train DeepAnnotation with cross-validation strategy
$ python ./code/DeepAnnotation_cv.py \
--training_steps 382 --early_stop 30 --kfold 2 --learning_rate 0.1 \
--genotype ./species/Duroc/Genotype_train_example.h5 \
--phenotype ./species/Duroc/phenotype_train_example.txt \
--work_path ./species/Duroc/Network_example \
--verbose True --save_prcocess True --save_earlymodel True \
--phe_flag 5 --batch_ratio 0.1 --max_gpu 0.9 --decay_steps 100 --decay_rate 0.9 \
--architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn \
--unit 110 --regularizer l2 --regularizer_rate 0.01 --momentum 0.9 --dropout 0.1 \
--bias_architecture fnn,none,none \
--unit_ratio 1,2,4,5,6,7 \
--calculate -1
# train DeepAnnotation with specific hyperparameters and make predictions for candidate samples
$ python ./code/DeepAnnotation_train.py \
--training_steps 382  --learning_rate 0.1 \
--genotype ./species/Duroc/Genotype_train_example.h5 \
--phenotype ./species/Duroc/phenotype_train_example.txt \
--work_path ./species/Duroc/Train \
--architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn \
--unit 110 --regularizer l2 --regularizer_rate 0.01 --momentum 0.9 --dropout 0.1 \
--bias_architecture fnn,none,none \
--unit_ratio 1,2,4,5,6,7 \
--calculate -1 --phe_flag 5 --verbose True \
--predict True --batch_ratio 0.1 --max_gpu 0.9 --decay_steps 100 --decay_rate 0.9 \
--genotype_test ./species/Duroc/Genotype_test_example.h5
# make predictions for candidate samples with trained DeepAnnotation model
$ python ./code/DeepAnnotation_predict.py \
--training_steps 382  --learning_rate 0.1 \
--genotype ./species/Duroc/Genotype_train_example.h5 \
--phenotype ./species/Duroc/phenotype_train_example.txt \
--work_path ./species/Duroc/Train \
--architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn \
--unit 110 --regularizer l2 --regularizer_rate 0.01 --momentum 0.9 --dropout 0.1 \
--bias_architecture fnn,none,none \
--unit_ratio 1,2,4,5,6,7 \
--calculate -1 --phe_flag 5 --batch_ratio 0.1 --max_gpu 0.9 \
--decay_steps 100 --decay_rate 0.9 \
--genotype_test ./species/Duroc/Genotype_test_example.h5 \
--model_dir ./species/Duroc/Train/final_model


# for genotype, SNP functional impact, and gene annotation data
# train DeepAnnotation with cross-validation strategy
$ python ./code/DeepAnnotation_cv.py \
--training_steps 382 --early_stop 30 --kfold 2 --learning_rate 0.1 \
--genotype ./species/Duroc/Genotype_train_example.h5 \
--phenotype ./species/Duroc/phenotype_train_example.txt \
--function ./species/Duroc/GeneFunc_example.h5  \
--work_path ./species/Duroc/Network_example \
--verbose True --save_prcocess True --save_earlymodel True \
--phe_flag 5 --batch_ratio 0.1 --max_gpu 0.9 --decay_steps 100 --decay_rate 0.9 \
--architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn \
--unit 110 --regularizer l2 --regularizer_rate 0.01 --momentum 0.9 --dropout 0.1 \
--bias_architecture fnn,fnn,none \
--unit_ratio 1,2,4,5,6,7 \
--calculate -1
# train DeepAnnotation with specific hyperparameters and make predictions for candidate samples
$ python ./code/DeepAnnotation_train.py \
--training_steps 382  --learning_rate 0.1 \
--genotype ./species/Duroc/Genotype_train_example.h5 \
--phenotype ./species/Duroc/phenotype_train_example.txt \
--function ./species/Duroc/GeneFunc_example.h5  \
--work_path ./species/Duroc/Train \
--architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn \
--unit 110 --regularizer l2 --regularizer_rate 0.01 --momentum 0.9 --dropout 0.1 \
--bias_architecture fnn,fnn,none \
--unit_ratio 1,2,4,5,6,7 \
--calculate -1 --phe_flag 5 --verbose True \
--predict True --batch_ratio 0.1 --max_gpu 0.9 --decay_steps 100 --decay_rate 0.9 \
--genotype_test ./species/Duroc/Genotype_test_example.h5
# make predictions for candidate samples with trained DeepAnnotation model
$ python ./code/DeepAnnotation_predict.py \
--training_steps 382  --learning_rate 0.1 \
--genotype ./species/Duroc/Genotype_train_example.h5 \
--phenotype ./species/Duroc/phenotype_train_example.txt \
--function ./species/Duroc/GeneFunc_example.h5  \
--work_path ./species/Duroc/Train \
--architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn \
--unit 110 --regularizer l2 --regularizer_rate 0.01 --momentum 0.9 --dropout 0.1 \
--bias_architecture fnn,fnn,none \
--unit_ratio 1,2,4,5,6,7 \
--calculate -1 --phe_flag 5 --batch_ratio 0.1 --max_gpu 0.9 \
--decay_steps 100 --decay_rate 0.9 \
--genotype_test ./species/Duroc/Genotype_test_example.h5 \
--model_dir ./species/Duroc/Train/final_model


# for all genotype, SNP functional impact, gene annotation, and metaterm regulatory module data
# train DeepAnnotation with cross-validation strategy
$ python ./code/DeepAnnotation_cv.py \
--training_steps 382 --early_stop 30 --kfold 2 --learning_rate 0.1 \
--genotype ./species/Duroc/Genotype_train_example.h5 \
--phenotype ./species/Duroc/phenotype_train_example.txt \
--function ./species/Duroc/GeneFunc_example.h5  \
--metagene ./species/Duroc/metagene_example.h5 \
--work_path ./species/Duroc/Network_example \
--verbose True --save_prcocess True --save_earlymodel True \
--phe_flag 5 --batch_ratio 0.1 --max_gpu 0.9 --decay_steps 100 --decay_rate 0.9 \
--architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn \
--unit 110 --regularizer l2 --regularizer_rate 0.01 --momentum 0.9 --dropout 0.1 \
--bias_architecture fnn,fnn,fnn \
--unit_ratio 1,2,4,5,6,7 \
--calculate -1
# train DeepAnnotation with specific hyperparameters and make predictions for candidate samples
$ python ./code/DeepAnnotation_train.py \
--training_steps 382  --learning_rate 0.1 \
--genotype ./species/Duroc/Genotype_train_example.h5 \
--phenotype ./species/Duroc/phenotype_train_example.txt \
--function ./species/Duroc/GeneFunc_example.h5  \
--metagene ./species/Duroc/metagene_example.h5 \
--work_path ./species/Duroc/Train \
--architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn \
--unit 110 --regularizer l2 --regularizer_rate 0.01 --momentum 0.9 --dropout 0.1 \
--bias_architecture fnn,fnn,fnn \
--unit_ratio 1,2,4,5,6,7 \
--calculate -1 --phe_flag 5 --verbose True \
--predict True --batch_ratio 0.1 --max_gpu 0.9 --decay_steps 100 --decay_rate 0.9 \
--genotype_test ./species/Duroc/Genotype_test_example.h5
# make predictions for candidate samples with trained DeepAnnotation model
$ python ./code/DeepAnnotation_predict.py \
--training_steps 382  --learning_rate 0.1 \
--genotype ./species/Duroc/Genotype_train_example.h5 \
--phenotype ./species/Duroc/phenotype_train_example.txt \
--function ./species/Duroc/GeneFunc_example.h5  \
--metagene ./species/Duroc/metagene_example.h5 \
--work_path ./species/Duroc/Train \
--architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn \
--unit 110 --regularizer l2 --regularizer_rate 0.01 --momentum 0.9 --dropout 0.1 \
--bias_architecture fnn,fnn,fnn \
--unit_ratio 1,2,4,5,6,7 \
--calculate -1 --phe_flag 5 --batch_ratio 0.1 --max_gpu 0.9 \
--decay_steps 100 --decay_rate 0.9 \
--genotype_test ./species/Duroc/Genotype_test_example.h5 \
--model_dir ./species/Duroc/Train/final_model
```
## Ask questions
Please use [DeepAnnotation/issues](https://github.com/mawenlong2016/DeepAnnotation/issues) for how to use DeepAnnotation and reporting bugs.

## Additional files
We downloaded the publicly available epigenome data from the work of [**Kern et al. (2021)**](https://doi.org/10.1038/s41467-021-22100-8). These files have been deposited in the Gene Expression Omnibus (GEO) and are available under accession [**GSE158414**]( https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158414). We downloaded these files at 13 May 2024 via the following two scripts:<br>
**download_conserved_Yorkshire.sh** - we used this file to download the **conserved regulatory elements of pig**.
**download_Yorkshire_ATAC.sh** - we used this file to download the **ATAC-Seq data of pig**.

## Comprehensive functional annotation
We used RNAfold, DeepSEA and easyMF approaches to build the the comprehensive functional annotation. The workflow for building those annotaions could be accessed via: **Script_for_building_the_comprehensive_functional_annotations.sh**
