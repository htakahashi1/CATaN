#!/bin/bash
#$ -S /bin/bash


CC=${1}
data_path=${2}


# ---- Set the following paths to match your environment ----
export PATH=$HOME/miniconda3/envs/renv/bin:${PATH}

cd ${data_path}/TFpeak_noexp/observed/${CC}
mkdir -p ${data_path}/CC_10percent_SNP_noexp_2



nrow_sum=`zcat all_target_all_SNP_sort.bed.gz | wc -l`
nrow_3=`echo $((nrow_sum/10))`
zcat all_target_all_SNP_sort.bed.gz |  sort -g -k 4,4 | tail -n ${nrow_3} | bgzip -c > ${data_path}/CC_10percent_SNP_noexp_2/observed_woL2_${CC}_top_sum_SNP.bed.gz
zcat all_target_all_SNP_sort.bed.gz |  sort -g -k 4,4 | head -n ${nrow_3} | bgzip -c > ${data_path}/CC_10percent_SNP_noexp_2/observed_woL2_${CC}_bottom_sum_SNP.bed.gz