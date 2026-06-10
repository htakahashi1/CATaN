#!/bin/bash
#$ -S /bin/bash


CC=${1}
data_path=${2}


export PATH=/home/imgtaka/miniconda3/envs/renv/bin:${PATH}
export PATH=/home/ha7477/tools/miniconda3/envs/atac1b/bin:${PATH}

cd ${data_path}/TFpeak_noexp_EAS/observed/${CC}
mkdir -p ${data_path}/CC_10percent_SNP_noexp_2_EAS



nrow_sum=`zcat all_target_all_SNP_sort.bed.gz | wc -l`
nrow_3=`echo $((nrow_sum/10))`
zcat all_target_all_SNP_sort.bed.gz |  sort -g -k 4,4 | tail -n ${nrow_3} | bgzip -c > ${data_path}/CC_10percent_SNP_noexp_2_EAS/observed_woL2_${CC}_top_sum_SNP.bed.gz
zcat all_target_all_SNP_sort.bed.gz |  sort -g -k 4,4 | head -n ${nrow_3} | bgzip -c > ${data_path}/CC_10percent_SNP_noexp_2_EAS/observed_woL2_${CC}_bottom_sum_SNP.bed.gz