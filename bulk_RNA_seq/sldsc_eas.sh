#!/bin/bash
#$ -S /bin/sh

dd=${1}
jobid=${SGE_TASK_ID}
export PATH=/home/imgishi/miniconda3/envs/ldsc/bin:${PATH}



data=$( cat -n $dd/info/data_list |
   awk -v jobid=$jobid '{if($1==jobid){print $2}}' )

 #suppress BLAS to use multiple CPUs
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

 #output directory

odir=${dd}/sldsc/EAS/${data}
mkdir -p $odir

 #main job
ldsc=/home/imgtaka/LDSC/ldsc/ldsc.py
weights2=/home/ha7477/reference/LDSCORE/1000G_Phase3_EAS_weights_hm3_no_MHC/weights.EAS.hm3_noMHC.
frqfile=/home/ha7477/reference/LDSCORE/1000G_Phase3_EAS_frq/1000G.EAS.QC.
baseline=/home/ha7477/reference/LDSCORE/baseline_eas_v1.2/baseline.
mainannot=${dd}/ldscore/EAS/${data}/$data.

for trait in Lupus_ard2021 RA_ishigaki ;do
   
   sumstats=/home/ha7477/reference/gwas_sumstats_ldsc/EAS/$trait/$trait.sumstats.gz
   
   python  $ldsc \
      --h2    $sumstats  \
      --ref-ld-chr  $mainannot,$baseline \
      --w-ld-chr    $weights2  \
      --overlap-annot         \
      --frqfile-chr   $frqfile  \
      --out  $odir/$trait \
      --print-coefficients
   
done
