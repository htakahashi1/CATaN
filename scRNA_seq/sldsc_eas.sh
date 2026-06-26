#!/bin/bash
#$ -S /bin/sh

jobid=${SGE_TASK_ID}
dd=${1}
LDSCORE_path=${2}
GWAS_path=${3}

# ---- Set the following paths to match your environment ----
export PATH=$HOME/miniconda3/envs/ldsc/bin:${PATH}

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
ldsc=$HOME/LDSC/ldsc/ldsc.py
weights2=${LDSCORE_path}/1000G_Phase3_EAS_weights_hm3_no_MHC/weights.EAS.hm3_noMHC.
frqfile=${LDSCORE_path}/1000G_Phase3_EAS_frq/1000G.EAS.QC.
baseline=${LDSCORE_path}/baseline_eas_v1.2/baseline.
mainannot=${dd}/ldscore/EAS/${data}/$data.

trait = "RA_ishigaki"  
sumstats=${GWAS_path}/EAS/$trait/$trait.sumstats.gz
python  $ldsc \
   --h2    $sumstats  \
   --ref-ld-chr  $mainannot,$baseline \
   --w-ld-chr    $weights2  \
   --overlap-annot         \
   --frqfile-chr   $frqfile  \
   --out  $odir/$trait \
   --print-coefficients
   

