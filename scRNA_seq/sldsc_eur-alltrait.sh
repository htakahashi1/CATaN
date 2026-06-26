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

export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

 #output directory

odir=${dd}/sldsc/EUR/${data}
mkdir -p $odir


ldsc=$HOME/LDSC/ldsc/ldsc.py
weights2=${LDSCORE_path}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
frqfile=${LDSCORE_path}/1000G_Phase3_EUR_frq/1000G.EUR.QC.
baseline=${LDSCORE_path}/baseline_v1.2/baseline.

mainannot=${dd}/ldscore/EUR/${data}/$data.

 #1, EUR-RA Ishigaki Nat Genet 2022
trait=RA_ishigaki
sumstats=${GWAS_path}/EUR/RA_ishigaki/RA_ishigaki.sumstats.gz

python  $ldsc \
   --h2    $sumstats  \
   --ref-ld-chr  $mainannot,$baseline \
   --w-ld-chr    $weights2  \
   --overlap-annot    \
   --frqfile-chr   $frqfile  \
   --out  $odir/$trait \
   --print-coefficients

 #2, EUR-SLE  Lupus_langefeld
trait=Lupus_langefeld
sumstats=${GWAS_path}/EUR/Lupus_langefeld/Lupus_langefeld.sumstats.gz

python  $ldsc \
   --h2    $sumstats  \
   --ref-ld-chr  $mainannot,$baseline \
   --w-ld-chr    $weights2  \
   --overlap-annot    \
   --frqfile-chr   $frqfile  \
   --out  $odir/$trait \
   --print-coefficients
   
        #3, Alkes lab
LIST=$( ls ${GWAS_path}/EUR/ALKES_LAB |
   grep ^PASS_ |
   sed -e "s/.sumstats.gz//g" )

for trait in $LIST;do
   
   sumstats=${GWAS_path}/EUR/ALKES_LAB/$trait.sumstats.gz
   
   python  $ldsc \
      --h2    $sumstats  \
      --ref-ld-chr  $mainannot,$baseline \
      --w-ld-chr    $weights2  \
      --overlap-annot   \
      --frqfile-chr   $frqfile  \
      --out  $odir/$trait \
      --print-coefficients
   
done