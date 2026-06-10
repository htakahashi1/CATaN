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

odir=${dd}/sldsc/EUR/${data}
mkdir -p $odir


ldsc=/home/imgtaka/LDSC/ldsc/ldsc.py
weights2=/home/ha7477/reference/LDSCORE/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
frqfile=/home/ha7477/reference/LDSCORE/1000G_Phase3_EUR_frq/1000G.EUR.QC.
baseline=/home/ha7477/reference/LDSCORE/baseline_v1.2/baseline.

mainannot=${dd}/ldscore/EUR/${data}/$data.

 #1, EUR-RA Ishigaki Nat Genet 2022
trait=RA_ishigaki
sumstats=/home/ha7477/reference/gwas_sumstats_ldsc/EUR/RA_ishigaki/RA_ishigaki.sumstats.gz

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
sumstats=/home/ha7477/reference/gwas_sumstats_ldsc/EUR/Lupus_langefeld/Lupus_langefeld.sumstats.gz

python  $ldsc \
   --h2    $sumstats  \
   --ref-ld-chr  $mainannot,$baseline \
   --w-ld-chr    $weights2  \
   --overlap-annot    \
   --frqfile-chr   $frqfile  \
   --out  $odir/$trait \
   --print-coefficients
   
        #3, Alkes lab
LIST=$( ls /home/ha7477/reference/gwas_sumstats_ldsc/EUR/ALKES_LAB |
   grep ^PASS_ |
   sed -e "s/.sumstats.gz//g" )

for trait in $LIST;do
   
   sumstats=/home/ha7477/reference/gwas_sumstats_ldsc/EUR/ALKES_LAB/$trait.sumstats.gz
   
   python  $ldsc \
      --h2    $sumstats  \
      --ref-ld-chr  $mainannot,$baseline \
      --w-ld-chr    $weights2  \
      --overlap-annot   \
      --frqfile-chr   $frqfile  \
      --out  $odir/$trait \
      --print-coefficients
   
done