#!/bin/sh
jobid=${SGE_TASK_ID}
dd=${1}
LDSCORE_path=${2}

# ---- Set the following paths to match your environment ----
export PATH=$HOME/miniconda3/envs/ldsc/bin:${PATH}

data=$( cat -n $dd/info/data_list |
   awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
   
 #suppress BLAS to use multiple CPUs
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

BEDFILE=$dd/${data}.bed.gz



 #output directory

odir=${dd}/ldscore/EAS/${data}
mkdir -p $odir
cd $odir

 #main job
for chr in $(seq 1 22);do
   
   #input bed file
   zcat $BEDFILE |
   awk -v chr=$chr 'BEGIN{OFS="\t"}{ if($1=="chr"chr){print $1,$2,$3}}'  |
   sort -k 1,1 -k 2,2n > $odir/tmp.$chr.bed
   
   #annotation
   $HOME/LDSC/make_annot.py  \
      --bimfile  ${LDSCORE_path}/1000G_Phase3_EAS_plinkfiles/1000G.EAS.QC.$chr.bim  \
      --bed-file  $odir/tmp.$chr.bed  \
      --annot-file $odir/$data.$chr.annot.gz
   
   #LD score
   $HOME/LDSC/ldsc/ldsc.py \
      --l2  \
      --thin-annot  \
      --bfile  ${LDSCORE_path}/1000G_Phase3_EAS_plinkfiles/1000G.EAS.QC.$chr   \
      --ld-wind-cm  1  \
      --annot  $odir/$data.$chr.annot.gz  \
      --out   $odir/$data.$chr  \
      --print-snps  ${LDSCORE_path}/baseline_eas_v1.2/$chr.snps
   
   rm -f tmp.$chr.bed
   
done
