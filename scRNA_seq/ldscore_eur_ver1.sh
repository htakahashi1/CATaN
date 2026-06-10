#!/bin/sh
jobid=${SGE_TASK_ID}
export PATH=/home/imgishi/miniconda3/envs/ldsc/bin:${PATH}

dd=${1}


data=$( cat -n $dd/info/data_list |
   awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
   
 #suppress BLAS to use multiple CPUs
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

BEDFILE=$dd/${data}.bed.gz



 #output directory

odir=${dd}/ldscore/EUR/${data}
mkdir -p $odir
cd $odir

 #main job
for chr in $(seq 1 22);do
   
   #input bed file
   zcat $BEDFILE |
   awk -v chr=$chr 'BEGIN{OFS="\t"}{ if($1=="chr"chr){print $1,$2,$3}}'  |
   sort -k 1,1 -k 2,2n > $odir/tmp.$chr.bed
   
   #annotation
   /home/imgtaka/LDSC/make_annot.py  \
      --bimfile  /home/ha7477/reference/LDSCORE/1000G_Phase3_EUR_plinkfiles/1000G.EUR.QC.$chr.bim \
      --bed-file  $odir/tmp.$chr.bed  \
      --annot-file $odir/$data.$chr.annot.gz
   
   #LD score
   /home/imgtaka/LDSC/ldsc/ldsc.py \
      --l2  \
      --thin-annot  \
      --bfile  /home/ha7477/reference/LDSCORE/1000G_Phase3_EUR_plinkfiles/1000G.EUR.QC.$chr   \
      --ld-wind-cm  1  \
      --annot  $odir/$data.$chr.annot.gz  \
      --out   $odir/$data.$chr  \
      --print-snps  /home/ha7477/reference/LDSCORE/baseline_v1.2/$chr.snps
   
   rm -f tmp.$chr.bed
   
done
