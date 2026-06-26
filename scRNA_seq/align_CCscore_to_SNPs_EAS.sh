#!/bin/bash
#$ -S /bin/bash


CC=${1}
data_path=${2}
save_path=${3}
TF_ChIP_path=${4}
1000G_EAS_plink_path=${5}


# ---- Set the following paths to match your environment ----
export PATH=$HOME/miniconda3/envs/renv/bin:${PATH}

ODIR=${save_path}/TFpeak_noexp_EAS/observed/${CC}
mkdir -p ${ODIR}

dd=${TF_ChIP_path}

data_list=`cat ${data_path}/result/raw_CC/${CC}_CCscore.txt | awk '{print $1}'`

for chr in `seq 22`
do
rm -f $ODIR/${chr}_tmp.all.bed.gz 
for data in $data_list
do
   
BEDFILE=$dd/$data.hg19.bed

score=$( cat ${data_path}/result/raw_CC/${CC}_CCscore.txt |
      awk -v data=$data '{if($1==data){print $2}}' )
  
    cat $BEDFILE |
   awk -v chr=$chr -v score=$score 'BEGIN{OFS="\t"}{
      if($1=="chr"chr){print $1,$2,$3,score}
   }'  | bgzip -c >> $ODIR/${chr}_tmp.all.bed.gz
   done

 #sort for bedtools
zcat $ODIR/${chr}_tmp.all.bed.gz |
sort -k1,1 -k2,2n - | bgzip -c > $ODIR/${chr}_tmp.all.sorted.bed.gz

#less $ODIR/${chr}_tmp.all.sorted.bed.gz
   
rm $ODIR/${chr}_tmp.all.bed.gz


cat ${1000G_EAS_plink_path}/1000G.EAS.QC.$chr.bim |
awk 'BEGIN{OFS="\t"}{print "chr"$1,$4-1,$4,$2}' |
sort -k1,1 -k2,2n - | bgzip -c > $ODIR/tmp.$chr.target.bed.gz

bedtools intersect \
   -a $ODIR/${chr}_tmp.all.sorted.bed.gz \
   -b $ODIR/tmp.$chr.target.bed.gz |
sort -k1,1 -k2,2n - | bgzip -c  > $ODIR/${chr}_tmp.target_with_score.bed.gz


bedtools merge -d -1 -c 4 -o sum \
   -i $ODIR/${chr}_tmp.target_with_score.bed.gz | sort -k1,1 -k2,2n - | bgzip -c > $ODIR/${chr}_target_with_score_sum.bed.gz

bedtools intersect -v \
   -a $ODIR/tmp.$chr.target.bed.gz \
   -b $ODIR/${chr}_target_with_score_sum.bed.gz | 
sort -k1,1 -k2,2n - | bgzip -c  > $ODIR/${chr}_target_withoutscore.bed.gz

rm -f $ODIR/${chr}_target_withoutscore_fillzero.bed.gz
zcat $ODIR/${chr}_target_withoutscore.bed.gz |
   awk 'BEGIN{OFS="\t"}{print $1,$2,$3,0}'  | bgzip -c >> $ODIR/${chr}_target_withoutscore_fillzero.bed.gz

cat $ODIR/${chr}_target_with_score_sum.bed.gz $ODIR/${chr}_target_withoutscore_fillzero.bed.gz > $ODIR/${chr}_target_all_SNP.bed.gz
zcat $ODIR/${chr}_target_all_SNP.bed.gz | sort -k1,1 -k2,2n - | bgzip -c  > $ODIR/${chr}_target_all_SNP_sort.bed.gz

rm $ODIR/${chr}_tmp.all.sorted.bed.gz
rm $ODIR/${chr}_tmp.target_with_score.bed.gz
rm $ODIR/tmp.$chr.target.bed.gz
rm $ODIR/${chr}_target_withoutscore.bed.gz
rm $ODIR/${chr}_target_with_score_sum.bed.gz
rm $ODIR/${chr}_target_withoutscore_fillzero.bed.gz
rm $ODIR/${chr}_target_all_SNP.bed.gz

done