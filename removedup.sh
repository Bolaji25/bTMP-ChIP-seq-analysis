#! /usr/bin/bash


base=$1
fq=(`ls /home/ubelix/izb/bi18k694/PhD_analysis/supercoiling/analysis/bTMPseq/trimmed2020110/bam/${base}*.bam`)
base=`basename $fq .bam`
echo "Sample name is $base"

align_sorted=./trimmed2020110/bam/${base}_sorted.bam
align_duprem=./trimmed2020110/bam/${base}.bam
align_qc=./trimmed2020110/bam/${base}_Duprem.txt
align_flag=./trimmed2020110/bam/${base}_dupflag.txt
align_sortge=./trimmed2020110/bam/${base}_cordsort.bam
align_flag2=./trimmed2020110/bam/${base}_dupflag2.txt

#mkdir -p ./reseq/intermediate_dups
intermediate_dups=./reseq/intermediate_dups

${picardDIR} = /home/ubelix/izb/bi18k694/PhD_analysis/supercoiling/analysis/bTMPseq/

#sort by query name for duplicate removal
java -Xms1g -Xmx8g -jar /home/ubelix/izb/bi18k694/PhD_analysis/supercoiling/analysis/bTMPseq/picard.jar SortSam I=$fq O=$align_sorted SORT_ORDER=queryname TMP_DIR=$intermediate_dups


# remove duplicates with picard (path to picard should be set in $PICARD variable in .bash_profile or in session)
# Note that to mark unmapped mates of mapped records and supplementary/secondary alignments as duplicates the bam
# file must be querysorted (by name) not by coordinate. 
java -Xms1g -Xmx8g -jar /home/ubelix/izb/bi18k694/PhD_analysis/supercoiling/analysis/bTMPseq/picard.jar MarkDuplicates I=$align_sorted O=$align_duprem M=$align_qc REMOVE_DUPLICATES=true REMOVE_SEQUENCING_DUPLICATES=true  
ASSUME_SORT_ORDER=queryname TMP_DIR=$intermediate_dups



rm $align_sorted

#get stats

samtools flagstat  $align_duprem > $align_flag

#sort by genomic coordinate
samtools sort -o $align_sortge $align_duprem
samtools flagstat  $align_sortge > $align_flag2

#index bam files
samtools index $align_sortge

