#! /usr/bin/bash

## Allocate resources
#SBATCH --time=1-24:00:00
#SBATCH --array=1

#SBATCH --mail-user=bolaji.isiaka@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="remdupl"
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --mem-per-cpu=20G
#SBATCH --tmp=32G

#module add vital-it;
#module add Development/java/11.0.6;
#module add UHTS/Analysis/samtools/1.10;

base=$1
fq=(`ls /home/ubelix/izb/bi18k694/PhD_analysis/supercoiling/analysis/bTMPseq/reseq/${base}*.bam`)
base=`basename $fq .bam`
echo "Sample name is $base"

align_sorted=./reseq/${base}_sorted.bam
align_duprem=./reseq/${base}_duprem.bam
align_qc=./reseq/${base}_Duprem.txt
align_flag=./reseq/${base}_dupflag.txt

mkdir -p ./reseq/intermediate_dups
intermediate_dups=./reseq/intermediate_dups

#${picardDIR} = /home/ubelix/izb/bi18k694/PhD_analysis/supercoiling/analysis/bTMPseq/

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
