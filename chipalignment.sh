#! /usr/bin/bash

## Allocate resources
#SBATCH --time=1-24:00:00
#SBATCH --array=1-8

#SBATCH --mail-user=bolaji.isiaka@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="alignentbwaprac"
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --mem-per-cpu=20G
#SBATCH --tmp=32G

module add vital-it;
module add UHTS/Quality_control/fastqc/0.11.7;
#module load gcc/6.2.0  
module add UHTS/Aligner/bowtie2/2.3.1;
module add UHTS/Analysis/samtools/1.10;
module add UHTS/Analysis/sambamba/0.7.1;

# get index for list of samples
#let i=$SLURM_ARRAY_TASK_ID-1

# This script takes a fastq file of ChIP-seq data, runs FastQC and outputs a BAM file for it that is ready for peak calling. Bowtie2 is the aligner used, and the outputted BAM file is sorted by genomic coordinates and has multi-mappers and duplicate reads removed using sambamba.
# USAGE: sh chipseq_analysis_on_input_file.sh <path to the fastq file>

# initialize a variable with an intuitive name to store the name of the input fastq file
base=$1
fq=(`ls /home/ubelix/izb/bi18k694/PhD_analysis/supercoiling/analysis/bTMPseq/trimmed2020110/${base}*.fq`) 

# grab base of filename for naming outputs
base=`basename $fq .fq`
echo "Sample name is $base"    

# directory with bowtie genome index
genome=./reference/c.elegans_270   

# make all of the output directories
# The -p option means mkdir will create the whole path if it 
# does not exist and refrain from complaining if it does exist
#mkdir -p ./reseq/fastqc
#mkdir -p ./reseq/intermediate_bams
#
## set up output filenames and locations
#fastqc_out=./reseq/fastqc/
#
### set up file names
align_out=./trimmed2020110/sam/${base}_unsorted.sam
align_bam=./trimmed2020110/bam/${base}_unsorted.bam
align_sorted=./trimmed2020110/bam/${base}_sorted.bam
align_filtered=./trimmed2020110/bam/${base}_aln.bam
align_flag=./trimmed2020110/${base}_aln.txt
align_flag2=./trimmed2020110/bam/text/${base}_aln2.txt
#
### set up more variables for 2 additional directoties to help clean up the results folder
#bowtie_results=./reseq
#intermediate_bams=./reseq/intermediate_bams
#
echo "Processing file $fq"
#
## Run FastQC and move output to the appropriate folder
#fastqc $fq
#
## Run bowtie2
bowtie2 -p 8 -q --local -x $genome -U $fq -S $align_out
#
## Create BAM from SAM
samtools view -h -S -b -@ 8 -o $align_bam $align_out
#
##check stat 
samtools flagstat  $align_bam > $align_flag
#
## Sort BAM file by genomic coordinates
samtools sort -@ 8 -o $align_sorted $align_bam

samtools flagstat $align_sorted > $align_flag2
#
## Filter out duplicates
##sambamba view -h -t 8 -f bam -F "[XS] == null and not unmapped " $align_bam > $align_filtered
#
## Create indices for all the bam files for visualization and QC
##samtools index $align_filtered
#
## Move intermediate files out of the bowtie2 directory
#mv $bowtie_results/${base}*sorted* $intermediate_bams
