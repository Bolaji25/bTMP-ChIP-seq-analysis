#! /usr/bin/bash


module add vital-it;
module add UHTS/Quality_control/fastqc/0.11.7;
#module load gcc/6.2.0  
module add UHTS/Aligner/bowtie2/2.3.1;
module add UHTS/Analysis/samtools/1.10;
module add UHTS/Analysis/sambamba/0.7.1;


# This script takes a fastq file of ChIP-seq data, outputs a BAM file for it that is ready for peak calling. Bowtie2 is the aligner used, and the outputted BAM file is sorted by genomic coordinates 

# make a variable with an intuitive name to store the name of the input fastq file
base=$1
fq=(`ls /home/ubelix/izb/bi18k694/PhD_analysis/supercoiling/analysis/bTMPseq/trimmed2020110/${base}*.fq`) 

# grab base of filename for naming outputs
base=`basename $fq .fq`
echo "Sample name is $base"    

# directory with bowtie genome index
genome=./reference/c.elegans_270   

# make output directories
# The -p option means mkdir will create the whole path if it does not exist and refrain from complaining if it does exist
mkdir -p ./reseq/fastqc
mkdir -p ./reseq/intermediate_bams

# set up file names
align_out=./trimmed2020110/sam/${base}_unsorted.sam
align_bam=./trimmed2020110/bam/${base}_unsorted.bam
align_sorted=./trimmed2020110/bam/${base}_sorted.bam
align_filtered=./trimmed2020110/bam/${base}_aln.bam
align_flag=./trimmed2020110/${base}_aln.txt
align_flag2=./trimmed2020110/bam/text/${base}_aln2.txt

echo "Processing file $fq"


# Run bowtie2
bowtie2 -p 8 -q --local -x $genome -U $fq -S $align_out

# Create BAM from SAM
samtools view -h -S -b -@ 8 -o $align_bam $align_out

#check stat 
samtools flagstat  $align_bam > $align_flag

# Sort BAM file by genomic coordinates
samtools sort -@ 8 -o $align_sorted $align_bam

samtools flagstat $align_sorted > $align_flag2
