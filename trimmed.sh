#! /usr/bin/bash

## Allocate resources
#SBATCH --time=1-12:00:00
#SBATCH --array=1

#SBATCH --mail-user=bolaji.isiaka@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="trimming"
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --mem-per-cpu=20G
#SBATCH --tmp=32G

module add vital-it;
module add UHTS/Quality_control/cutadapt/2.5;

module load R/3.5.1;
module add UHTS/Quality_control/fastqc/0.11.7;      
module add UHTS/Analysis/trimmomatic/0.36;
module add UHTS/Analysis/samtools/1.8;
module add UHTS/Analysis/picard-tools/2.18.11;
module add UHTS/Quality_control/qualimap/2.2.1;
module add UHTS/Aligner/bwa/0.7.17;
module add UHTS/Analysis/bamtools/2.4.1;



# Install Trim Galore
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.5.tar.gz -o trim_galore.tar.gz
tar xvzf trim_galore.tar.gz
#PATH=$PATH: /home/ubelix/izb/bi18k694/PhD_analysis/supercoiling/analysis/scripts/TrimGalore-0.6.5
#export PATH
# Run Trim Galore
#./TrimGalore-0.6.5/trim_galore

#mkdir ./trimmedreseq/
#mkdir ./qcreseq/

#gunzip ../rawData/reseq2020110/*.fastq.gz
#fastqc /home/ubelix/izb/bi18k694/PhD_analysis/supercoiling/analysis/rawData/reseq2020110/*.fastq -o ./qcreseq/
#the clip is to clip off 50bp from the 100bp reads to make it similar to the previous, omit this option for 50bp read length reads
./TrimGalore-0.6.5/trim_galore -q 20 --three_prime_clip_R1 50  /home/ubelix/izb/bi18k694/PhD_analysis/supercoiling/analysis/rawData/reseq2020110/*.fastq -o ./trimmedreseq/
fastqc ./trimmedreseq/*.fq -o ./qcreseq/




