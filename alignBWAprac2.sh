#!  /usr/bin/bash

## Allocate resources
#SBATCH --time=1-40:00:00
#SBATCH --array=1-24

#SBATCH --mail-user=bolaji.isiaka@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="alignentbwa"
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --mem-per-cpu=20G
#SBATCH --tmp=32G

module add vital-it;
module add UHTS/Quality_control/cutadapt/2.5;

module load R/3.6.1;
module add UHTS/Quality_control/fastqc/0.11.7;
module add UHTS/Analysis/trimmomatic/0.36;
module add UHTS/Analysis/samtools/1.8;
module add UHTS/Analysis/picard-tools/2.18.11;
module add UHTS/Quality_control/qualimap/2.2.1;
module add UHTS/Aligner/bwa/0.7.17;
module add UHTS/Analysis/bamtools/2.4.1;

fastqc

###############################
########### VARIABLES #########
###############################
# variable are set in the varSettings.sh file and values read in from command l$

#source ./variablesettings.sh

#the start of the basename of the fastq files
bname=$1

#the biological test group that you will later compare (used for preparing the $
#testGroup=$2

# number of threads
numThreads=$3
# get foward and reverse read files for this sample
fileList=( `ls ./trimmed/${bname}*.fq` )


#######################################################
## get read stats                            ##
#######################################################

#run fastqc on sequences
mkdir -p ./alignedBWAprac2/fastQC/

fastqc -t ${numThreads} ${fileList[@]} -o ./alignedBWAprac2/fastQC/

# Use precompiled binaries (recommended)
curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.0pre2/bwa-mem$
  | tar jxf -
bwa-mem2-2.0pre2_x64-linux/bwa-mem2 index /home/ubelix/izb/bi18k694/genomeversion/ws270/c_elegans.PRJNA13758.WS270.genomic.fa
bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem /home/ubelix/izb/bi18k694/genomeversion/ws270/c_elegans.PRJNA13758.WS270.genomic.fa ${fileList[@]} > ./alignedBWAprac2/

