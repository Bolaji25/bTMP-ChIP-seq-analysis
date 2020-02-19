#! /usr/bin/bash

## Allocate resources
#SBATCH --time=1-20:00:00
#SBATCH --array=1

#SBATCH --mail-user=bolaji.isiaka@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="mappabilityandBEADS"
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --mem-per-cpu=20G
#SBATCH --tmp=32G

module add vital-it;

module load R/3.6.1;
module add UHTS/Analysis/bamtools/2.4.1;
module add UHTS/Analysis/RGT/0.11.4;

gemDIRold=/home/ubelix/izb/bi18k694/PhD_analysis/supercoiling/analysis/scripts/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin
gemDIR=/home/ubelix/izb/bi18k694/PhD_analysis/supercoiling/analysis/scripts/gem3-mapper/bin
# create a folder to store all results
basefolder="/home/ubelix/izb/bi18k694/PhD_analysis/supercoiling/analysis/scripts/gemresults"
mkdir -p ${basefolder} && cd ${basefolder}
 
pref="ws270"
reference="./c_elegans.PRJNA13758.WS270.genomic.fa"
idxpref="ws270_index"
#thr=8; # use 8 cores
 
${gemDIR}/gem-indexer -i ${reference} -o ${idxpref} 
${gemDIR}/gem-indexer  -i ${reference} -o ${idxpref}
${gemDIRold}/gem-indexer  -i ${reference} -o ${idxpref}
${gemDIRold}/gem-indexer  -i ${reference} -o ${idxpref}
