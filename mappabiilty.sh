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

gemDIR=/home/ubelix/izb/bi18k694/PhD_analysis/supercoiling/analysis/scripts/gem3-mapper/bin
gemDIRold=/home/ubelix/izb/bi18k694/PhD_analysis/supercoiling/analysis/scripts/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin

#git clone --recursive https://github.com/smarco/gem3-mapper.git gem3-mapper
cd /home/ubelix/izb/bi18k694/PhD_analysis/supercoiling/analysis/scripts/gem3-mapper-3.3.0 
./configure
make

cd /home/ubelix/izb/bi18k694/PhD_analysis/supercoiling/analysis/scripts/

gemDIR3.3=/home/ubelix/izb/bi18k694/PhD_analysis/supercoiling/analysis/scripts/gem3-mapper/bin

#script to make mappabilty track for 50bp read length
${gemDIRold}/gem-indexer -i /home/ubelix/izb/bi18k694/genomeversion/ws270/c_elegans.PRJNA13758.WS270.genomic.fa -o ws270indexed 
${gemDIRold}/gem-mappability -I ws270.gem -l 50 -o ws270 
${gemDIRold}/gem-2-wig -I ws270.gem -i ws270.map -o ws270
wigToBigWig ws270.wig ws270.chrom.sizes ws270.bw

#wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
#cd ~
#mkdir tmp
#cd tmp
#https://repo.continuum.io/archive/Anaconda3<release>.sh
#bash Anaconda3-5.2.0-Linux-x86_64.sh

#source ${HOME}/.bashrc



# to be able to activate environments from inside a slurm script you need to add
# the path to the activate script to .bashrc:
#export CONDA_ACTIVATE=/home/ubelix/izb/bi18k694/anaconda3/bin/activate
# then in the script use
#source $CONDA_ACTIVATE
#conda install -c bioconda ucsc-wigtobigwig
#conda update ucsc-wigtobigwig
#conda activate ucsc-wigtobigwig

#wigToBigWig ws270.wig ws270.chrom.sizes ws270.bw 

Rscript ./BEADS.R
