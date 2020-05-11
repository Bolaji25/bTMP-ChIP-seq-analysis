#! /usr/bin/bash

## Allocate resources
#SBATCH --time=1-5:00:00
#SBATCH --array=1

#SBATCH --mail-user=bolaji.isiaka@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="normalization"
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --mem-per-cpu=20G
#SBATCH --tmp=32G

module add vital-it;

module load R/3.6.1;
module add UHTS/Analysis/bamtools/2.4.1;

Rscript ./Normalization+8.R
