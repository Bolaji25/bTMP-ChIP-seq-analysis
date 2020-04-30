#! /usr/bin/bash

## Allocate resources
#SBATCH --time=1-20:00:00
#SBATCH --array=4

#SBATCH --mail-user=bolaji.isiaka@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="alignment"
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --mem-per-cpu=20G
#SBATCH --tmp=32G

module add vital-it;

module load R/3.6.1;

Rscript ./quasralign.R
