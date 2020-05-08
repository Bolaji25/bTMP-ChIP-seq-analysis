#! /usr/bin/bash

## Allocate resources
#SBATCH --time=0-5:00:00
#SBATCH --array=2,6,8

#SBATCH --mail-user=bolaji.isiaka@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="alignment"
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --mem-per-cpu=8G
#SBATCH --tmp=32G

module add vital-it;
module add UHTS/Quality_control/fastqc/0.11.7;
#module load gcc/6.2.0
module add UHTS/Aligner/bowtie2/2.3.1;
module add UHTS/Analysis/samtools/1.10;
module add UHTS/Analysis/sambamba/0.7.1;


# read in the run specific settings
sampleName=( bFlp1IP bFlp1ipt bFlp2IP bFlp2ipt bN2not1IP bN2no1Inpt bN2not2IP bN2no2Inpt bTrp1input bTrp1IP)

# get index for list of samples
let i=$SLURM_ARRAY_TASK_ID-1


source ./chipalignment.sh ${sampleName[$i]} $SLURM_CPUS_PER_TASK
        echo "running"

