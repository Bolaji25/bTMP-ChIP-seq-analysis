#! /usr/bin/bash

## Allocate resources
#SBATCH --time=1-20:00:00
#SBATCH --array=1-26%6

#SBATCH --mail-user=bolaji.isiaka@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="subsampling"
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --mem-per-cpu=8G
##SBATCH --tmp=32G

module load vital-it
module add UHTS/Analysis/samtools/1.8;
 
sourceBamDIR=$PWD/alignedBWA/picdup
 
#ls ${sourceBamDIR}/*.bam > bamList.txt
# then manually add a column with fraction to be sub sampled
 
destBamDIR=${sourceBamDIR%/picdup}
 
destBamDIR=${destBamDIR}/subBam
 
mkdir $destBamDIR
 
 
bamListFile=./bamList.txt
 
bamFiles=(`cut -f1 ${bamListFile}`)
fractions=(`cut -f2 ${bamListFile}`)
 
#numFiles=25

## to run as a loop
#for i in $(seq 0 ${numFiles})
#do
#	echo $i
#	if [[ "${fractions[$i]}" != "1" ]]
#	outFile=`basename ${bamFiles[$i]}`
#	then
#		samtools view -s ${fractions[$i]} -b -o  ${destBamDIR}/${outFile} ${bamFiles[$i]}
#		echo "subsampling " $outFile
#	else
#		cp ${bamFiles[$i]}  ${destBamDIR}/${outFile}
#		echo "copying " $outFile
#	fi
#done


let i=${SLURM_ARRAY_TASK_ID}-1

echo $i
outFile=`basename ${bamFiles[$i]}`
if [[ "${fractions[$i]}" != "1" ]]
then
	samtools view -s ${fractions[$i]} -b -o  ${destBamDIR}/${outFile} ${bamFiles[$i]}
	echo "subsampling " $outFile
else
	cp ${bamFiles[$i]}  ${destBamDIR}/${outFile}
	echo "copying " $outFile
fi
 
