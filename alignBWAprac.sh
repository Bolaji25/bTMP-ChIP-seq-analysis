#! /usr/bin/bash

## Allocate resources
#SBATCH --time=1-40:00:00
#SBATCH --array=1

#SBATCH --mail-user=bolaji.isiaka@izb.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="alignentbwaprac"
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





# Use precompiled binaries (recommended)
curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.0pre2/bwa-mem2-2.0pre2_x64-linux.tar.bz2 \
  | tar jxf -
mkdir ./alnBWApracpipe
mkdir ./alnBWApracpipe/bam

bwa-mem2-2.0pre2_x64-linux/bwa-mem2 index /home/ubelix/izb/bi18k694/genomeversion/ws270/c_elegans.PRJNA13758.WS270.genomic.fa
bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem /home/ubelix/izb/bi18k694/genomeversion/ws270/c_elegans.PRJNA13758.WS270.genomic.fa  ./trimmed/*.fq | samtools view -b -o ./alnBWApracpipe/bam/
 



mkdir ./alnBWApracpipe/fastQC
mkdir ./alnBWApracpipe/fastQC/1
mkdir ./alnBWApracpipe/sorted
samtools sort -@ ./alnBWApracpipe/bam/*.bam -o  | samtools flagstat > ./alnBWApracpipe/fastQC/1/ 

# keep only reads that have Q>=30, both are mapped and in a FR or RF orientation.
#mkdir ./alnBWAprac/filtered
	bamtools filter -in ./alnBWAprac/sorted/${bname}.sorted.bam -out ./alnBWAprac/filtered/${bname}_filtered.bam

# get alignment stats again post-filtering
#mkdir ./alnBWAprac/fastQC/2
	samtools flagstat ./alnBWAprac/filtered/${bname}.filtered.bam  > ./alnBWAprac/fastQC/2/${bname}.filtered.txt

#sort by query name for duplicate removal
#mkdir ./alnBWAprac/qsort
	curl -fsSL https://github.com/broadinstitute/picard/releases/download/2.21.8/picard.jar -o picard.jar

	java -Xms1g -Xmx8g -jar ./picard.jar SortSam I=./alnBWAprac/filtered/${bname}.filtered.bam O=./alnBWAprac/qsort/${bname}.qsort.bam SORT_ORDER=queryname TMP_DIR=${TMPDIR}


# remove duplicates with picard (path to picard should be set in $PICARD variable in .bash_profile or in session)
# Note that to mark unmapped mates of mapped records and supplementary/secondary alignments as duplicates the bam
# file must be querysorted (by name) not by coordinate. 
#mkdir ./alnBWAprac/noDup
	java -Xms1g -Xmx8g -jar ./picard.jar MarkDuplicates I=./alnBWAprac/qsort/${bname}.qsort.bam O=./alnBWAprac/noDup/${bname}.noDup.bam M=./alnBWAprac/fastQC/${bname}_picard.txt REMOVE_DUPLICATES=false REMOVE_SEQUENCING_DUPLICATES=true  ASSUME_SORT_ORDER=queryname TMP_DIR=${TMP}

# get alignment stats again post-duplicate removal
#mkdir ./alnBWAprac/3
	samtools flagstat ./alnBWAprac/noDup/${bname}.noDup.bam > ./alnBWAprac/fastQC/3/${bname}.noDUP.txt

done
