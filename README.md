# bTMP-ChIP-seq-analysis
Pipeline for preprocessing and analysing ChIP-seq data. 
Adjust all scripts according to preference and change the paths to files accordingly.
Run trimming.sh script to remove adaptor sequences and clip where necessary (fastQC included before and after)
Run torunchipalignments.sh which runs chipalignments.sh to align data to reference genome using bowtie2
Run torunremdup.sh which runs removedup.sh to remove dupilcate reads using picard
Run normalization.sh which runs Normalization+8.R to normalize the data.
