#!/bin/bash
#SBATCH -c 16
#SBATCH -p normal
#SBATCH --mem=64g
#SBATCH --time=12:00:00

sample=$1

# Remove linear splice junction reads ("overlap" intron)
samtools view -h -F 4 ./$sample/$sample.star.bam | awk '$6 !~ /N/ || $1 ~ /^@/' | samtools view -bS - > ./$sample/$sample.noSpliceReads.bam

# Create Intron Coverage and Gene Coverage .BAM files
bedtools intersect -a ./$sample/$sample.noSpliceReads.bam -b /path/hg38v33_ctat_introns_noOverlap.bed -s > ./$sample/$sample.intron.stranded.bam
bedtools intersect -a ./$sample/$sample.star.bam -b /path/hg38v33_ctat_genes_noOverlap.bed -s > ./$sample/$sample.gene.stranded.bam

# Index .bam files
samtools index -b ./$sample/$sample.intron.stranded.bam
samtools index -b ./$sample/$sample.gene.stranded.bam

# create bigwig coverage file
bamCoverage --binSize 1 -p "max" -b ./$sample/$sample.intron.stranded.bam --outFileFormat bigwig -o ./$sample/$sample.intron.stranded.bw
bamCoverage --binSize 1 -p "max" -b ./$sample/$sample.gene.stranded.bam --outFileFormat bigwig -o ./$sample/$sample.gene.stranded.bw

# Remove temporary files
rm ./$sample/$sample.noSpliceReads.bam
rm ./$sample/$sample.intron.stranded.bam
rm ./$sample/$sample.intron.stranded.bam.bai
rm ./$sample/$sample.gene.stranded.bam
rm ./$sample/$sample.gene.stranded.bam.bai

# Create track id file used in Rscript 
echo "$sample NMIBC /path/$sample/$sample.intron.stranded.bw" >> /path/track_data_intron_stranded.txt

echo "$sample NMIBC /path/$sample/$sample.gene.stranded.bw" >> /path/track_data_gene_stranded.txt
