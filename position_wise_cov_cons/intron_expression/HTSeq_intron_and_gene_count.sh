#!/bin/bash
#SBATCH -c 16
#SBATCH -p normal
#SBATCH --mem=64g
#SBATCH --time=12:00:00


sample=$1

genomeGTF=/path/hg38v33_ctat_introns_noOverlap.gtf


htseq-count -m intersection-nonempty -i gene_id -r pos -s yes ./${sample}/${sample}.star.bam $genomeGTF > ./${sample}/${sample}.htseq.introns.tsv

