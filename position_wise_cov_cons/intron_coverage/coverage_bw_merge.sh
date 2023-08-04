#!/bin/bash
#SBATCH -c 16
#SBATCH -p normal
#SBATCH --mem=64g
#SBATCH --time=24:00:00

## STEP 1
# extract bw paths from /path/track_data_intron_stranded.txt
wiggletools mean $(</path/bw_paths.txt) > /path/mean_intron_coverage.bedGraph

## STEP 2
sort -k1,1 -k2,2n /path/mean_intron_coverage.bedGraph > /path/mean_intron_coverage_sorted.bedGraph

# for /path/hg38.chrom.sizes see https://github.com/astamr/ciRNA_Uromol_trackHub/blob/master/hg38/hg38.chrom.sizes 
bedGraphToBigWig /path/mean_intron_coverage_sorted.bedGraph /path/hg38.chrom.sizes /path/mean_intron_coverage_sorted.bw



