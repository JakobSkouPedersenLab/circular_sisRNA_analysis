#!/bin/bash
#SBATCH -c 16
#SBATCH -p normal
#SBATCH --mem=64g
#SBATCH --time=12:00:00

#_____________________________________________________
# Use UCSC/Kent tools to find summary stats of bw file
#cd /path/

# normalize (if relevant)
#awk '{ $4 =  ($4*1000000/total_read_depth)/nr_samples} 1' /path/mean_intron_coverage_sorted.bedGraph > /path/mean_intron_coverage_sorted_norm_perSample.bedGraph

#bedGraphToBigWig /path/mean_intron_coverage_sorted_norm_perSample.bedGraph /path/hg38.chrom.sizes /path/mean_intron_coverage_sorted_norm_perSample.bw

# summary stats
#bigWigAverageOverBed mean_intron_coverage_sorted_norm_perSample.bw /path/hg38v33_ctat_introns_noOverlap.bed mean_intron_coverage_sorted_norm_perSample.tab

#_____________________________________________________
# Extract conservation/coverage matching region for each position
##conda activate bwtool_env

cd /path/

bwtool extract bed /path/hg38v33_ctat_introns_noOverlap.bed mean_intron_coverage_sorted.bw mean_intron_coverage_sorted.bed

# Note: bwtools only have .00 precision, so using CPM makes everything 0.00! Instead, use raw count data and then normalize read-depth afterwadrs (as 10^6/total_read_depth)/nr_samples)

#______________________________________________________