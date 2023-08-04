#!/bin/bash
#SBATCH -c 16
#SBATCH -p normal
#SBATCH --mem=64g
#SBATCH --time=04:00:00

# run per sample
sample=$1
read_1=$2
read_2=$3

genomeDir=/path/hg38/software_indexes/star_fusion/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx
genome=${genomeDir}/../ref_genome.fa

cd /output/path/folder
mkdir -p ./${sample}
cd ./${sample}

# run STAR
/path/miniconda3/envs/circexplorer2_2/bin/runStar.pl --inFile1 ${read_1} --inFile2 ${read_2} --genomeDir ${genomeDir} --maxMismatch 0.02 --outPrefix ./${sample}

# run circRNA_finder
/path/miniconda3/envs/circexplorer2_2/bin/postProcessStarAlignment.pl --starDir ./ --minLen 50 --outDir ./