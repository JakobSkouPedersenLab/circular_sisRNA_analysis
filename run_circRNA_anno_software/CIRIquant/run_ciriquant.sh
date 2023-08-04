#!/bin/bash
#SBATCH -c 16
#SBATCH -p normal
#SBATCH --mem=64g
#SBATCH --time=6:00:00

# run per sample
sample=$1
read_1=$2
read_2=$3
id=$4

cd /path/hg38/230113_CIRIquant_circexplorer2
mkdir ./${sample}

CIRIquant -t 64 -1 ${read_1} -2 ${read_2} --config /path/circ_detection/circexplorer2_CIRIquant/hg38v33.yml -o ./${sample} -p ${sample}.${id} --circ /path/hg38/200805_circexplorer2/${sample}/circularRNA_known.txt --tool CIRCexplorer2

# ./circularRNA_known.txt is the output of CIRCexplorer2