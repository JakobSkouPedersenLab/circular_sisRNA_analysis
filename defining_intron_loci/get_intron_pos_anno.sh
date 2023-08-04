#!/bin/bash
#SBATCH --mem 32G
#SBATCH -c 4
#SBATCH -t 42:00:00

Rscript /path/defining_intron_loci/get_intron_pos_anno.R
