#!/bin/bash
#SBATCH -c 16
#SBATCH -p normal
#SBATCH --mem=64g
#SBATCH --time=12:00:00

# Generate anno & exp matrices
Rscript /path/generate_exp_matrices.R