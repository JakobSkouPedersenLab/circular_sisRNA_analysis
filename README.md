# circular_sisRNA_analysis
[![DOI](https://zenodo.org/badge/665570589.svg)](https://zenodo.org/badge/latestdoi/665570589)

*Overview of R code used for circular sisRNA analysis*

The code in this repository represents a subset used to generate and analyse the results presented in: 

Asta M. Rasmussen†, Trine Line H. Okholm†, Michael Knudsen, Søren Vang, Lars Dyrskjøt, Thomas B. Hansen, and Jakob S. Pedersen. Circular Stable Intronic RNAs possess distinct biological features and are deregulated in bladder cancer. NAR Cancer (2023). 

For further requests please contact the corresponding author.

## Brief overview

All analyses are performed using the human hg38 reference genome.
Note, that most code are implemented in R.

Circular stable intronic sequence RNA (sisRNA) and circular intronic long non-coding RNA (ciRNA) are used interchangeably throughout the scripts.

Please see the workflow files in each subfolder for walkthroughs. 

### Subfolders

**./defining_intron_loci**

Code for extracting and representing intronic loci.

**./run_circRNA_anno_software**

Code used for running existing tools for circular sisRNA (intronic circular RNA) detection and annotation from total RNA-sequencing data.

Tools used:
- CIRCexplorer2
- CIRIquant
- DCC
- circRNA_finder

**./position_wise_cov_cons**

Tutorial and code for extracting position-wise intronic read-depth (and sequence conservation (phyloP-score))
as well as generating an intronic read-count matrix (using HTSeq-count).
