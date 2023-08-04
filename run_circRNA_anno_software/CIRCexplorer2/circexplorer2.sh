#!/bin/bash
#SBATCH -c 16
#SBATCH -p normal
#SBATCH --mem=64g
#SBATCH --time=06:00:00

# sbatch circexplorer.sh dataset_name /path/to/SRR_1 /path/to/SRR_2

date1=$(date +"%s")

echo "========= Job started  at `date` =========="
echo "Command line options:"
echo $@


sample=$1
read_1=$2
read_2=$3
genomeDir=/path/ctat_hg38/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx
genome=${genomeDir}/../ref_genome.fa
refflat=/path/hg38_kg.txt
scratch=/scratch/$SLURM_JOBID


destination=$(pwd -P)
mkdir -p ${scratch}/${sample}
cd ${scratch}/${sample}
cp -p -r $0 .

# STAR
STAR --chimSegmentMin 10 --runThreadN 16 --genomeDir ${genomeDir} --readFilesIn ${read_1} ${read_2} --readFilesCommand zcat
cat Aligned.out.sam | samtools sort -T ./ --threads 4 -m 2G -o Aligned.out.bam - && rm -f Aligned.out.sam

# CE2
CIRCexplorer2 parse -t STAR Chimeric.out.junction > CIRCexplorer2_parse.log
CIRCexplorer2 annotate -r ${refflat} -g ${genome} -b back_spliced_junction.bed -o circularRNA_known.txt > CIRCexplorer2_annotate.log

pigz Chimeric.out.junction

# Copy to folder
cp -a ${scratch}/${sample} ${destination}/${sample}_tmp
mv ${destination}/${sample}_tmp ${destination}/${sample}

exit 0


