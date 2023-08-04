#!/bin/bash
#SBATCH -c 16
#SBATCH -p normal
#SBATCH --mem=64g
#SBATCH --time=2:00:00


# Run STAR and DCC per sample
# Possible to first run STAR and then DCC.
sample=$1
read_1=$2
read_2=$3
genomeDir=/path/star_fusion/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx
genome=${genomeDir}/../ref_genome.fa
refflat=/path/DCC/gencode.v33.annotation.gtf
repeats=/path/DCC/hg38_repeats.gtf

cd /path/out/221222_DCC

## STAR ##
mkdir -p ./${sample}_pairs ./${sample}_mate1 ./${sample}_mate2

STAR --runThreadN 16 --genomeDir ${genomeDir} --outSAMtype BAM Unsorted --readFilesIn ${read_1} ${read_2} --readFilesCommand zcat --outFileNamePrefix ./${sample}_pairs/${sample}_pairs_ --outReadsUnmapped Fastx --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2  --chimSegmentMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15 

echo "${sample} pairs done" 

STAR --runThreadN 16 --genomeDir ${genomeDir}  --outSAMtype None --readFilesIn ${read_1} --readFilesCommand zcat --outFileNamePrefix ./${sample}_mate1/${sample}_mate1_ --outReadsUnmapped Fastx --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2 --chimSegmentMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15 

echo "${sample} mate1 done" 

STAR --runThreadN 16 --genomeDir ${genomeDir} --outSAMtype None --readFilesIn ${read_2}  --readFilesCommand zcat --outFileNamePrefix ./${sample}_mate2/${sample}_mate2_ --outReadsUnmapped Fastx --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2 --chimSegmentMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15 

echo "${sample} mate2 done"

echo "/path/out/221222_DCC/${sample}_pairs/${sample}_pairs_Chimeric.out.junction" > ./input_DCC/${sample}.samplesheet
echo "/path/out/221222_DCC/${sample}_mate1/${sample}_mate1_Chimeric.out.junction" > ./input_DCC/${sample}.mate1
echo "/path/out/221222_DCC/${sample}_mate2/${sample}_mate2_Chimeric.out.junction" > ./input_DCC/${sample}.mate2


## DCC ##
mkdir ./${sample}_DCC_out
cd ./${sample}_DCC_out

DCC @/path/221222_DCC/input_DCC/${sample}.samplesheet -mt1 @/path/221222_DCC/input_DCC/${sample}.mate1 -mt2 @/path/221222_DCC/input_DCC/${sample}.mate2 -D -N -R ${repeats} -F -M -Pi -an ${refflat} -Nr 1 1 -fg 

echo "DCC ${sample} done"
