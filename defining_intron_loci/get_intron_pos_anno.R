#############################################
##### Setting up intron annotation data #####
#############################################

# NOTE: Due to the time needed to extract every intron using this script, it was run multiple times, by splitting and looping over the gene annotations in 3-4 parts.

# In terminal: conda activate r_env

#########################
##### Load packages #####
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
#library(data.table)

#####################
##### Load data #####

# hg38 v33 uses Ensembl 99
# For CIRCexplorer2 anno origin see: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/ 
# GENCODE releases: https://www.gencodegenes.org/human/releases.html 

##  CTAT Plug-n-play hg38 v33 annoations ##
gtf33_v33 <- rtracklayer::import("/path/ctat_hg38/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf")
gtf33_v33_df <- as.data.frame(gtf33_v33)


###########################################
##### Extract intron for each isoform #####

intron_anno <- gtf33_v33_df[, -c(6,8,9,13:15,19:26)]
intron_anno <- unique.data.frame(intron_anno)
dim(intron_anno)

intron_anno <- intron_anno %>% filter(type == "exon")
intron_anno2 <- as.data.frame(NULL)

# Replace NA (so logical argument can be used) with "missing"
intron_anno[is.na(intron_anno)] = "missing"
dim(intron_anno) 

transcripts <- unique(intron_anno$transcript_id) # length: 225,618, half: 112809
length(transcripts)

for (id in transcripts[56405:(112810-1)]) { # CHANGE HERE TO ONLY LOOP OVER PART OF THE GENE ANNOTATIONS. v2: 112810:length(transcripts), v1_1: 1:(56405-1), v1_2: 56405:(112810-1)
  exon_subset <- intron_anno[which(intron_anno$transcript_id == id),]
  if (nrow(exon_subset) > 1 & exon_subset$strand[1] == "+"){
    for (idx in 1:(nrow(exon_subset)-1)) {
      
      intron_anno2 <- rbind(intron_anno2, t(c(exon_subset$seqnames[idx],
                                              exon_subset$end[idx],
                                              exon_subset$start[idx+1],
                                              exon_subset$strand[idx],
                                              exon_subset$transcript_id[idx],
                                              exon_subset$transcript_name[idx],
                                              exon_subset$transcript_type[idx],
                                              exon_subset$gene_id[idx],
                                              exon_subset$gene_name[idx],
                                              exon_subset$gene_type[idx])))
    }
  } 
    
  if (nrow(exon_subset) > 1 & exon_subset$strand[1] == "-"){
    for (idx in 1:(nrow(exon_subset)-1)) {
      
      intron_anno2 <- rbind(intron_anno2, t(c(exon_subset$seqnames[idx],
                                              exon_subset$start[idx],
                                              exon_subset$end[idx+1],
                                              exon_subset$strand[idx],
                                              exon_subset$transcript_id[idx],
                                              exon_subset$transcript_name[idx],
                                              exon_subset$transcript_type[idx],
                                              exon_subset$gene_id[idx],
                                              exon_subset$gene_name[idx],
                                              exon_subset$gene_type[idx])))
    }
  }
    
  pos <- which(transcripts == id)
  len <- length(transcripts)
  print((pos-56405)/(len-(112810-56405)))
}

save(intron_anno2, file = "/path/CTAT_intron_pos_anno_ALL_1_2.Rdata")


######################################################################################################################################################
################################################## RUN AFTER SBATCH JOB WITH ABOVE CODE! #############################################################
######################################################################################################################################################

# Comment out if running above script on cluster using sbatch


#########################
##### Setting up DF #####

load("/path/CTAT_intron_pos_anno_ALL_1_1.Rdata")
intron_anno1_1 <- intron_anno2
load("/path/CTAT_intron_pos_anno_ALL_1_2.Rdata")
intron_anno1_2 <- intron_anno2
load("/path/CTAT_intron_pos_anno_ALL_2.Rdata")
intron_anno <- rbind(intron_anno1_1, intron_anno1_2, intron_anno2)

colnames(intron_anno) <- c("chr", "start", "end", "strand", "transcript_id", "transcript_name", "transcript_type", "gene_id", "gene_name", "gene_type")

cols <- c("start", "end", "chr", "strand")
intron_anno[, cols] <- sapply(intron_anno[, cols], as.character)
cols <- c("start", "end")
intron_anno[, cols] <- sapply(intron_anno[, cols], as.numeric)

intron_anno$chr[which(intron_anno$chr == "24")] <- "Y"
intron_anno$chr[which(intron_anno$chr == "23")] <- "X"
intron_anno$strand[which(intron_anno$strand == "1")] <- "+"
intron_anno$strand[which(intron_anno$strand == "2")] <- "-"

intron_anno$chr <- paste("chr", sep = "", intron_anno$chr)
intron_anno$start <- intron_anno$start+1
intron_anno$end <- intron_anno$end-1

intron_anno[which((intron_anno$end - intron_anno$start) < 0), c(2,3)] <- intron_anno[which((intron_anno$end - intron_anno$start) < 0), c(3,2)]
intron_anno <- intron_anno[which((intron_anno$end - intron_anno$start) != 0),]

intron_anno <- unnest(intron_anno, transcript_id)

# Export
intron_ranges <- makeGRangesFromDataFrame(intron_anno, 
                                                seqnames.field = "chr",
                                                start.field = "start",
                                                end.field = "end",
                                                strand.field = "strand",
                                                keep.extra.columns = T)

rtracklayer::export(intron_ranges, con = "/path/CTAT_intron_pos_ALL.gtf", format = "gtf")


