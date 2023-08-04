#############################################
##### Setting up intron annotation data #####
#############################################

# Intron anno are updated so circular sisRNAs that can be produced from smaller introns are assigned that intron and its corresponding annotations (e.g. gene and transcript type)
# Multiple gene annotation files are used due to difference in how comprehensive the annotations are as well as ensuring consistency with annotations used by CIRCexplorer2.


#########################
##### Load packages #####
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)


#####################
##### Load data #####

# User mount path:
mount_path <- "/path"

##  CTAT Plug-n-play hg38 v33 annotations ##
gtf33_v33 <- rtracklayer::import(paste(mount_path, "/path/star_fusion/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf", sep = ""))
gtf33_v33_df <- as.data.frame(gtf33_v33)

## NMIBC ##
load(paste(mount_path, "/path/200805_circexplorer2/CIRCexplorer2_NMIBC.Rdata", sep = ""))

## Fractionated BC cell lines ##
## HepG2 ##
## K562 ##
## ENCODE ##


###########################################
##### Rename DF so code can be reused #####

## NMIBC ##
ciRNA_anno <- ciRNAs_circAnno_NMIBC 
ciRNA_exp <- ciRNAs_circTable.mat_NMIBC
ciRNA_exp_norm <- ciRNAs_circTable.mat_NMIBC_norm


#####################################
##### Add gene type annotations #####

## Match transcript ID and gene name ##
gtf_anno <- gtf33_v33_df[, c(11:12, 16:18)]
gtf_anno <- unique.data.frame(gtf_anno)
colnames(gtf_anno)[2:3] <- c("geneName", "isoformName")
df_anno1 <- left_join(ciRNA_anno, gtf_anno[,3:5], by = c("isoformName"))
df_anno2 <- left_join(df_anno1, gtf_anno[,1:2], by = c("geneName"))
df_anno2 <- unique.data.frame(df_anno2)
dim(df_anno2)

## Collapse multiple annotations into one, so one row for each ciRNA ##
df_anno2 <- df_anno2 %>% group_by(name) %>% # Group by ciRNA
  mutate(
    gene_type = paste(gene_type, collapse = ", "),
    transcript_type = paste(transcript_type, collapse = ", ")
  ) %>%
  ungroup() 
df_anno2 <- unique.data.frame(df_anno2)
dim(df_anno2) 

ciRNA_anno <- df_anno2


##################################################
##### Find smallest ciRNA containing introns #####

# Some flanking introns called by circexplorer2 are clearly not the likely ciRNA producing introns (sometimes more than 5x longer than ciRNA, with other intron more fitting given what we know about 3' SS to BP distance)
# Update using same annotations as used in CIRCexplorer2 pipeline (might just randomly annotate to one isoform if multiple can produce the ciRNA)

gtf33 <- rtracklayer::import(paste(mount_path, "/path/hg38_kg.gtf", sep = "")) # GTF file containing annotations used as input in circexplorer2
gtf33_df <- as.data.frame(gtf33)
gtf_ctat_v33 <- file.path(paste(mount_path, "/path/hg38_kg.gtf", sep = ""))
txdb <- makeTxDbFromGFF(file = gtf_ctat_v33, format = "gtf")

## Add intron position columns (represent circexplorer2 called introns) ##
ciRNA_anno <- ciRNA_anno %>% separate(flankIntron, c("intron_chr", "intron_pos"), sep = ":", remove = F)
ciRNA_anno <- ciRNA_anno %>% separate(intron_pos, c("intron_start", "intron_end"), sep = "-")
cols <- c("start", "end", "intron_start", "intron_end", "sumReads", "max", "sumExp", "median", "sumReads_lin")
ciRNA_anno[, cols] <- sapply(ciRNA_anno[, cols], as.numeric)

colnames(ciRNA_anno)[14] <- "width" 
ciRNA_anno$intron_width <- ciRNA_anno$intron_end - ciRNA_anno$intron_start +1 

## Find smallest overlapping introns ##

# Find smallest GENCODE intron that overlaps ciRNA
intron_ciRNA_ranges <- makeGRangesFromDataFrame(ciRNA_anno[, -c(6, 9, 11:24)], 
                                                seqnames.field = "chr",
                                                start.field = "start",
                                                end.field = "end",
                                                strand.field = "strand",
                                                keep.extra.columns = T)

intronic_parts <- intronsByTranscript(txdb) # ALL introns from all isoforms
intronic_parts <- unlist(intronic_parts)

hits <- findOverlaps(intron_ciRNA_ranges, intronic_parts, type = "within")

intron_ciRNA_ranges_df <- as.data.frame(intron_ciRNA_ranges)
intronic_parts_df <- as.data.frame(intronic_parts, row.names = NULL)

intron_ciRNA_ranges_df <- intron_ciRNA_ranges_df[queryHits(hits),] # few ciRNAs not included
intronic_parts_df <- intronic_parts_df[subjectHits(hits),] # match ciRNA and intron anno rows

intron_ranges <- cbind(intron_ciRNA_ranges_df, intronic_parts_df) # bind the two df together, as the rows match.

colnames(intron_ranges)[c(1, 10:14)] <- c("chr",  "intron_chr", "intron_start", "intron_end", "intron_width", "intron_strand")

intron_ranges$keep <- ifelse((intron_ranges$strand == "-" & intron_ranges$end == intron_ranges$intron_end), "yes", 
                             ifelse((intron_ranges$strand == "+" & intron_ranges$start == intron_ranges$intron_start), "yes", "no")) 
intron_ranges <- intron_ranges[which(intron_ranges$keep == "yes"),] # Remove all intron anno that does not share a 5' SS with the overlapped ciRNA.
dim(intron_ranges)
length(unique(intron_ranges$name)) # NMIBC: 11326, 16 ciRNAs removed by this step (doesn't have intron anno in ctat gtf file with same 5' SS)

intron_ranges <- intron_ranges %>% group_by(name, strand) %>% # Group by ciRNA location and strand
  summarise(intron_start = intron_start[which(intron_width == min(intron_width))], # Only include smallest intron anno
            intron_end = intron_end[which(intron_width == min(intron_width))],
            intron_width = intron_width[which(intron_width == min(intron_width))]
  ) %>%
  ungroup() 
intron_ranges <- unique.data.frame(intron_ranges)
dim(intron_ranges) # NMIBC: 11326

intron_ranges <- intron_ranges[match(ciRNA_anno$name, intron_ranges$name),] # Match ciRNAs in old and new intron anno df
dim(intron_ranges)

ciRNA_anno <- cbind(ciRNA_anno, intron_ranges)
colnames(ciRNA_anno)[27:29] <- c("gtf_intron_start", "gtf_intron_end", "gtf_intron_width") 
ciRNA_anno <- ciRNA_anno[,-c(25:26)]

ciRNA_anno_sub <- ciRNA_anno %>% filter(gtf_intron_width < intron_width) %>% # Update intron anno where ciRNA can be found in smaller intron
  mutate(intron_width = gtf_intron_width,
         intron_start =  gtf_intron_start,
         intron_end = gtf_intron_end,
         flankIntron = paste0(chr, sep = ":", gtf_intron_start, sep = "-", gtf_intron_end)) 
dim(ciRNA_anno_sub) 

## Update intron annotations to match newly assigned introns ##

# Find matches from extracted introns from ref_annot.gtf in plug-n-play CTAT. Here, we find exact matches.
# See ./get_intron_pos_anno.R for constructing the file CTAT_intron_pos_ALL.gtf with all intron positions and annotations from the hg38 v33 CTAT gtf file.
gtf_intron_anno <- rtracklayer::import(paste(mount_path, "/path/CTAT_intron_pos_ALL.gtf", sep = ""))
gtf_intron_anno_df <- as.data.frame(gtf_intron_anno)
gtf_intron_anno_df <- gtf_intron_anno_df[, -c(6:9)]
gtf_intron_anno_df$name <- paste0(gtf_intron_anno_df$seqnames, sep = ":", gtf_intron_anno_df$start, sep = "-", gtf_intron_anno_df$end) # NMIBC: 431 equal to updated flankIntron

gtf_intron_anno_df_sub <- gtf_intron_anno_df[which(gtf_intron_anno_df$gene_name %in% ciRNA_anno_sub$geneName),] # check for exact matches 
length(which(ciRNA_anno_sub$flankIntron %in% gtf_intron_anno_df$name)) # NMIBC: 431, not all introns have exact matches in the gtf manually extracted intron positions

# Try to match gtf intron df with updated ciRNA-introns
intron_ranges2 <- makeGRangesFromDataFrame(gtf_intron_anno_df, 
                                           seqnames.field = "seqnames",
                                           start.field = "start",
                                           end.field = "end",
                                           strand.field = "strand",
                                           keep.extra.columns = T)

ciRNA_ranges_sub <- makeGRangesFromDataFrame(ciRNA_anno_sub[,-c(1:3, 6:7, 9, 14:27)], 
                                             seqnames.field = "intron_chr",
                                             start.field = "intron_start",
                                             end.field = "intron_end",
                                             strand.field = "strand",
                                             keep.extra.columns = T)

hits <- findOverlaps(ciRNA_ranges_sub, intron_ranges2, type = "any") # Find all overlaps, trim later.
intron_ranges_df2 <- as.data.frame(intron_ranges2)

intron_ranges_df2 <- intron_ranges_df2[subjectHits(hits),] # Match rows
ciRNA_anno_sub2 <- ciRNA_anno_sub[queryHits(hits),] # NMIBC: 8788

dist1 <- abs(ciRNA_anno_sub2$intron_start - intron_ranges_df2$start) # difference in 3' and 5' ends of overlapping flankIntrons and gtf introns
dist2 <- abs(ciRNA_anno_sub2$intron_end - intron_ranges_df2$end)

keep <- ifelse(dist1 <= 5 & dist2 <= 5, T, F) # restraint on how much matching introns can differ at 3' and 5' ends.
length(which(keep == T)) # NMIBC: T = 3052 (ones matching close enough to be assumed the same, thus having same annotations)

ciRNA_anno_sub2 <- ciRNA_anno_sub2[keep,] # NMIBC: 3052
intron_ranges_df2 <- intron_ranges_df2[keep,]
length(unique(ciRNA_anno_sub2$name)) # 927

# No need to update intron annotation, when all overlapping introns (that only differs at a few bases in intron start and end) are from the same gene
# But we still need to update intron region.
dont_update <- ciRNA_anno_sub2$name[which(ciRNA_anno_sub2$flankIntron == intron_ranges_df2$name & ciRNA_anno_sub2$geneName == intron_ranges_df2$gene_name)] # NMIBC: 413 unique flankIntrons
dont_update <- which(ciRNA_anno_sub2$name %in% dont_update)
ciRNA_update1 <- ciRNA_anno_sub2[dont_update,] # Only update flankIntron region, but keep annotation, as updated intron are from the same gene as original intron region.
ciRNA_anno_sub2 <- ciRNA_anno_sub2[-dont_update,]
intron_ranges_df2 <- intron_ranges_df2[-dont_update,]

update <- c(rep(NA, length(ciRNA_anno_sub2$name))) # if ciRNA flankIntron have multiple matches, where one originates from same gene as old flankIntron, we keep the original annotation.
for (name in unique(ciRNA_anno_sub2$name)){
  genes <- intron_ranges_df2$gene_name[which(ciRNA_anno_sub2$name == name)] # extract all potential gene annotations for ciRNA
  update[which(ciRNA_anno_sub2$name == name)] <- ifelse(length(unique(genes)) == 1 & ciRNA_anno_sub2$geneName[which(ciRNA_anno_sub2$name == name)] == intron_ranges_df2$gene_name[which(ciRNA_anno_sub2$name == name)], "no", "yes")
}

update <- which(update == "yes")
ciRNA_update1 <- rbind(ciRNA_update1, ciRNA_anno_sub2[-update,]) # flankIntron where we keep old anno for new intron region. Note, that removing transcript anno leaves a DF with one row for each ciRNA.
ciRNA_anno_sub2 <- ciRNA_anno_sub2[update,] # flankIntron where annotations needs to be updated.
intron_ranges_df2 <- intron_ranges_df2[update,]
intron_ranges_df2$name_id <- paste(intron_ranges_df2$name, sep = "_", intron_ranges_df2$transcript_id)

# Update remaining intron annotation if closest intron matching new intron positions are from a different gene.
keep_intron_anno <- c(rep(NA, length(intron_ranges_df2$name)))
for (name in unique(ciRNA_anno_sub2$name)){
  idx <- which(ciRNA_anno_sub2$name == name) # Index rows with same ciRNA
  diff_width <- abs(intron_ranges_df2$width[idx] - unique(ciRNA_anno_sub2$gtf_intron_width[idx])) # How closely do the remaining gtf introns match new flankIntron width for each ciRNA.
  closest_match <- min(diff_width) # Which matches closely (note that diff gtf Introns likely represent different isoforms or non-coding species)
  idx_min <- which(diff_width == closest_match) # Find row(s) with anno
  if (length(idx_min) > 1){
    # Rank overlapping closest match gtf introns by (which anno we want to use): 1) From same gene, 2) from protein-coding gene, 3) From none of previous mentioned categories.
    intron_keep <- ifelse(intron_ranges_df2$gene_name[idx][idx_min] == ciRNA_anno_sub2$geneName[idx][idx_min], "same", ifelse(intron_ranges_df2$gene_type[idx][idx_min] == "protein_coding", "PC", "none"))
    keep_intron_anno[idx] <- ifelse("same" %in% intron_keep, intron_ranges_df2$name_id[idx][idx_min][which(intron_keep == "same")[1]], ifelse("PC" %in% intron_keep, intron_ranges_df2$name_id[idx][idx_min][which(intron_keep == "PC")[1]], intron_ranges_df2$name_id[idx][idx_min][which(intron_keep == "none")[1]]))     
  } 
  if (length(which(diff_width == closest_match)) == 1){
    keep_intron_anno[idx] <- intron_ranges_df2$name_id[idx][idx_min] # Use closest gtf intron match for annotations 
  }
}

ciRNA_anno_sub2 <- ciRNA_anno_sub2[keep_intron_anno == intron_ranges_df2$name_id,] # NMIBC: 46, Only have one row per ciRNA (one set of annotations)
intron_ranges_df2 <- intron_ranges_df2[keep_intron_anno == intron_ranges_df2$name_id,]
ciRNA_anno_sub2[duplicated(ciRNA_anno_sub2$name),] # NMIBC: none

## make final DFs with positions and annotations that need to be updated ##
# Only update flankIntron 
ciRNA_update1 <- unique.data.frame(ciRNA_update1[, -c(8, 21,22)])
dim(ciRNA_update1) # NMIBC: 881

# Update flankIntron and gene annotation
ciRNA_anno_sub2[, c(7, 23)] <- intron_ranges_df2[, c(10, 11)]

## Update intron annotation for relevant ciRNAs ##
ciRNA_anno[match(ciRNA_update1$name, ciRNA_anno$name),10] <- ciRNA_update1[,9]
ciRNA_anno[match(ciRNA_anno_sub2$name, ciRNA_anno$name),] <- ciRNA_anno_sub2

ciRNA_anno <- ciRNA_anno[, -c(6, 8, 9, 11:13, 21:22, 24:27)]


###########################################
##### Save new annotation data frames #####

path = "/path/200805_circexplorer2/"

# NMIBC
CIRCexplorer2_NMIBC = c("CIRCexplorer2_NMIBC", "ciRNA_anno", "ciRNA_exp", "sampleAnno_NMIBC", "lin_exp", "ciRNA_exp_norm", "lin_exp_norm")
save(list=CIRCexplorer2_NMIBC, file=paste0(path, "Updated_CIRCexplorer2_NMIBC.Rdata"))


