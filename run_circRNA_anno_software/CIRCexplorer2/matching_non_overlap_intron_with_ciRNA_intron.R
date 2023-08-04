########################################################################################################################
########################################################################################################################
##### Matching htseq intron count with ciRNA anno data #################################################################
########################################################################################################################
########################################################################################################################


#########################
##### Load packages #####
library(tidyverse)
library(GenomicRanges)


#####################
##### Load data #####

mount_path <- "/path" # Change to match user's mount path

## NMIBC ##
load(paste(mount_path, "/path/200805_circexplorer2/Updated_CIRCexplorer2_NMIBC.Rdata", sep = ""))
load(paste(mount_path, "/path/Introns_NMIBC.Rdata", sep = ""))


###########################################
##### Rename DF so code can be reused #####

## NMIBC ##
sample_anno <- sampleAnno_NMIBC
ciRNA_anno_old <- ciRNA_anno

intron_anno <- intronAnno_NMIBC
intron_exp <- intron_counts.mat_NMIBC
intron_exp_norm <- intron_counts.mat_NMIBC_norm


###########################
##### Add annotations #####

## Extract intron position and and width ##

# Non-collapsed data
ciRNA_anno <- ciRNA_anno %>% separate(flankIntron, c("intron_chr", "intron_pos"), sep = ":", remove = F)
ciRNA_anno <- ciRNA_anno %>% separate(intron_pos, c("intron_start", "intron_end"), sep = "-")
cols <- c("start", "end", "intron_start", "intron_end", "sumReads", "max", "sumExp", "median", "sumReads_lin")
ciRNA_anno[, cols] <- sapply(ciRNA_anno[, cols], as.numeric)


############################################################
##### Matching htseq intron count with ciRNA anno data #####
############################################################

# Reads counted using gtf file with intron ranges not overlapping exons

## Non-collapsed ciRNA data ##
intron_anno <- intron_anno[, -c(6:10)]

intron_ranges <- makeGRangesFromDataFrame(intron_anno, 
                                          seqnames.field = "seqnames",
                                          start.field = "start",
                                          end.field = "end",
                                          strand.field = "strand",
                                          keep.extra.columns = T)

ciRNA_ranges <- makeGRangesFromDataFrame(ciRNA_anno[,-c(6:18)], 
                                         seqnames.field = "chr",
                                         start.field = "start",
                                         end.field = "end",
                                         strand.field = "strand",
                                         keep.extra.columns = T)

hits <- findOverlaps(ciRNA_ranges, intron_ranges, type = "any")
intron_ranges_df <- as.data.frame(intron_ranges)

intron_ranges_df <- intron_ranges_df[subjectHits(hits),] # Match rows
ciRNA_anno <- ciRNA_anno[queryHits(hits),] # NMIBC: 9830

dist1 <- abs(ciRNA_anno$start - intron_ranges_df$start)
dist2 <- abs(ciRNA_anno$end - intron_ranges_df$end)

keep <- ifelse(dist1 <= 10 | dist2 <= 10, T, F) # 5' SS should be within 10 bp for introns to match.
length(which(keep == T)) # NMIBC: keep = T = 9157

ciRNA_anno <- ciRNA_anno[keep,] # keep ciRNAs where we can match introns within 10 bp of 5' splice site.
intron_ranges_df <- intron_ranges_df[keep,]

intron_ranges_df <- cbind(intron_ranges_df, ciRNA_anno$name) # Match intron to ciRNA
colnames(intron_ranges_df)[c(15)] <- c("name")
intron_ranges_df$intron_name <- paste0(intron_ranges_df$seqnames, sep = ":", intron_ranges_df$start, sep = "-", intron_ranges_df$end)
dim(intron_ranges_df)

ciRNA_anno <- cbind(ciRNA_anno, intron_ranges_df$intron_name) # Match intron to ciRNA
colnames(ciRNA_anno)[c(19)] <- c("intron_name")

length(unique(intron_ranges_df$name)) # NMIBC: 9080; 2340 ciRNAs doesn't have matching intron (given the way we have defined our introns)

dim(ciRNA_anno[duplicated(ciRNA_anno$name),]) # NMIBC: 77 ciRNAs overlap min 2 introns

test <- lapply(1:nrow(intron_ranges_df), function(x) max(intron_ranges_df$width[which(intron_ranges_df$name == intron_ranges_df$name[x])])) # When long intron overlaps >1 intron, we find the largest of the two introns and use as representative.
test <- unlist(test)
keep <- intron_ranges_df$width == test # Only keep largest intron, when ciRNA overlaps multiple introns (remember the introns are non-overlapping)

ciRNA_anno <- ciRNA_anno[keep,] # NMIBC: 9080
intron_ranges_df <- intron_ranges_df[keep,]
dim(ciRNA_anno[duplicated(ciRNA_anno$name),]) # NMIBC: 0 duplicates
ciRNA_anno$intron_id <- intron_ranges_df$intron_id

ciRNA_anno_old$intron_id <- ciRNA_anno$intron_id[match(ciRNA_anno_old$name, ciRNA_anno$name)] 
ciRNA_anno <- ciRNA_anno_old

# NOTE: ONLY INCLUDE CIRNA-INTRON MATCHES WHERE THEY OVERLAP AND EITHER START OR END ARE WITHIN 10 BP FROM EACHOTHER


###############################################
#################### Save #####################

path = "/path/200805_circexplorer2/"

# NMIBC
CIRCexplorer2_NMIBC = c("CIRCexplorer2_NMIBC", "ciRNA_anno", "ciRNA_exp", "sample_anno", "lin_exp", "ciRNA_exp_norm", "lin_exp_norm")
save(list=CIRCexplorer2_NMIBC, file=paste0(path, "Updated_CIRCexplorer2_NMIBC.Rdata"))



