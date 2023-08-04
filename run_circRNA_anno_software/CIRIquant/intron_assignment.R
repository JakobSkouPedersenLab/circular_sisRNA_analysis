##########################################################################################
################# match CIRIquant-CE2 ciRNAs to non-overlapping introns ##################
##########################################################################################
# 24/01/2023

# Need to properly annotate circRNAs to identify sisRNA subset.

##### Load data and packages #####
library(tidyverse)
library(rtracklayer)

##### Load general data #####

# User mount path:
mount_path <- "/Users/path/to/mount"
path <- paste(mount_path, "/path/to/OutputFolder/out_CIRIquant_circexplorer2", sep = "")

# circRNA_Finder
load(paste(mount_path, "/path/to/CE2_CIRIquant_NMIBC.Rdata", sep = ""))
hg38_introns <- rtracklayer::import(paste(mount_path, "/iRNA2020/data/hg38/hg38v33_ctat_introns_noOverlap_withGeneId.gtf", sep = ""))


#############################################################
############# collapse circRNA sharing 5' end ###############

# add basic annotations
circCounts <- apply(circCounts, 2, as.numeric) #as.numeric(as.matrix(circCounts))
circCounts_fw <- apply(circCounts_fw, 2, as.numeric) 
circCounts_norm <- apply(circCounts_norm, 2, as.numeric) 
circCounts_ratio <- apply(circCounts_ratio, 2, as.numeric) 
circCounts_norm_old = t(t(circCounts)/sample_anno$totalReads)*10^6

rownames(circCounts) <- circTable$name
rownames(circCounts_fw) <- circTable$name
rownames(circCounts_norm) <- circTable$name
rownames(circCounts_ratio) <- circTable$name
rownames(circCounts_norm_old) <- circTable$name

circTable$sumReads = apply(circCounts, 1, sum)
circTable$mean = apply(circCounts, 1, mean)
circTable$max = apply(circCounts, 1, max)
circTable$sumExp = apply(circCounts, 1, function(x) sum(x > 0))
circTable$median = apply(circCounts, 1, median)

# width
#circTable$width <- circTable$end - circTable$start +1

##### Try making intron width a bit broader to capture circles #####
hg38_introns_df <- as.data.frame(hg38_introns)
hg38_introns_df <- hg38_introns_df[, -c(6:10)]
hg38_introns_df$intron_name <- paste0(hg38_introns_df$seqnames, sep = ":", hg38_introns_df$start, sep = "-", hg38_introns_df$end)
hg38_introns_df$tmp <- paste0(hg38_introns_df$start, sep = "-", hg38_introns_df$end)
hg38_introns_df$start <- hg38_introns_df$start -10
hg38_introns_df$end <- hg38_introns_df$end +10
hg38_introns <-  makeGRangesFromDataFrame(hg38_introns_df, 
                                          seqnames.field = "seqnames",
                                          start.field = "start",
                                          end.field = "end",
                                          strand.field = "strand",
                                          keep.extra.columns = T)


################################################
############# assign intron anno ###############

# Make circRNA GRanges
circRNA_ranges <- makeGRangesFromDataFrame(circTable, 
                                         seqnames.field = "chr",
                                         start.field = "start",
                                         end.field = "end",
                                         strand.field = "strand",
                                         keep.extra.columns = T)

hits <- findOverlaps(circRNA_ranges, hg38_introns, type = "within")
#hits <- findOverlaps(circRNA_ranges, hg38_introns, type = "any") # anyIntronOverlap
intron_ranges_df <- as.data.frame(hg38_introns)

intron_ranges_df <- intron_ranges_df[subjectHits(hits),] # Match rows
circRNA_anno <- circTable
circRNA_anno <- circRNA_anno[queryHits(hits),] # NMIBC: 
dim(intron_ranges_df)
dim(circRNA_anno)

#____________
# Reestablish correct intron coordinates
intron_ranges_df$start <- intron_ranges_df$start + 10
intron_ranges_df$end <- intron_ranges_df$end - 10
#____________

dist1 <- abs(circRNA_anno$start - intron_ranges_df$start)
dist2 <- abs(circRNA_anno$end - intron_ranges_df$end)

pot_sisRNA <- ifelse(dist1 <= 10 | dist2 <= 10, T, F) # 5' SS should be within 10 bp for introns to match.
length(which(pot_sisRNA == T)) # mostly sisRNA from ncRNA genes
circRNA_anno$pot_sisRNA <- pot_sisRNA # sisRNAs where we can match introns within 10 bp of 5' splice site.

idx <- which(circRNA_anno$pot_sisRNA)
sisRNA_df <- circRNA_anno[idx,]
intron_ranges_df <- intron_ranges_df[idx,]

intron_ranges_df$intron_name <- paste0(intron_ranges_df$seqnames, sep = ":", intron_ranges_df$start, sep = "-", intron_ranges_df$end)
#____
intron_ranges_df <- intron_ranges_df[,-11]
#____

intron_ranges_df$name_sisRNA <- sisRNA_df$name
sisRNA_df$intron_name <- intron_ranges_df$intron_name
sisRNA_df$intron_id <- intron_ranges_df$intron_id
dim(intron_ranges_df)

length(unique(intron_ranges_df$name_sisRNA)) 
length(unique(sisRNA_df$intron_name)) 

dim(sisRNA_df[duplicated(sisRNA_df$name),]) # NMIBC: 0 ciRNAs overlap on both strands

sisRNA_df <- sisRNA_df[,-c(7:10,19)]

## Collapse ciRNA data ##

# Collapse based on intron cluster definition

#merged_ciRNA_anno$intron_id <- ciRNA_anno$intron_id[match(merged_ciRNA_anno$name, ciRNA_anno$name)]
length(unique(sisRNA_df$intron_name)) 
merged_sisRNA_df <- sisRNA_df %>% group_by(intron_name) %>% summarise(name = name[which(sumReads == max(sumReads))[1]],
                                                                      chr = chr[which(sumReads == max(sumReads))[1]],
                                                                      start = start[which(sumReads == max(sumReads))[1]],
                                                                      end = end[which(sumReads == max(sumReads))[1]],
                                                                      strand = strand[which(sumReads == max(sumReads))[1]],
                                                                      width = width[which(sumReads == max(sumReads))[1]],
                                                                      collapsed = paste(unique(name), collapse = ","),
                                                                      geneName = paste(unique(gene_name), collapse = "-"),
                                                                      gene_id = paste(unique(gene_id), collapse = "-"),
                                                                      gene_type = paste(unique(gene_type), collapse = "-"),
                                                                      sumReads = sum(sumReads),
                                                                      mean = NA,
                                                                      max = NA,
                                                                      sumExp = NA,
                                                                      median = NA,
                                                                      pot_sisRNA = (T %in% pot_sisRNA),
                                                                      intron_id = intron_id[which(sumReads == max(sumReads))[1]]) %>% ungroup()
merged_sisRNA_df <- unique.data.frame(merged_sisRNA_df)
dim(merged_sisRNA_df)
length(unique(merged_sisRNA_df$intron_name))


## exp matrices ##
sisRNA_introns <- unique(sisRNA_df$intron_name)
length(sisRNA_introns)

sisRNA_exp <- as.data.frame(circCounts[match(sisRNA_df$name, rownames(circCounts)),])
sisRNA_exp_norm <- as.data.frame(circCounts_norm[match(sisRNA_df$name, rownames(circCounts_norm)),])
sisRNA_exp_norm_old <- as.data.frame(circCounts_norm_old[match(sisRNA_df$name, rownames(circCounts_norm_old)),])
sisRNA_exp_fw <- as.data.frame(circCounts_fw[match(sisRNA_df$name, rownames(circCounts_fw)),])
sisRNA_exp_ratio <- as.data.frame(circCounts_ratio[match(sisRNA_df$name, rownames(circCounts_ratio)),])

circCounts_fw_norm_old = t(t(circCounts_fw)/sample_anno$totalReads)*10^6
sisRNA_exp_fw_norm_old <- as.data.frame(circCounts_fw_norm_old[match(sisRNA_df$name, rownames(circCounts_fw_norm_old)),])

#identical(rownames(sisRNA_exp), sisRNA_df$name) # T
#length(which(rownames(sisRNA_exp) != sisRNA_df$name)) # 0
#length(which(rownames(sisRNA_exp_norm) != sisRNA_df$name)) # 0


## merged ##
merged_sisRNA_exp <- lapply(1:length(sisRNA_df$name), function(x) colSums(sisRNA_exp[which(sisRNA_df$intron_name == sisRNA_df$intron_name[x]),], na.rm = T))
merged_sisRNA_exp <- do.call(rbind, merged_sisRNA_exp)
merged_sisRNA_exp2 <- as.data.frame(merged_sisRNA_exp)
merged_sisRNA_exp2$intron_name <- sisRNA_df$intron_name
merged_sisRNA_exp2 <- unique.data.frame(merged_sisRNA_exp2)
rownames(merged_sisRNA_exp2) <- merged_sisRNA_df$name[match(merged_sisRNA_exp2$intron_name, merged_sisRNA_df$intron_name)]
merged_sisRNA_exp2 <- merged_sisRNA_exp2[,-458]
merged_sisRNA_exp2 <- merged_sisRNA_exp2[match(merged_sisRNA_df$name, rownames(merged_sisRNA_exp2)),]
identical(merged_sisRNA_df$name, rownames(merged_sisRNA_exp2))

merged_sisRNA_exp_norm <- lapply(1:length(sisRNA_df$name), function(x) colSums(sisRNA_exp_norm[which(sisRNA_df$intron_name == sisRNA_df$intron_name[x]),], na.rm = T))
merged_sisRNA_exp_norm <- do.call(rbind, merged_sisRNA_exp_norm)
merged_sisRNA_exp_norm2 <- as.data.frame(merged_sisRNA_exp_norm)
merged_sisRNA_exp_norm2$intron_name <- sisRNA_df$intron_name
merged_sisRNA_exp_norm2 <- unique.data.frame(merged_sisRNA_exp_norm2)
rownames(merged_sisRNA_exp_norm2) <- merged_sisRNA_df$name[match(merged_sisRNA_exp_norm2$intron_name, merged_sisRNA_df$intron_name)]
merged_sisRNA_exp_norm2 <- merged_sisRNA_exp_norm2[,-458]
merged_sisRNA_exp_norm2 <- merged_sisRNA_exp_norm2[match(merged_sisRNA_df$name, rownames(merged_sisRNA_exp_norm2)),]
identical(merged_sisRNA_df$name, rownames(merged_sisRNA_exp_norm2))

merged_sisRNA_exp_norm_old <- lapply(1:length(sisRNA_df$name), function(x) colSums(sisRNA_exp_norm_old[which(sisRNA_df$intron_name == sisRNA_df$intron_name[x]),], na.rm = T))
merged_sisRNA_exp_norm_old <- do.call(rbind, merged_sisRNA_exp_norm_old)
merged_sisRNA_exp_norm_old2 <- as.data.frame(merged_sisRNA_exp_norm_old)
merged_sisRNA_exp_norm_old2$intron_name <- sisRNA_df$intron_name
merged_sisRNA_exp_norm_old2 <- unique.data.frame(merged_sisRNA_exp_norm_old2)
rownames(merged_sisRNA_exp_norm_old2) <- merged_sisRNA_df$name[match(merged_sisRNA_exp_norm_old2$intron_name, merged_sisRNA_df$intron_name)]
merged_sisRNA_exp_norm_old2 <- merged_sisRNA_exp_norm_old2[,-458]
merged_sisRNA_exp_norm_old2 <- merged_sisRNA_exp_norm_old2[match(merged_sisRNA_df$name, rownames(merged_sisRNA_exp_norm_old2)),]
identical(merged_sisRNA_df$name, rownames(merged_sisRNA_exp_norm_old2))

merged_sisRNA_exp_fw <- lapply(1:length(sisRNA_df$name), function(x) colSums(sisRNA_exp_fw[which(sisRNA_df$intron_name == sisRNA_df$intron_name[x]),], na.rm = T))
merged_sisRNA_exp_fw <- do.call(rbind, merged_sisRNA_exp_fw)
merged_sisRNA_exp_fw2 <- as.data.frame(merged_sisRNA_exp_fw)
merged_sisRNA_exp_fw2$intron_name <- sisRNA_df$intron_name
merged_sisRNA_exp_fw2 <- unique.data.frame(merged_sisRNA_exp_fw2)
rownames(merged_sisRNA_exp_fw2) <- merged_sisRNA_df$name[match(merged_sisRNA_exp_fw2$intron_name, merged_sisRNA_df$intron_name)]
merged_sisRNA_exp_fw2 <- merged_sisRNA_exp_fw2[,-458]
merged_sisRNA_exp_fw2 <- merged_sisRNA_exp_fw2[match(merged_sisRNA_df$name, rownames(merged_sisRNA_exp_fw2)),]
identical(merged_sisRNA_df$name, rownames(merged_sisRNA_exp_fw2))

merged_sisRNA_exp_fw_norm_old <- lapply(1:length(sisRNA_df$name), function(x) colSums(sisRNA_exp_fw_norm_old[which(sisRNA_df$intron_name == sisRNA_df$intron_name[x]),], na.rm = T))
merged_sisRNA_exp_fw_norm_old <- do.call(rbind, merged_sisRNA_exp_fw_norm_old)
merged_sisRNA_exp_fw_norm_old2 <- as.data.frame(merged_sisRNA_exp_fw_norm_old)
merged_sisRNA_exp_fw_norm_old2$intron_name <- sisRNA_df$intron_name
merged_sisRNA_exp_fw_norm_old2 <- unique.data.frame(merged_sisRNA_exp_fw_norm_old2)
rownames(merged_sisRNA_exp_fw_norm_old2) <- merged_sisRNA_df$name[match(merged_sisRNA_exp_fw_norm_old2$intron_name, merged_sisRNA_df$intron_name)]
merged_sisRNA_exp_fw_norm_old2 <- merged_sisRNA_exp_fw_norm_old2[,-458]
merged_sisRNA_exp_fw_norm_old2 <- merged_sisRNA_exp_fw_norm_old2[match(merged_sisRNA_df$name, rownames(merged_sisRNA_exp_fw_norm_old2)),]
identical(merged_sisRNA_df$name, rownames(merged_sisRNA_exp_fw_norm_old2))

merged_sisRNA_exp_ratio <- lapply(1:length(sisRNA_df$name), function(x) colSums(sisRNA_exp_ratio[which(sisRNA_df$intron_name == sisRNA_df$intron_name[x]),], na.rm = T))
merged_sisRNA_exp_ratio <- do.call(rbind, merged_sisRNA_exp_ratio)
merged_sisRNA_exp_ratio2 <- as.data.frame(merged_sisRNA_exp_ratio)
merged_sisRNA_exp_ratio2$intron_name <- sisRNA_df$intron_name
merged_sisRNA_exp_ratio2 <- unique.data.frame(merged_sisRNA_exp_ratio2)
rownames(merged_sisRNA_exp_ratio2) <- merged_sisRNA_df$name[match(merged_sisRNA_exp_ratio2$intron_name, merged_sisRNA_df$intron_name)]
merged_sisRNA_exp_ratio2 <- merged_sisRNA_exp_ratio2[,-458]
merged_sisRNA_exp_ratio2 <- merged_sisRNA_exp_ratio2[match(merged_sisRNA_df$name, rownames(merged_sisRNA_exp_ratio2)),]
identical(merged_sisRNA_df$name, rownames(merged_sisRNA_exp_ratio2))

dim(merged_sisRNA_exp2)
dim(merged_sisRNA_exp_norm2)
dim(merged_sisRNA_exp_norm_old2)
dim(merged_sisRNA_exp_fw2)
dim(merged_sisRNA_exp_fw_norm_old2)
dim(merged_sisRNA_exp_ratio2)


## anno ##
merged_sisRNA_df$sumReads2 = apply(merged_sisRNA_exp2, 1, sum) # check if match
merged_sisRNA_df$mean = apply(merged_sisRNA_exp2, 1, mean)
merged_sisRNA_df$max = apply(merged_sisRNA_exp2, 1, max)
merged_sisRNA_df$sumExp = apply(merged_sisRNA_exp2, 1, function(x) sum(x > 0))
merged_sisRNA_df$median = apply(merged_sisRNA_exp2, 1, median)

identical(merged_sisRNA_df$sumReads, merged_sisRNA_df$sumReads2) # F
length(which(merged_sisRNA_df$sumReads != merged_sisRNA_df$sumReads2)) # 0
#length(merged_sisRNA_df$name != rownames(merged_sisRNA_exp)) #0
length(which(merged_sisRNA_df$name != rownames(merged_sisRNA_exp2)))


#######################################
############# save data ###############

sisRNA_anno <- sisRNA_df
merged_sisRNA_anno <- merged_sisRNA_df
intron_anno <- intron_ranges_df
merged_sisRNA_exp <- merged_sisRNA_exp2
merged_sisRNA_exp_norm <- merged_sisRNA_exp_norm2
merged_sisRNA_exp_norm_old <- merged_sisRNA_exp_norm_old2
merged_sisRNA_exp_fw <- merged_sisRNA_exp_fw2
merged_sisRNA_exp_ratio <- merged_sisRNA_exp_ratio2
merged_sisRNA_exp_fw_norm_old <- merged_sisRNA_exp_fw_norm_old2



CE2_CIRIquant_NMIBC = c("sisRNA_anno", "sisRNA_exp", "sisRNA_exp_norm", "sisRNA_exp_norm_old", "sisRNA_exp_ratio", "sisRNA_exp_fw", "sisRNA_exp_fw_norm_old", 
              "merged_sisRNA_anno", "merged_sisRNA_exp", "merged_sisRNA_exp_norm", "merged_sisRNA_exp_norm_old", "merged_sisRNA_exp_ratio", 
              "merged_sisRNA_exp_fw", "merged_sisRNA_exp_fw_norm_old", "sample_anno", "intron_anno")

save(list=CE2_CIRIquant_NMIBC, file=paste0(path, "/Updated_CE2_CIRIquant_NMIBC.Rdata"))


