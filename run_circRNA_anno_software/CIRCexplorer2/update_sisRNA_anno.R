#################################################################################
##### Redefine merged sisRNAs as clusters - needs to share same intron anno #####
#################################################################################

# Intron anno are updated so sisRNAs that can be produced from smaller introns are assigned that intron and its corresponding annotations (e.g. gene and transcript type)
# Also update so only merge sisRNAs derived from same intron; sisRNA cluster.


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

load(paste(mount_path, "/path/200805_circexplorer2/Updated_CIRCexplorer2_NMIBC.Rdata", sep = ""))


#####################################
##### Update merged sisRNA anno #####

# Based on shared flankIntron anno
ciRNA_anno$SS5 <- ifelse(ciRNA_anno$strand == "+", ciRNA_anno$start, ciRNA_anno$end)
#ciRNA_anno$SS3 <- ifelse(ciRNA_anno$strand == "-", ciRNA_anno$start, ciRNA_anno$end) # after for cluster anno; add to new_anno
ciRNA_anno$width <- as.numeric(ciRNA_anno$width)

merged_sisRNA_anno <- ciRNA_anno %>% group_by(flankIntron, SS5) %>% summarise(name2 = name[which(sumReads == max(sumReads))[order(width[which(sumReads == max(sumReads))], decreasing = T)][1]],
                                                                  chr = chr[which(sumReads == max(sumReads))[order(width[which(sumReads == max(sumReads))], decreasing = T)][1]],
                                                                  start = start[which(sumReads == max(sumReads))[order(width[which(sumReads == max(sumReads))], decreasing = T)][1]],
                                                                  end = end[which(sumReads == max(sumReads))[order(width[which(sumReads == max(sumReads))], decreasing = T)][1]],
                                                                  strand = strand[which(sumReads == max(sumReads))[order(width[which(sumReads == max(sumReads))], decreasing = T)][1]],
                                                                  collapsed = paste(unique(name), collapse = ","),
                                                                  nr_collapsed = n(),
                                                                  #dist_3SS = max(SS3) - min(SS3), # after for cluster anno; add to new_anno
                                                                  geneName = paste(unique(geneName), collapse = ","),
                                                                  width = width[which(sumReads == max(sumReads))[order(width[which(sumReads == max(sumReads))], decreasing = T)][1]],
                                                                  sumReads = sum(sumReads),
                                                                  mean = NA,
                                                                  max = NA,
                                                                  sumExp = NA,
                                                                  median = NA,
                                                                  sumReads_lin = sum(sumReads_lin),
                                                                  gene_type = paste(unique(gene_type), collapse = ","),
                                                                  intron_id2 = intron_id[which(sumReads == max(sumReads))[order(width[which(sumReads == max(sumReads))], decreasing = T)][1]],
                                                                  intron_id_all = paste(unique(intron_id), collapse = ",")) %>% ungroup()

merged_sisRNA_anno <- unique.data.frame(merged_sisRNA_anno)
merged_sisRNA_anno <- merged_sisRNA_anno[,-2]


## exp matrices ##

# circular junction
ciRNA_exp_df <- as.data.frame(ciRNA_exp)
merged_sisRNA_exp <- lapply(1:length(ciRNA_anno$name), function(x) colSums(ciRNA_exp_df[which(ciRNA_anno$flankIntron == ciRNA_anno$flankIntron[x]),], na.rm = T))
merged_sisRNA_exp <- do.call(rbind, merged_sisRNA_exp)
#identical(merged_sisRNA_anno$sumReads, rowSums(merged_sisRNA_exp)) # check: T
merged_sisRNA_exp2 <- as.data.frame(merged_sisRNA_exp)
merged_sisRNA_exp2$flankIntron <- ciRNA_anno$flankIntron
merged_sisRNA_exp2 <- unique.data.frame(merged_sisRNA_exp2)
rownames(merged_sisRNA_exp2) <- merged_sisRNA_anno$name2[match(merged_sisRNA_exp2$flankIntron, merged_sisRNA_anno$flankIntron)]
merged_sisRNA_exp2 <- merged_sisRNA_exp2[,-458]
merged_sisRNA_exp2 <- merged_sisRNA_exp2[match(merged_sisRNA_anno$name2, rownames(merged_sisRNA_exp2)),]
identical(merged_sisRNA_anno$name2, rownames(merged_sisRNA_exp2))

# exp_norm
#length(which(sample_anno$UniqueID != colnames(ciRNA_exp_norm)))
merged_sisRNA_exp_norm = t(t(merged_sisRNA_exp2)/sample_anno$totalReads)*10^6


## Update sisRNA_anno ##
merged_sisRNA_anno$sumReads2 = apply(merged_sisRNA_exp2, 1, sum) # check if match
merged_sisRNA_anno$mean = apply(merged_sisRNA_exp2, 1, mean)
merged_sisRNA_anno$max = apply(merged_sisRNA_exp2, 1, max)
merged_sisRNA_anno$sumExp = apply(merged_sisRNA_exp2, 1, function(x) sum(x > 0))
merged_sisRNA_anno$median = apply(merged_sisRNA_exp2, 1, median)

identical(merged_sisRNA_anno$sumReads, merged_sisRNA_anno$sumReads2)
length(which(merged_sisRNA_anno$sumReads != merged_sisRNA_anno$sumReads2))
#length(merged_sisRNA_anno$name != rownames(merged_sisRNA_exp)) #0
length(which(merged_sisRNA_anno$name2 != rownames(merged_sisRNA_exp2)))

merged_sisRNA_anno <- merged_sisRNA_anno[,c(3:5, 2, 6,7,8,9,1,10:19)]
colnames(merged_sisRNA_anno)[c(4,18)] <- c("name", "intron_id")


############################
##### Save updated DFs #####

merged_ciRNA_exp <- merged_sisRNA_exp2
merged_ciRNA_anno <- merged_sisRNA_anno
merged_ciRNA_exp_norm <- merged_sisRNA_exp_norm

path <- paste(mount_path, "/output/path/200805_circexplorer2/", sep = "")
list_collapsed_df = c("list_collapsed_df", "merged_ciRNA_exp", "merged_ciRNA_anno", "merged_ciRNA_exp_norm","sample_anno")


# NMIBC
save(list = list_collapsed_df, file = paste0(path, "Updated_merged_ciRNAs_NMIBC_c.Rdata")) 

