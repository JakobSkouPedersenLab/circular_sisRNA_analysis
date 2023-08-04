##########################################################################################
################ match circRNA_Finder circRNAs to non-overlapping introns ################
##########################################################################################
# 17/01/2023

# Need to properly annotate circRNAs to identify sisRNA subset.

##### Load data and packages #####
library(tidyverse)
library(rtracklayer)

##### Load general data #####

# User mount path:
mount_path <- "/path"
path <- paste(mount_path, "/output/path", sep = "")

# circRNA_Finder
load(paste(mount_path, "/path/circRNA_Finder_NMIBC.Rdata", sep = ""))
hg38_introns <- rtracklayer::import(paste(mount_path, "/path/hg38v33_ctat_introns_noOverlap_withGeneId.gtf", sep = ""))


#############################################################
############# collapse circRNA sharing 5' end ###############

# add basic annotations
circTable$sumReads = apply(circCounts, 1, sum)
circTable$mean = apply(circCounts, 1, mean)
circTable$max = apply(circCounts, 1, max)
circTable$sumExp = apply(circCounts, 1, function(x) sum(x > 0))
circTable$median = apply(circCounts, 1, median)

# width
circTable$width <- circTable$end - circTable$start +1

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

hits <- findOverlaps(circRNA_ranges, hg38_introns, type = "within") # strict definition
intron_ranges_df <- as.data.frame(hg38_introns)

intron_ranges_df <- intron_ranges_df[subjectHits(hits),] # Match rows
circRNA_anno <- circTable
circRNA_anno <- circRNA_anno[queryHits(hits),] 
#dim(intron_ranges_df)
#dim(circRNA_anno)

#____________
# Reestablish correct intron coordinates
intron_ranges_df$start <- intron_ranges_df$start + 10
intron_ranges_df$end <- intron_ranges_df$end - 10
#____________

dist1 <- abs(circRNA_anno$start - intron_ranges_df$start)
dist2 <- abs(circRNA_anno$end - intron_ranges_df$end)

pot_sisRNA <- ifelse(dist1 <= 10 | dist2 <= 10, T, F) # 5' SS should be within 10 bp for introns to match.
length(which(pot_sisRNA == T)) 

circRNA_anno$pot_sisRNA <- pot_sisRNA # sisRNAs where we can match introns within 10 bp of 5' splice site.
intron_ranges_df$intron_name <- paste0(intron_ranges_df$seqnames, sep = ":", intron_ranges_df$start, sep = "-", intron_ranges_df$end)
#____
intron_ranges_df <- intron_ranges_df[,-11] # Be careful with column idx!
#____

sisRNA_df <- cbind(circRNA_anno, intron_ranges_df[,c(6,7,8,9,10)]) # Match intron to ciRNA # Be careful with column idx!

intron_ranges_df$name_sisRNA <- sisRNA_df$name
dim(intron_ranges_df)

#length(unique(intron_ranges_df$name_sisRNA)) 
#length(unique(sisRNA_df$intron_name)) 

#dim(sisRNA_df[duplicated(sisRNA_df$name),]) # NMIBC: 8 ciRNAs overlap on both strands
idx <- which(duplicated(sisRNA_df$name))
idx2 <- which(sisRNA_df$name %in% sisRNA_df$name[idx])

if (length(idx2) > 0){
  sisRNA_df_dup <- sisRNA_df[idx2,]
  sisRNA_df <- sisRNA_df[-idx2,]
  sisRNA_df_dup <- sisRNA_df_dup %>%  group_by(chr, start, end,  strand, tmp, name, width, sumReads, mean, max, sumExp, median, pot_sisRNA) %>% 
                                      summarise(gene_id = paste(unique(gene_id), collapse = ","), gene_type = paste(unique(gene_type), collapse = ","), gene_name = paste(unique(gene_name), collapse = ","),
                                                intron_id = paste(unique(intron_id), collapse = ","), intron_name = paste(unique(intron_name), collapse = ",")) %>% ungroup()
  sisRNA_df_dup <- unique.data.frame(sisRNA_df_dup)
  sisRNA_df <- rbind(sisRNA_df, sisRNA_df_dup) 
}

## Collapse ciRNA data ##

# Collapse based on intron cluster definition; note that circRNA_Finder generates many sisRNA candidates not matching canonical splicing!

#length(unique(sisRNA_df$intron_name))
merged_sisRNA_df <- sisRNA_df %>% group_by(intron_name) %>% summarise(name = name[which(sumReads == max(sumReads))[1]],
                                                                      chr = chr[which(sumReads == max(sumReads))[1]],
                                                                      start = start[which(sumReads == max(sumReads))[1]],
                                                                      end = end[which(sumReads == max(sumReads))[1]],
                                                                      strand = strand[which(sumReads == max(sumReads))[1]],
                                                                      collapsed = paste(unique(name), collapse = ","),
                                                                      geneName = paste(unique(gene_name), collapse = "-"),
                                                                      junctionType = paste(unique(tmp), collapse = "-"),
                                                                      gene_id = paste(unique(gene_id), collapse = "-"),
                                                                      gene_type = paste(unique(gene_type), collapse = "-"),
                                                                      sumReads = sum(sumReads),
                                                                      mean = NA,
                                                                      max = NA,
                                                                      sumExp = NA,
                                                                      median = NA,
                                                                      width = width[which(sumReads == max(sumReads))[1]],
                                                                      pot_sisRNA = (T %in% pot_sisRNA),
                                                                      intron_id = intron_id[which(sumReads == max(sumReads))[1]]) %>% ungroup()
merged_sisRNA_df <- unique.data.frame(merged_sisRNA_df)
#dim(merged_sisRNA_df)
#length(unique(merged_sisRNA_df$intron_name))
#length(which(merged_sisRNA_df$pot_sisRNA == T))

## exp matrices ##
sisRNA_introns <- unique(sisRNA_df$intron_name)
length(sisRNA_introns)

sisRNA_exp <- circCounts[match(sisRNA_df$name, rownames(circCounts)),]
sisRNA_exp_norm <- as.data.frame(circCounts_norm)
sisRNA_exp_norm <- sisRNA_exp_norm[match(sisRNA_df$name, rownames(sisRNA_exp_norm)),]
##identical(rownames(sisRNA_exp), sisRNA_df$name) # T
##length(which(rownames(sisRNA_exp) != sisRNA_df$name)) # 0
##length(which(rownames(sisRNA_exp_norm) != sisRNA_df$name)) # 0

sisRNA_s_exp <- circCounts_s[which(rownames(circCounts_s) %in% sisRNA_df$name),]
sisRNA_s_exp_norm <- as.data.frame(circCounts_s_norm)
sisRNA_s_exp_norm <- sisRNA_s_exp_norm[match(rownames(sisRNA_s_exp), rownames(sisRNA_s_exp_norm)),]
identical(rownames(sisRNA_s_exp_norm), rownames(sisRNA_s_exp))

sisRNA_fw_exp <- circCounts_fw[match(rownames(sisRNA_s_exp), rownames(circCounts_fw)),]
sisRNA_fw_exp_norm <- as.data.frame(circCounts_fw_norm)
sisRNA_fw_exp_norm <- sisRNA_fw_exp_norm[match(rownames(sisRNA_fw_exp), rownames(sisRNA_fw_exp_norm)),]
identical(rownames(sisRNA_s_exp_norm), rownames(sisRNA_s_exp))

# merged
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

#dim(merged_sisRNA_exp2)
#dim(merged_sisRNA_exp_norm2)

# merged _s & _fw
merged_Inames <- rownames(sisRNA_s_exp) # needed for _fw, whereas _sissubset of full count set
anno_tmp <- sisRNA_df[match(merged_Inames, sisRNA_df$name),]
merged_sisRNA_df_s_fw <- anno_tmp %>% group_by(intron_name) %>% summarise(name = name[which(sumReads == max(sumReads))[1]],
                                                                      chr = chr[which(sumReads == max(sumReads))[1]],
                                                                      start = start[which(sumReads == max(sumReads))[1]],
                                                                      end = end[which(sumReads == max(sumReads))[1]],
                                                                      strand = strand[which(sumReads == max(sumReads))[1]],
                                                                      collapsed = paste(unique(name), collapse = ","),
                                                                      geneName = paste(unique(gene_name), collapse = "-"),
                                                                      junctionType = paste(unique(tmp), collapse = "-"),
                                                                      gene_id = paste(unique(gene_id), collapse = "-"),
                                                                      gene_type = paste(unique(gene_type), collapse = "-"),
                                                                      sumReads = sum(sumReads),
                                                                      mean = NA,
                                                                      max = NA,
                                                                      sumExp = NA,
                                                                      median = NA,
                                                                      width = width[which(sumReads == max(sumReads))[1]],
                                                                      pot_sisRNA = (T %in% pot_sisRNA),
                                                                      intron_id = intron_id[which(sumReads == max(sumReads))[1]]) %>% ungroup()
merged_sisRNA_df_s_fw <- unique.data.frame(merged_sisRNA_df_s_fw)
#dim(merged_sisRNA_df_s_fw)

merged_sisRNA_s_exp <- lapply(1:length(anno_tmp$name), function(x) colSums(sisRNA_s_exp[which(anno_tmp$intron_name == anno_tmp$intron_name[x]),], na.rm = T))
merged_sisRNA_s_exp <- do.call(rbind, merged_sisRNA_s_exp)
merged_sisRNA_s_exp2 <- as.data.frame(merged_sisRNA_s_exp)
merged_sisRNA_s_exp2$intron_name <- anno_tmp$intron_name
merged_sisRNA_s_exp2 <- unique.data.frame(merged_sisRNA_s_exp2)
rownames(merged_sisRNA_s_exp2) <- merged_sisRNA_df_s_fw$name[match(merged_sisRNA_s_exp2$intron_name, merged_sisRNA_df_s_fw$intron_name)]
merged_sisRNA_s_exp2 <- merged_sisRNA_s_exp2[,-458]
merged_sisRNA_s_exp2 <- merged_sisRNA_s_exp2[match(merged_sisRNA_df_s_fw$name, rownames(merged_sisRNA_s_exp2)),]
identical(merged_sisRNA_df_s_fw$name, rownames(merged_sisRNA_s_exp2))

merged_sisRNA_fw_exp <- lapply(1:length(anno_tmp$name), function(x) colSums(sisRNA_fw_exp[which(anno_tmp$intron_name == anno_tmp$intron_name[x]),], na.rm = T))
merged_sisRNA_fw_exp <- do.call(rbind, merged_sisRNA_fw_exp)
merged_sisRNA_fw_exp2 <- as.data.frame(merged_sisRNA_fw_exp)
merged_sisRNA_fw_exp2$intron_name <- anno_tmp$intron_name
merged_sisRNA_fw_exp2 <- unique.data.frame(merged_sisRNA_fw_exp2)
rownames(merged_sisRNA_fw_exp2) <- merged_sisRNA_df_s_fw$name[match(merged_sisRNA_fw_exp2$intron_name, merged_sisRNA_df_s_fw$intron_name)]
merged_sisRNA_fw_exp2 <- merged_sisRNA_fw_exp2[,-458]
merged_sisRNA_fw_exp2 <- merged_sisRNA_fw_exp2[match(merged_sisRNA_df_s_fw$name, rownames(merged_sisRNA_fw_exp2)),]
identical(merged_sisRNA_df_s_fw$name, rownames(merged_sisRNA_fw_exp2))

#dim(merged_sisRNA_s_exp2)
#dim(merged_sisRNA_fw_exp2)

merged_sisRNA_s_exp_norm <- lapply(1:length(anno_tmp$name), function(x) colSums(sisRNA_s_exp_norm[which(anno_tmp$intron_name == anno_tmp$intron_name[x]),], na.rm = T))
merged_sisRNA_s_exp_norm <- do.call(rbind, merged_sisRNA_s_exp_norm)
merged_sisRNA_s_exp_norm2 <- as.data.frame(merged_sisRNA_s_exp_norm)
merged_sisRNA_s_exp_norm2$intron_name <- anno_tmp$intron_name
merged_sisRNA_s_exp_norm2 <- unique.data.frame(merged_sisRNA_s_exp_norm2)
rownames(merged_sisRNA_s_exp_norm2) <- merged_sisRNA_df_s_fw$name[match(merged_sisRNA_s_exp_norm2$intron_name, merged_sisRNA_df_s_fw$intron_name)]
merged_sisRNA_s_exp_norm2 <- merged_sisRNA_s_exp_norm2[,-458]
merged_sisRNA_s_exp_norm2 <- merged_sisRNA_s_exp_norm2[match(merged_sisRNA_df_s_fw$name, rownames(merged_sisRNA_s_exp_norm2)),]
identical(merged_sisRNA_df_s_fw$name, rownames(merged_sisRNA_s_exp_norm2))

merged_sisRNA_fw_exp_norm <- lapply(1:length(anno_tmp$name), function(x) colSums(sisRNA_fw_exp_norm[which(anno_tmp$intron_name == anno_tmp$intron_name[x]),], na.rm = T))
merged_sisRNA_fw_exp_norm <- do.call(rbind, merged_sisRNA_fw_exp_norm)
merged_sisRNA_fw_exp_norm2 <- as.data.frame(merged_sisRNA_fw_exp_norm)
merged_sisRNA_fw_exp_norm2$intron_name <- anno_tmp$intron_name
merged_sisRNA_fw_exp_norm2 <- unique.data.frame(merged_sisRNA_fw_exp_norm2)
rownames(merged_sisRNA_fw_exp_norm2) <- merged_sisRNA_df_s_fw$name[match(merged_sisRNA_fw_exp_norm2$intron_name, merged_sisRNA_df_s_fw$intron_name)]
merged_sisRNA_fw_exp_norm2 <- merged_sisRNA_fw_exp_norm2[,-458]
merged_sisRNA_fw_exp_norm2 <- merged_sisRNA_fw_exp_norm2[match(merged_sisRNA_df_s_fw$name, rownames(merged_sisRNA_fw_exp_norm2)),]
identical(merged_sisRNA_df_s_fw$name, rownames(merged_sisRNA_fw_exp_norm2))

## anno ##
merged_sisRNA_df$sumReads2 = apply(merged_sisRNA_exp2, 1, sum) # check if match
merged_sisRNA_df$mean = apply(merged_sisRNA_exp2, 1, mean)
merged_sisRNA_df$max = apply(merged_sisRNA_exp2, 1, max)
merged_sisRNA_df$sumExp = apply(merged_sisRNA_exp2, 1, function(x) sum(x > 0))
merged_sisRNA_df$median = apply(merged_sisRNA_exp2, 1, median)

identical(merged_sisRNA_df$sumReads, merged_sisRNA_df$sumReads2)
length(which(merged_sisRNA_df$sumReads != merged_sisRNA_df$sumReads2))
#length(merged_sisRNA_df$name != rownames(merged_sisRNA_exp)) #0
length(which(merged_sisRNA_df$name != rownames(merged_sisRNA_exp2)))

# anno _s_fw
merged_sisRNA_df_s_fw$sumReads2 = apply(merged_sisRNA_s_exp2, 1, sum) # check if match
merged_sisRNA_df_s_fw$mean = apply(merged_sisRNA_s_exp2, 1, mean)
merged_sisRNA_df_s_fw$max = apply(merged_sisRNA_s_exp2, 1, max)
merged_sisRNA_df_s_fw$sumExp = apply(merged_sisRNA_s_exp2, 1, function(x) sum(x > 0))
merged_sisRNA_df_s_fw$median = apply(merged_sisRNA_s_exp2, 1, median)
merged_sisRNA_df_s_fw$sumReads_fw = apply(merged_sisRNA_fw_exp2, 1, sum)

identical(merged_sisRNA_df_s_fw$sumReads, merged_sisRNA_df_s_fw$sumReads2)
length(which(merged_sisRNA_df_s_fw$sumReads != merged_sisRNA_df_s_fw$sumReads2))
#length(merged_sisRNA_df_s_fw$name != rownames(merged_sisRNA_s_exp)) #0
length(which(merged_sisRNA_df_s_fw$name != rownames(merged_sisRNA_s_exp2)))


#######################################
############# save data ###############

sisRNA_anno <- sisRNA_df
sisRNA_anno_s_fw <- merged_sisRNA_df_s_fw
merged_sisRNA_anno <- merged_sisRNA_df
intron_anno <- intron_ranges_df
merged_sisRNA_exp <- merged_sisRNA_exp2
merged_sisRNA_exp_norm <- merged_sisRNA_exp_norm2
merged_sisRNA_s_exp <- merged_sisRNA_s_exp2
merged_sisRNA_s_exp_norm <- merged_sisRNA_s_exp_norm2
merged_sisRNA_fw_exp <- merged_sisRNA_fw_exp2
merged_sisRNA_fw_exp_norm <- merged_sisRNA_fw_exp_norm2


circRNA_Finder_NMIBC = c("sisRNA_anno", "sisRNA_anno_s_fw", "sisRNA_exp", "sisRNA_exp_norm", "sisRNA_s_exp", "sisRNA_s_exp_norm", "sisRNA_fw_exp", "sisRNA_fw_exp_norm", 
              "merged_sisRNA_anno", "merged_sisRNA_exp", "merged_sisRNA_exp_norm", "merged_sisRNA_s_exp", "merged_sisRNA_s_exp_norm", 
              "merged_sisRNA_fw_exp", "merged_sisRNA_fw_exp_norm", "sample_anno", "intron_anno")

save(list=circRNA_Finder_NMIBC, file=paste0(path, "/Updated_circRNA_Finder_NMIBC.Rdata"))

#save(list=circRNA_Finder_NMIBC, file=paste0(path, "/Updated_circRNA_Finder_NMIBC_anyIntronOverlap.Rdata"))


