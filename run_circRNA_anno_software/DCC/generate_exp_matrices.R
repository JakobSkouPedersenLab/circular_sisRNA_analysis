##############################################################################
##############################################################################
############# Make expression and annotation DFs for DCC output ##############
##############################################################################
##############################################################################
# Generated: 05/01/2023


###############################################
############### load material #################

# run on cluster using
# srun --mem 32G -t 6:00:00 --pty /bin/bash
# conda activate X
# cd ./path/script
# Rscript generate_exp_matrices.R

##### Load packages #####
library(tidyverse)

##### Load data #####
mount_path <- "/path"
path <- paste(mount_path, "/output/path/221222_DCC", sep = "")


###############################################
############# generate matrices ###############

## expression matrix ##
dirs <- list.files(path = path, recursive = TRUE, include.dirs = TRUE, pattern="CircCoordinates") # change according to output name
#dirs
length(dirs)

samples <- sub("/.*","",dirs)
#samples
samples <- sapply(stringr::str_split(samples, pattern = "_"), `[`, 1)
#length(samples)
#length(unique(samples))

# Only keep NMIBC samples
idx <- which(samples %in% sample_anno$UniqueID)
length(idx)
dirs <- dirs[idx]
samples <- samples[idx]

df <- read.table(paste0(path, "/", dirs[1]), as.is = TRUE, fill=TRUE, header=FALSE, skip =1)
col_names <- c("chr", "start", "end", "geneName", "junctionType","strand", "start_end_region", "region")
colnames(df) <- col_names
df$name = paste(df$chr, ":", df$start, "-", df$end, sep="")

for (i in 2:length(dirs)){
  df2<- read.table(paste0(path,"/",dirs[i]), as.is = TRUE, fill=TRUE, header=FALSE, skip =1)
  col_names <- c("chr", "start", "end", "geneName", "junctionType","strand", "start_end_region", "region")
  colnames(df2) <- col_names
  df2$name = paste(df2$chr, ":", df2$start, "-", df2$end, sep="")
  df <- merge(df, df2, by.x = c(col_names, "name"), by.y = c(col_names, "name"), all=TRUE)
}
circTable <- df
dim(circTable)

circTable <- circTable %>% filter(grepl("chr", chr)) %>% group_by(name, strand) %>% summarise(chr = chr,
                                                                                      start = start,
                                                                                      end = end,
                                                                                      #strand = strand,
                                                                                      geneName = paste(unique(geneName), collapse = ","),
                                                                                      junctionType = paste(unique(junctionType), collapse = ","),
                                                                                      start_end_region = paste(unique(start_end_region), collapse = ","),
                                                                                      region = paste(unique(unlist(stringr::str_split(region, pattern = ","))), collapse = ",")) %>% ungroup()

circTable <- unique.data.frame(circTable)
#dim(circTable)
#length(which(circTable$start_end_region == "intron-intron"))
#length(which(circTable$geneName == "HNRNPK"))

## count/expression table ##
# Only keep NMIBC samples
idx <- which(samples %in% sample_anno$UniqueID)
length(idx)
dirs <- dirs[idx]
samples <- samples[idx]

dirs2 <- list.files(path = path, recursive = TRUE, include.dirs = TRUE, pattern="CircRNACount")
dirs2 <- dirs2[idx]

df <- read.table(paste0(path, "/", dirs2[1]), as.is = TRUE, fill=TRUE, header=FALSE, skip =1)
col_names <- c("chr", "start", "end", samples[1])
colnames(df) <- col_names
df$name = paste(df$chr, ":", df$start, "-", df$end, sep="")
df <- df[, c(1:3,5,4)]

for (i in 2:length(dirs2)){ 
  df2<- read.table(paste0(path, "/",dirs2[i]), as.is = TRUE, fill=TRUE, header=FALSE, skip =1)
  col_names <- c("chr", "start", "end", samples[i])
  colnames(df2) <- col_names
  df2$name = paste(df2$chr, ":", df2$start, "-", df2$end, sep="")
  df2 <- df2[, c(1:3,5,4)]
  df <- merge(df, df2, by.x = c("chr", "start", "end", "name"), by.y = c("chr", "start", "end", "name"), all=TRUE)
}
circCounts <- df
circCounts[is.na(circCounts)] = 0
circCounts <- circCounts %>% filter(grepl("chr", chr))
dim(circCounts)

# DCC uses 0-based bed format. Add 1 to start coordinate
circCounts$start = circCounts$start+1
rownames(circCounts) = make.names(circCounts$name, unique=TRUE)

# filter circRNA only found in NMIBC (and not MIBC, CL)
circTable <- circTable %>% filter(name %in% circCounts$name)
idx <- which(duplicated(circTable$name))
idx2 <- which(circTable$name %in% circTable$name[idx])
circTable_dup <- circTable[idx2,]
circTable2 <- circTable[-idx2,]
circTable_dup <- circTable_dup %>% filter(geneName != "not_annotated" & region != "not_annotated") # filter odd annotations (from duplicated names)
circTable2 <- rbind(circTable2, circTable_dup)

circCounts2 <- circCounts %>% filter(name %in% circTable2$name)
#dim(circCounts2)
#dim(circTable2)

# combine remaining duplicated circRNA anno rows
circTable2 <- circTable2 %>% group_by(name) %>% summarise(chr = chr,
                                                          start = start,
                                                          end = end,
                                                          strand = paste(unique(strand), collapse = ";"),
                                                          geneName = paste(unique(geneName), collapse = ";"),
                                                          junctionType = paste(unique(junctionType), collapse = ";"),
                                                          start_end_region = paste(unique(start_end_region), collapse = ";"),
                                                          region = paste(unique(region), collapse = ";")) %>% ungroup()

circTable2 <- unique.data.frame(circTable2)
#dim(circTable2)

# rearrange and match row order
circTable2 <- circTable2[order(circTable2$chr, circTable2$start),]
circCounts2 <- circCounts2[order(circCounts2$chr, circCounts2$start),]
circCounts2 <- circCounts2[match(circTable2$name,circCounts2$name),]
identical(circCounts2$name, circTable2$name) # T

sample_anno <- sample_anno[match(colnames(circCounts2[5:461]), sample_anno$UniqueID),]
identical(as.character(sample_anno$UniqueID), colnames(circCounts2[5:461])) # T
#length(which(sample_anno2$UniqueID != colnames(circCounts2[5:461])))
rownames(circCounts2) <- circCounts2$name
circCounts2 <- circCounts2[5:461]

# Normalized exp matrix
circCounts_norm = t(t(circCounts2)/sample_anno$totalReads)*10^6

##### SAVE DATA ######
circCounts <- circCounts2
circTable <- circTable2

DCC_NMIBC = c("circTable", "circCounts", "circCounts_norm", "sample_anno")
save(list=DCC_NMIBC, file=paste0(path, "/DCC_NMIBC.Rdata"))




