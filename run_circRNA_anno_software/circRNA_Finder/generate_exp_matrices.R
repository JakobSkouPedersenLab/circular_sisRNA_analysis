#########################################################################################
#########################################################################################
############# Make expression and annotation DFs for circRNA_Finder output ##############
#########################################################################################
#########################################################################################
# Generated: 17/01/2023


###############################################
############### load material #################

# run on cluster using
# use bash script with same name


##### Load packages #####
library(tidyverse)
library(rtracklayer)

##### Load data #####
mount_path <- "/path"

path <-paste(mount_path, "/output/path", sep = "")

hg38_introns <- rtracklayer::import(paste(mount_path, "/path/hg38v33_ctat_introns_noOverlap_withGeneId.gtf", sep = ""))


###############################################
############# generate matrices ###############

## expression matrix ##
dirs <- list.files(path = path, recursive = TRUE, include.dirs = TRUE, pattern="filteredJunctions.bed") # change according to output name
dirs2 <- list.files(path = path, recursive = TRUE, include.dirs = TRUE, pattern="s_filteredJunctions.bed")
dirs3 <- list.files(path = path, recursive = TRUE, include.dirs = TRUE, pattern="s_filteredJunctions_fw.bed")
#dirs
length(dirs)
length(dirs2)
dirs <- dirs[-which(dirs %in% dirs2)]

#samples
samples <- sub("/.*","",dirs)

# Only keep NMIBC samples
idx <- which(samples %in% sample_anno$UniqueID)
length(idx)
dirs <- dirs[idx]
dirs2 <- dirs2[idx]
dirs3 <- dirs3[idx]
samples <- samples[idx]

df <- read.table(paste0(path, "/", dirs[1]), as.is = TRUE, fill=TRUE, header=FALSE) # , skip =1
df_s <- read.table(paste0(path, "/", dirs2[1]), as.is = TRUE, fill=TRUE, header=FALSE) 
df_fw <- read.table(paste0(path, "/", dirs3[1]), as.is = TRUE, fill=TRUE, header=FALSE)
col_names <- c("chr", "start", "end", "tmp", samples[1],"strand")
colnames(df) <- col_names
colnames(df_s) <- col_names
colnames(df_fw) <- col_names
df$name = paste(df$chr, ":", df$start, "-", df$end, sep="")
df_s$name = paste(df_s$chr, ":", df_s$start, "-", df_s$end, sep="")
df_fw$name = paste(df_fw$chr, ":", df_fw$start, "-", df_fw$end, sep="")
df <- df[,c(1:3,6,4,7,5)]
df_s <- df_s[,c(1:3,6,4,7,5)]
df_fw <- df_fw[,c(1:3,6,4,7,5)]

# Make circRNA GRanges
circRNA_ranges <- makeGRangesFromDataFrame(df, 
                                           seqnames.field = "chr",
                                           start.field = "start",
                                           end.field = "end",
                                           strand.field = "strand",
                                           keep.extra.columns = T)

hits <- findOverlaps(circRNA_ranges, hg38_introns, type = "any") # due to many circRNA calls; only keep ones overlapping intronic loci, to reduce matrix size
df <- df[queryHits(hits),]

circRNA_ranges <- makeGRangesFromDataFrame(df_s, 
                                           seqnames.field = "chr",
                                           start.field = "start",
                                           end.field = "end",
                                           strand.field = "strand",
                                           keep.extra.columns = T)

hits <- findOverlaps(circRNA_ranges, hg38_introns, type = "any") # due to many circRNA calls; only keep ones overlapping intronic loci, to reduce matrix size
df_s <- df_s[queryHits(hits),]

circRNA_ranges <- makeGRangesFromDataFrame(df_fw, 
                                           seqnames.field = "chr",
                                           start.field = "start",
                                           end.field = "end",
                                           strand.field = "strand",
                                           keep.extra.columns = T)

hits <- findOverlaps(circRNA_ranges, hg38_introns, type = "any") # due to many circRNA calls; only keep ones overlapping intronic loci, to reduce matrix size
df_fw <- df_fw[queryHits(hits),]

for (i in 2:length(dirs)){
  col_names <- c("chr", "start", "end", "tmp", samples[i],"strand")
  
  df2 <- read.table(paste0(path,"/",dirs[i]), as.is = TRUE, fill=TRUE, header=FALSE)
  colnames(df2) <- col_names
  df2$name = paste(df2$chr, ":", df2$start, "-", df2$end, sep="")
  df2 <- df2[,c(1:3,6,4,7,5)]
  
  # Make circRNA GRanges
  circRNA_ranges <- makeGRangesFromDataFrame(df2, 
                                             seqnames.field = "chr",
                                             start.field = "start",
                                             end.field = "end",
                                             strand.field = "strand",
                                             keep.extra.columns = T)
  
  hits <- findOverlaps(circRNA_ranges, hg38_introns, type = "any") # due to many circRNA calls; only keep ones overlapping intronic loci, to reduce matrix size
  df2 <- df2[queryHits(hits),]
  df <- merge(df, df2, by.x = c("chr", "start", "end", "strand", "tmp", "name"), by.y = c("chr", "start", "end", "strand", "tmp", "name"), all=TRUE)
  
  #______________________________________________________________________________________________________
  
  df_s2 <- read.table(paste0(path,"/",dirs2[i]), as.is = TRUE, fill=TRUE, header=FALSE)
  colnames(df_s2) <- col_names
  df_s2$name = paste(df_s2$chr, ":", df_s2$start, "-", df_s2$end, sep="")
  df_s2 <- df_s2[,c(1:3,6,4,7,5)]
  
  # Make circRNA GRanges
  circRNA_ranges <- makeGRangesFromDataFrame(df_s2, 
                                             seqnames.field = "chr",
                                             start.field = "start",
                                             end.field = "end",
                                             strand.field = "strand",
                                             keep.extra.columns = T)
  
  hits <- findOverlaps(circRNA_ranges, hg38_introns, type = "any") # due to many circRNA calls; only keep ones overlapping intronic loci, to reduce matrix size
  df_s2 <- df_s2[queryHits(hits),]
  df_s <- merge(df_s, df_s2, by.x = c("chr", "start", "end", "strand", "tmp", "name"), by.y = c("chr", "start", "end", "strand", "tmp", "name"), all=TRUE)
  
  #______________________________________________________________________________________________________
  
  df_fw2 <- read.table(paste0(path,"/",dirs3[i]), as.is = TRUE, fill=TRUE, header=FALSE)
  colnames(df_fw2) <- col_names
  df_fw2$name = paste(df_fw2$chr, ":", df_fw2$start, "-", df_fw2$end, sep="")
  df_fw2 <- df_fw2[,c(1:3,6,4,7,5)]
  
  # Make circRNA GRanges
  circRNA_ranges <- makeGRangesFromDataFrame(df_fw2, 
                                             seqnames.field = "chr",
                                             start.field = "start",
                                             end.field = "end",
                                             strand.field = "strand",
                                             keep.extra.columns = T)
  
  hits <- findOverlaps(circRNA_ranges, hg38_introns, type = "any") # due to many circRNA calls; only keep ones overlapping intronic loci, to reduce matrix size
  df_fw2 <- df_fw2[queryHits(hits),]
  df_fw <- merge(df_fw, df_fw2, by.x = c("chr", "start", "end", "strand", "tmp", "name"), by.y = c("chr", "start", "end", "strand", "tmp", "name"), all=TRUE)
}
circTable <- df[,c(1:6)]
circCounts <- df # count/expression table - full
circCounts_s <- df_s # count/expression table - GT/AG SS
circCounts_fw <- df_fw # count/expression - foward junction
dim(circTable) 
unique(circTable$chr)

# filter(grepl("chr", chr)) %>%
circTable <- circTable %>%  group_by(chr, start, end,  strand, name) %>% summarise(tmp = paste(unique(tmp), collapse = ",")) %>% ungroup() # , width = width
circTable$start = circTable$start+1
circTable$width <- circTable$end - circTable$start
circTable <- circTable %>% filter(width > 50, width < 10000)
circTable <- unique.data.frame(circTable)
dim(circTable)
circTable <- circTable[,c(1:4,6,5,7)] # BE CAREUL WITH COLUMN IDX

# fix unresolved strands
idx <- which(duplicated(circTable$name))
idx2 <- which(circTable$name %in% circTable$name[idx])
if (length(idx2) > 0){
  circTable_dup <- circTable[idx2,]
  circTable <- circTable[-idx2,]
  circTable_dup$strand <- "*" 
  circTable_dup <- circTable_dup %>%  group_by(chr, start, end,  strand, name, width) %>% summarise(tmp = paste(unique(tmp), collapse = ",")) %>% ungroup()
  circTable_dup <- unique.data.frame(circTable_dup)
  circTable_dup <- circTable_dup[,c(1:4,7,5,6)] # BE CAREUL WITH COLUMN IDX
  circTable <- rbind(circTable, circTable_dup) 
}

# circRNA_Finder uses 0-based bed format. Add 1 to start coordinate
circCounts$start = circCounts$start+1
circCounts_s$start = circCounts_s$start+1
circCounts_fw$start = circCounts_fw$start+1
rownames(circCounts) = make.names(circCounts$name, unique=TRUE)
rownames(circCounts_s) = make.names(circCounts_s$name, unique=TRUE)
rownames(circCounts_fw) = make.names(circCounts_fw$name, unique=TRUE)

# order table
circTable <- circTable[order(circTable$chr, circTable$start),]

# filter circRNA only found in NMIBC & filtered circTable
circCounts <- circCounts %>% filter(grepl("chr", chr))
circCounts <- circCounts %>% filter(name %in% circTable$name)
circCounts[is.na(circCounts)] = 0

circCounts_s <- circCounts_s %>% filter(grepl("chr", chr))
circCounts_s <- circCounts_s %>% filter(name %in% circTable$name)
circCounts_s[is.na(circCounts_s)] = 0

circCounts_fw <- circCounts_fw %>% filter(grepl("chr", chr))
circCounts_fw <- circCounts_fw %>% filter(name %in% circTable$name)
circCounts_fw[is.na(circCounts_fw)] = 0

# combine duplicated circRNA anno rows for the count tables
circCounts2 <- lapply(1:length(circTable$name), function(x) colSums(circCounts[which(circCounts$name == circTable$name[x]),-c(1:6)], na.rm = T))
circCounts2 <- do.call(rbind, circCounts2)
circCounts2 <- as.data.frame(circCounts2)
rownames(circCounts2) <- circTable$name
identical(circTable$name, rownames(circCounts2)) # T

m_names <- unique(circCounts_s$name)
circCounts_s2 <- lapply(1:length(m_names), function(x) colSums(circCounts_s[which(circCounts_s$name == m_names[x]),-c(1:6)], na.rm = T))
circCounts_s2 <- do.call(rbind, circCounts_s2)
circCounts_s2 <- as.data.frame(circCounts_s2)
rownames(circCounts_s2) <- m_names

m_names2 <- unique(circCounts_fw$name)
circCounts_fw2 <- lapply(1:length(m_names2), function(x) colSums(circCounts_fw[which(circCounts_fw$name == m_names2[x]),-c(1:6)], na.rm = T))
circCounts_fw2 <- do.call(rbind, circCounts_fw2)
circCounts_fw2 <- as.data.frame(circCounts_fw2)
rownames(circCounts_fw2) <- m_names2
circCounts_fw2 <- circCounts_fw2[rownames(circCounts_s2),]

# get sample anno
sample_anno <- sample_anno[match(samples, sample_anno$UniqueID),]
identical(as.character(sample_anno$UniqueID), samples) # T
identical(as.character(sample_anno$UniqueID), colnames(circCounts_fw2))

# Normalized exp matrix
circCounts_norm = t(t(circCounts2)/sample_anno$totalReads)*10^6
circCounts_s_norm = t(t(circCounts_s2)/sample_anno$totalReads)*10^6
circCounts_fw_norm = t(t(circCounts_fw2)/sample_anno$totalReads)*10^6

##### SAVE DATA ######
circCounts <- circCounts2
circCounts_s <- circCounts_s2
circCounts_fw <- circCounts_fw2


circRNA_Finder_NMIBC = c("circTable", "circCounts", "circCounts_s", "circCounts_fw", "circCounts_norm", "circCounts_s_norm", "circCounts_fw_norm", "sample_anno")
save(list=circRNA_Finder_NMIBC, file=paste0(path, "/circRNA_Finder_NMIBC.Rdata"))




