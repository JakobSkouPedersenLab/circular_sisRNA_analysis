#########################################################################################
#########################################################################################
############# Make expression and annotation DFs for CE2-CIRIquant output ###############
#########################################################################################
#########################################################################################
# Generated: 23/01/2023


###############################################
############### load material #################

##### Load packages #####
library(tidyverse)
library(rtracklayer)

##### Load data #####
mount_path <- "/Users/path/to/mount"
path <- paste(mount_path, "/path/to/outputFolder/out_CIRIquant_circexplorer2", sep = "")

sample_ids <- read.table(paste(mount_path, "/path/input.lst", sep = ""))


###########################################################
############# gather output for each sample ###############

dirs <- list.files(path = path, recursive = TRUE, include.dirs = TRUE, pattern=".gtf")
idx <- grepl("out|cov", dirs) # remove elements containing cov & out
dirs <- dirs[-which(idx)]

samples <- unique(sample_ids$V1)
samples <- samples[which(samples %in% sample_anno$UniqueID)]

dirs_new <- paste(samples, "/", samples, ".gtf", sep = "")
idx <- which(dirs %in% dirs_new) # remove updated data
dirs <- dirs[-idx]

# check for missing samples!
sample_ids$V5 <- paste(sample_ids$V1, "/", sample_ids$V1, ".", sample_ids$V4, ".gtf", sep ="")
idx <- which(sample_ids$V5 %in% dirs)
sample_ids$V5[-idx]

for (s in samples[1:length(samples)]){ 
  dirs_s <- dirs[which(grepl(s, dirs))]
  s_lst <- list()
  
  for (i in 1:length(dirs_s)){
    tmp <- rtracklayer::import(paste(path, "/", dirs_s[i], sep = ""))
    tmp <- as.data.frame(tmp)
    s_lst[[i]] <- tmp
  }
  
  #which(samples == s) # can't load empty gtf files
  
  s_df <- do.call(rbind, s_lst)
  s_df <- s_df %>% filter(circ_type != "exon")
  colnames(s_df) <- c("chr", "start", "end", "width", "strand", "source", "type", "CPM", "phase", "circ_id", "circ_type", "count", "count_fw", "junc_ratio", "gene_id", "gene_name", "gene_type")
  s_df <- s_df %>% group_by(chr, start, end, strand, circ_id, width) %>% 
    summarise(type = paste(unique(type), collapse = ","),
              CPM = sum(CPM),
              circ_type = paste(unique(circ_type), collapse = ","),
              count = sum(as.numeric(count)),
              count_fw = sum(as.numeric(count_fw)),
              junc_ratio = mean(as.numeric(junc_ratio)), # recalculate as 2 * bsj / ( 2 * bsj + fsj)
              gene_id = paste(unique(gene_id), collapse = ","),
              gene_name = paste(unique(gene_name), collapse = ","),
              gene_type = paste(unique(gene_type), collapse = ",")) %>% ungroup()
  
  s_df <- unique.data.frame(s_df)
  
  if(dim(s_df)[1] == 0){
    s_df[1,] <- NA
    s_df <- as.matrix(s_df)
    s_df[1,] <- c("chr0", 000, 111, "*", "id", 200, "dummy", 1, "dummy", 1, 1, 1, "dummy", "dummy", "dummy")
    s_df <- as.data.frame(s_df)
  } 
  rtracklayer::export(s_df, con = paste(path, "/", s, "/", s, ".gtf", sep = ""), format = "gtf")

}


###############################################
############# generate matrices ###############

## expression matrix ##

#dirs
length(dirs_new)

# generate exp DF
df <- as.data.frame(rtracklayer::import(paste0(path, "/", dirs_new[1])))
df <- df[,-c(8,9)]
colnames(df)[1] <- "chr"
df$name = paste(df$chr, ":", df$start, "-", df$end, sep="")
df <- df[,c(1:3,5,17, 4, 6:8, 10, 14:16, 9, 11, 12, 13)] # Be careful with indexing columns

df_CPM <- df[,-c(15,16,17)] # CPM
df_fw <- df[,-c(14,15,17)] # count_fw
df_ratio <- df[,-c(14,15,16)] # junc_ratio
df_count <- df[,-c(14,16,17)] # count

colnames(df_CPM)[14] <- samples[1]
colnames(df_fw)[14] <- samples[1]
colnames(df_ratio)[14] <- samples[1]
colnames(df_count)[14] <- samples[1]

for (i in 2:length(dirs_new)){ 
  
  df2 <- as.data.frame(rtracklayer::import(paste0(path, "/", dirs_new[i])))
  df2 <- df2[,-c(8,9)]
  colnames(df2)[1] <- "chr"
  df2$name = paste(df2$chr, ":", df2$start, "-", df2$end, sep="")
  df2 <- df2[,c(1:3,5,17, 4, 6:8, 10, 14:16, 9, 11, 12, 13)]
  
  df_CPM2 <- df2[,-c(15,16,17)] # CPM
  df_fw2 <- df2[,-c(14,15,17)] # count_fw
  df_ratio2 <- df2[,-c(14,15,16)] # junc_ratio
  df_count2 <- df2[,-c(14,16,17)] # count
  
  colnames(df_CPM2)[14] <- samples[i]
  colnames(df_fw2)[14] <- samples[i]
  colnames(df_ratio2)[14] <- samples[i]
  colnames(df_count2)[14] <- samples[i]
  
  df_CPM <- merge(df_CPM, df_CPM2, by.x = c("chr", "start", "end", "strand", "name", "width", "source", "type", "circ_id", "circ_type", "gene_id", "gene_name", "gene_type"), 
                  by.y = c("chr", "start", "end", "strand", "name", "width", "source", "type", "circ_id", "circ_type", "gene_id", "gene_name", "gene_type"), all=TRUE)
  df_fw <- merge(df_fw, df_fw2, by.x = c("chr", "start", "end", "strand", "name", "width", "source", "type", "circ_id", "circ_type", "gene_id", "gene_name", "gene_type"), 
                  by.y = c("chr", "start", "end", "strand", "name", "width", "source", "type", "circ_id", "circ_type", "gene_id", "gene_name", "gene_type"), all=TRUE)
  df_ratio <- merge(df_ratio, df_ratio2, by.x = c("chr", "start", "end", "strand", "name", "width", "source", "type", "circ_id", "circ_type", "gene_id", "gene_name", "gene_type"), 
                  by.y = c("chr", "start", "end", "strand", "name", "width", "source", "type", "circ_id", "circ_type", "gene_id", "gene_name", "gene_type"), all=TRUE)
  df_count <- merge(df_count, df_count2, by.x = c("chr", "start", "end", "strand", "name", "width", "source", "type", "circ_id", "circ_type", "gene_id", "gene_name", "gene_type"), 
                  by.y = c("chr", "start", "end", "strand", "name", "width", "source", "type", "circ_id", "circ_type", "gene_id", "gene_name", "gene_type"), all=TRUE)
  
  
}
circTable <- df_CPM[,c(1:13)]
circCounts <- df_count # count/expression table - intron
circCounts_norm <- df_CPM # count/expression table - CPM
circCounts_fw <- df_fw # count/expression - foward junction
circCounts_ratio <- df_ratio # count/expression - circ/lin junction ratio
dim(circTable) 
unique(circTable$chr)




circTable <- circTable %>% filter(chr != "chr0") 
length(unique(circTable$name))
#circTable$start = circTable$start+1 # is already one-based when importing gtf file.
circTable <- unique.data.frame(circTable)
dim(circTable)


# row names
rownames(circCounts) = circCounts$name #make.names(circCounts$name, unique=TRUE)
rownames(circCounts_norm) = circCounts_norm$name #make.names(circCounts_norm$name, unique=TRUE)
rownames(circCounts_fw) = circCounts_fw$name #make.names(circCounts_fw$name, unique=TRUE)
rownames(circCounts_ratio) = circCounts_ratio$name #make.names(circCounts_ratio$name, unique=TRUE)

# order table
circTable <- circTable[order(circTable$chr, circTable$start),]

# filter circRNA only found in NMIBC & filtered circTable
circCounts <- circCounts %>% filter(chr != "chr0")
circCounts <- circCounts %>% filter(name %in% circTable$name)
circCounts[is.na(circCounts)] = 0
circCounts <- circCounts[match(circTable$name, circCounts$name),]
identical(circTable$name,  rownames(circCounts)) 
circCounts <- circCounts[,-c(1:13)]

circCounts_norm <- circCounts_norm %>% filter(chr != "chr0")
circCounts_norm <- circCounts_norm %>% filter(name %in% circTable$name)
circCounts_norm[is.na(circCounts_norm)] = 0
circCounts_norm <- circCounts_norm[match(circTable$name, circCounts_norm$name),]
identical(circTable$name,  rownames(circCounts_norm)) 
circCounts_norm <- circCounts_norm[,-c(1:13)]

circCounts_fw <- circCounts_fw %>% filter(chr != "chr0")
circCounts_fw <- circCounts_fw %>% filter(name %in% circTable$name)
circCounts_fw[is.na(circCounts_fw)] = 0
circCounts_fw <- circCounts_fw[match(circTable$name, circCounts_fw$name),]
identical(circTable$name,  rownames(circCounts_fw)) 
circCounts_fw <- circCounts_fw[,-c(1:13)]

circCounts_ratio <- circCounts_ratio %>% filter(chr != "chr0")
circCounts_ratio <- circCounts_ratio %>% filter(name %in% circTable$name)
circCounts_ratio[is.na(circCounts_ratio)] = 0
circCounts_ratio <- circCounts_ratio[match(circTable$name, circCounts_ratio$name),]
identical(circTable$name,  rownames(circCounts_ratio)) 
circCounts_ratio <- circCounts_ratio[,-c(1:13)]

identical(colnames(circCounts), colnames(circCounts_fw)) # T

# get sample anno
sample_anno <- sample_anno[match(samples, sample_anno$UniqueID),]
identical(as.character(sample_anno$UniqueID), samples) # T
identical(as.character(sample_anno$UniqueID), colnames(circCounts_norm)) # T

# Normalized exp matrix
circCounts2 <- apply(circCounts, 2, as.numeric) #as.numeric(as.matrix(circCounts))
circCounts_norm_old = t(t(circCounts2)/sample_anno$totalReads)*10^6

##### SAVE DATA ######

circRNA_Finder_NMIBC = c("circTable", "circCounts", "circCounts_norm", "circCounts_fw", "circCounts_ratio", "circCounts_norm_old", "sample_anno")
save(list=circRNA_Finder_NMIBC, file=paste0(path, "/CE2_CIRIquant_NMIBC.Rdata"))

# sample anno, circTable & circCounts in the 230113_CIRIquant_circexplorer2 folder.
# then use as input for intron_assignment.R also in CIRIquant_circexplorer2 folder.



