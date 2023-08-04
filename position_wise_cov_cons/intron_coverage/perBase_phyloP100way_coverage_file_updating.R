###################################################################################
####### Match coverage and conservation files to sisRNA and intron anno data ######
###################################################################################

# Necessary to reduce file size, match naming by adding Updated_ as prefix for updated files.
# Also compute summary statistics here, because data is too large to handle locally.

# When running bash script: conda activate r_env

##### Load data and packages #####
library(tidyverse)
library(rtracklayer)
library(GenomicFeatures) 

##### Load general data #####

# User mount path:
mount_path <- "/path/to/mount"

# Conservation
conservation <- read_tsv(paste(mount_path, "/path/phyloP100Way/hg38.phyloP100way.intron.bed", sep = ""), col_names = F)
mean_conservation <- read_tsv(paste(mount_path, "/path/phyloP100Way/hg38.phyloP100way.intron.tab", sep = ""), col_names = F)

# NMIBC Coverage
coverage <- read_tsv(paste(mount_path, "/path/mean_intron_coverage_sorted.bed", sep = ""), col_names = F) 
mean_coverage <- read_tsv(paste(mount_path, "/path/mean_intron_coverage_sorted_norm_perSample.tab", sep = ""), col_names = F)

# sisRNA annotation 
load(paste(mount_path, "/path/200805_circexplorer2/Updated_merged_ciRNAs_NMIBC.Rdata", sep = ""))

total_read_depth <- 1000000 # Add correct value here
nr_samples <- 100 # Add correct value here


######################################################
########## Align branch points from sisRNAs ##########
######################################################

# get non-exon-non-intron overlapping intron positions
merged_ciRNA_anno$tmp1 <- unlist(lapply(merged_ciRNA_anno$intron_id, function(x) unlist(str_split(x, pattern = ":"))[2]))
merged_ciRNA_anno$intron_start2 <- unlist(lapply(merged_ciRNA_anno$tmp1, function(x) unlist(str_split(x, pattern = "-"))[1]))
merged_ciRNA_anno$intron_start2 <- as.numeric(merged_ciRNA_anno$intron_start2)
merged_ciRNA_anno$tmp2 <- unlist(lapply(merged_ciRNA_anno$tmp1, function(x) unlist(str_split(x, pattern = "-"))[2]))
merged_ciRNA_anno$intron_end2 <- unlist(lapply(merged_ciRNA_anno$tmp2, function(x) unlist(str_split(x, pattern = "_"))[1]))
merged_ciRNA_anno$intron_end2 <- as.numeric(merged_ciRNA_anno$intron_end2)

merged_idx_coding <- merged_ciRNA_anno$gene_type == "protein_coding"
merged_idx_circle <- merged_ciRNA_anno$name == merged_ciRNA_anno$flankIntron

merged_ciRNA_anno$coding <- ifelse(merged_idx_coding, "protein_coding", "non_coding")
merged_ciRNA_anno$circle <- ifelse(merged_idx_circle, "circle", "lariat")

merged_ciRNA_anno$BP_3SS_dist <- ifelse(merged_ciRNA_anno$strand == "-", 
                                        merged_ciRNA_anno$start-merged_ciRNA_anno$intron_start2, 
                                        merged_ciRNA_anno$intron_end2-merged_ciRNA_anno$end)

# keep PC, ones with length >= 0. for circle set NA before BP position and only use values after (avoid exonic positions).
quantile(na.omit(merged_ciRNA_anno$BP_3SS_dist))
length(which(na.omit(merged_ciRNA_anno$BP_3SS_dist) < 0)) # 300
length(which(na.omit(merged_ciRNA_anno$BP_3SS_dist) < 25))

merged_ciRNA_anno <- merged_ciRNA_anno %>% filter(coding == "protein_coding", BP_3SS_dist >= 0)


######################################################
############# Conservation data update ###############
######################################################

conservation <- as.data.frame(conservation)
mean_conservation <- as.data.frame(mean_conservation)

# Add colnames, merge and remove identical columns
colnames(conservation) <- c("chr", "start", "end", "name", "width", "strand", "nr_bases", "cons_per_base")
identical(conservation$width, conservation$nr_bases) # T
conservation <- conservation[, -7]

colnames(mean_conservation) <- c("name", "width_U", "nr_bases_cons_U", "sum_U", "mean0_U", "mean_U")
identical(mean_conservation$width, mean_conservation$nr_bases_cons) # F
identical(mean_conservation$width, conservation$width) # T
mean_conservation <- mean_conservation[, -2]

identical(mean_conservation$name, conservation$name) # T
conservation <- cbind(conservation, mean_conservation[, c(2:5)])


##### Make summary statistic for conservation location within region #####

# cons_per_base column is too large, compute summary stat and then remove it.

# 50 bp 3' SS and 5' SS conservation/coverage function
perBase_stat <- function(stat_per_base, cov, type, strand){
  vector <- as.numeric(unlist(strsplit(stat_per_base, split = ",")))
  len <- 100
  idx <- 50
  
  if(cov == "cov"){
    vector <- (vector*1000000/total_read_depth)/nr_samples # Normalize from raw counts to CPM per sample
  }
  
  if(type == "5SS" & length(vector) >= len){
      vector <- vector[1:(length(vector)-10)] # trim 10 nt at 3' end (remove potential 3' SS signal)
      col_out <- vector[1:idx]
    return(col_out) # return summary stat list
  } 
  
  if(type == "3SS" & length(vector) >= len){ 
      vector <- vector[10:length(vector)] # trim 10 nt at 5' end (remove potential 5' SS signal)
      col_out <- vector[length(vector):(length(vector)-(idx-1))]
    return(col_out) # return summary stat list
  }
  
  return(rep(NA, idx))
} 

# Add per base statistics to DF
# 3' SS
stat_out <- lapply(1:length(conservation$cons_per_base), function(x) perBase_stat(conservation$cons_per_base[x], "cons", "3SS", conservation$strand[x]))
df_out <- NA
df_out <- as.data.frame(do.call(rbind, stat_out))
rownames(df_out) <- conservation$name
colnames(df_out) <- 1:50
df_out <- df_out[rowSums(is.na(df_out)) != ncol(df_out), ]
# 5' SS
stat_out <- lapply(1:length(conservation$cons_per_base), function(x) perBase_stat(conservation$cons_per_base[x], "cons", "5SS", conservation$strand[x]))
df_out_5SS <- NA
df_out_5SS <- as.data.frame(do.call(rbind, stat_out))
rownames(df_out_5SS) <- conservation$name
colnames(df_out_5SS) <- 1:50
df_out_5SS <- df_out_5SS[rowSums(is.na(df_out_5SS)) != ncol(df_out_5SS), ]

##### Save data #####
save(df_out, file=paste(mount_path, "/path/phyloP100way_intron_perBase_3SS.Rdata", sep = ""))
save(df_out_5SS, file=paste(mount_path, "/path/phyloP100way_intron_perBase_5SS.Rdata", sep = ""))


######################################################
############### Coverage data update #################
######################################################

coverage <- as.data.frame(coverage)
mean_coverage <- as.data.frame(mean_coverage)

# Add colnames, merge and remove identical columns
colnames(coverage) <- c("chr", "start", "end", "name", "width", "strand", "nr_bases", "cov_per_base")
identical(coverage$width, coverage$nr_bases) # T
coverage <- coverage[, -7]

colnames(mean_coverage) <- c("name", "width_U", "nr_bases_cov_U", "sum_U", "mean0_U", "mean_U")
identical(mean_coverage$width, mean_coverage$nr_bases_cov) # F
identical(mean_coverage$width, coverage$width) # T
mean_coverage <- mean_coverage[, -2]

identical(mean_coverage$name, coverage$name) # T
coverage <- cbind(coverage, mean_coverage[, c(2:5)])


##### Make summary statistic for coverage location within region #####

# Add per base statistics to DF
stat_out <- lapply(1:length(coverage$cov_per_base), function(x) perBase_stat(coverage$cov_per_base[x], "cov", "3SS", conservation$strand[x]))
df_out_cov <- NA
df_out_cov <- as.data.frame(do.call(rbind, stat_out))
rownames(df_out_cov) <- coverage$name
colnames(df_out_cov) <- 1:50
df_out_cov <- df_out_cov[rowSums(is.na(df_out_cov)) != ncol(df_out_cov), ]

stat_out <- lapply(1:length(coverage$cov_per_base), function(x) perBase_stat(coverage$cov_per_base[x], "cov", "5SS", conservation$strand[x]))
df_out_cov_5SS <- NA
df_out_cov_5SS <- as.data.frame(do.call(rbind, stat_out))
rownames(df_out_cov_5SS) <- coverage$name
colnames(df_out_cov_5SS) <- 1:50
df_out_cov_5SS <- df_out_cov_5SS[rowSums(is.na(df_out_cov_5SS)) != ncol(df_out_cov_5SS), ]


##### Save data #####
save(df_out_cov, file=paste(mount_path, "/path/coverage_intron_perBase_3SS_norm_perSample.Rdata", sep = ""))
save(df_out_cov_5SS, file=paste(mount_path, "/path/coverage_intron_perBase_5SS_norm_perSample.Rdata", sep = ""))


print("perBase 5SS and 3SS finished")


