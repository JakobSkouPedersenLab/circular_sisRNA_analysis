########################################################################################################################
############################### Make Rdata frames from reads counts in all introns hg38 ################################
########################################################################################################################
# 27/07/20
# Run this script to generate expression matrices from for all introns.

# args[4]

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

############### To run on cluster ###############
# Run the script:
#     srun --mem=64g -c 1 --time=5:0:0 --pty bash     # Get a node for 5 hours
#     conda activate run-R      # Activate conda environment where I have installed R
#     Rscript --vanilla /path/intron_counts_output.R /path/htseq/output/


########################################################################################################################
######################################### Set path and load samples ####################################################
########################################################################################################################
# The output files are named [sampleName].htseq.introns.tsv.gz

# At the cluster I will provide the path as an input in the following format: 
path = args[1]
dirs <- list.files(path = path, recursive = TRUE, include.dirs = TRUE, pattern="htseq.introns.tsv.gz") 
length(dirs)
head(dirs)

samples <- sub("/.*","",dirs)
samples
length(samples)
length(unique(samples))


########################################################################################################################
######################################### Generate data frame with reads counts  #######################################
########################################################################################################################
# Use first sample as reference
df1 <- read.table(gzfile(paste0(path,dirs[1])), as.is = TRUE, fill=TRUE, header = FALSE)

head(df1)
dim(df1)
tail(df1)
colnames(df1) = c("name", samples[1])

# Extract counts from all samples
count <- 1 
# length(dirs)
for (i in 2:length(dirs)){
  count <- count+1
  cat(count," ",samples[i],dim(df1),"\n")
  df2 <- read.table(gzfile(paste0(path,dirs[i])), as.is = TRUE, fill=TRUE, header=FALSE)
  colnames(df2) = c("position", samples[i])
  df1 <- merge(df1,df2,by.x=c("name"),by.y=c("position"),all=TRUE) 
}

introns_counts = df1
head(introns_counts)
dim(introns_counts)

# set rownames as names
rownames(introns_counts) = introns_counts$name

# generate matrix
intron_counts.mat = as.matrix(introns_counts[, 2:ncol(introns_counts)])
head(intron_counts.mat)


########################################################################################################################
################################################# Annotate introns #####################################################
########################################################################################################################

# load gtf file with intron positions
#BiocManager::install(c("rtracklayer"))
library(rtracklayer)

hg38_introns <- rtracklayer::import("/path/hg38v33_ctat_introns_noOverlap_withGeneId.gtf") 
head(hg38_introns)
hg38_introns = as.data.frame(hg38_introns)
dim(hg38_introns)
head(hg38_introns)

# check
hg38_introns[grep("HIPK3", hg38_introns$gene_name), ]

hg38_introns[grep("HNRNPK", hg38_introns$gene_name), ]
hg38_introns[grep("83976994", hg38_introns$end), ] # grep for HNRNPK intron that correspond to ciHNRNPK; chr9 83975722 83976994
hg38_introns[grep("48598958", hg38_introns$start), ] # WDR13; chrX 48598958 48599352

# order intron_counts.mat according to hg38_introns
intron_counts.mat = intron_counts.mat[match(hg38_introns$intron_id, rownames(intron_counts.mat)), ]
head(intron_counts.mat)
identical(rownames(intron_counts.mat), hg38_introns$intron_id)

# Rename type-column to introns
table(hg38_introns$type)
hg38_introns$type = "intron"

################ expression information ################
hg38_introns$sumReads = apply(intron_counts.mat, 1, sum)
head(hg38_introns)
hg38_introns[which.max(hg38_introns$sumReads), ]
quantile(hg38_introns$sumReads)

hg38_introns$mean = apply(intron_counts.mat, 1, mean)
hg38_introns$max = apply(intron_counts.mat, 1, max)
hg38_introns$sumExp = apply(intron_counts.mat, 1, function(x) sum(x > 0))
hg38_introns$median = apply(intron_counts.mat, 1, median)
head(hg38_introns)

hg38_introns[grep("HNRNPK", hg38_introns$gene_name), ]
hg38_introns[grep("83976994", hg38_introns$end), ]
hg38_introns[grep("WDR13", hg38_introns$gene_name), ]

nrow(hg38_introns)
nrow(intron_counts.mat)


########################################################################################################################
######################################### Save data frames as R data ###################################################
########################################################################################################################
Intron_counts_list = c("Intron_counts_list", "intron_counts.mat", "hg38_introns")
save(list=Intron_counts_list, file=paste0(path, "Intron_counts.Rdata"))
