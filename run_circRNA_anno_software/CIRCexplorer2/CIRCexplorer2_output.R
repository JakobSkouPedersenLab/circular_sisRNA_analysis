########################################################################################################################
######################################### Make Rdata frames from CIRCexplorer2 output ##################################
########################################################################################################################
# 21/07/20
# Run this script to generate circRNA- and ciRNA expression matrices from CIRCexplorer2 output.

# args[4]

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

############### To run on cluster ###############
# Run the script:
#     srun --mem=64g -c 1 --time=10:0:0 --pty bash       # Get a node for 5 hours
#     conda activate run-R      # Activate conda environment where I have installed R
#     Rscript --vanilla /path/CIRCexplorer2_output.R /path/to/CIRCexplorer2_output/


########################################################################################################################
######################################### Set path and load samples ####################################################
########################################################################################################################
# The CIRCexplorer2 output files are named circularRNA_known.txt 

# At the cluster I will provide the path as an input in the following format: 
path = args[1]
dirs <- list.files(path = path, recursive = TRUE, include.dirs = TRUE, pattern="circularRNA_known.txt") # change according to output name
dirs
length(dirs)

samples <- sub("/.*","",dirs)
samples
length(samples)
length(unique(samples))


########################################################################################################################
######################################### Generate data frame with circular counts  ####################################
########################################################################################################################
# Use first sample as reference
cc1 <- read.table(paste0(path,dirs[1]), as.is = TRUE, fill=TRUE, header = FALSE)

head(cc1) # information about columns are found here: https://circexplorer2.readthedocs.io/en/latest/modules/annotate/
dim(cc1)
colnames(cc1) = c("chr", "start", "end", "name", "score", "strand", "thickstart", "thickend", "itemRgb", "exonCount", "exonSizes", "exonOffsets", "readNumber", "circType", "geneName", "isoformName", "index", "flankIntron")
tail(cc1)
table(cc1$circType)
quantile(cc1$readNumber)

# extract certain columns
cc <- cc1[,c("chr", "start", "end", "name", "readNumber", "strand", "circType", "geneName", "isoformName", "index", "flankIntron", "exonCount", "exonSizes", "exonOffsets")] # circType: type of circular RNA, e.g. "ciRNA" or "circRNA", index: index of exon or intron, flankIntron: position of entire intron for intronic ciRNAs, exonCount: Number of exons, exonSizes: Size of transcript (end-start) for both intronic and exonic circles, exonOffset: Exon offsets (?)
head(cc)
colnames(cc)[5] <- samples[1] # change "junction_reads"-column name to the name of the samples
head(cc)
cc$name = paste(cc$chr, ":", cc$start, "-", cc$end, sep="")

# Extract counts from all samples
count <- 1
for (i in 2:length(dirs)){
  count <- count+1
  cat(count," ",samples[i],dim(cc),"\n")
  df2 <- read.table(paste0(path,dirs[i]), as.is = TRUE, fill=TRUE, header=FALSE)
  colnames(df2) = c("chr", "start", "end", "name", "score", "strand", "thickstart", "thickend", "itemRgb", "exonCount", "exonSizes", "exonOffsets", "readNumber", "circType", "geneName", "isoformName", "index", "flankIntron")
  df2 <- df2[,c("chr", "start", "end", "name", "readNumber", "strand", "circType", "geneName", "isoformName", "index", "flankIntron", "exonCount", "exonSizes", "exonOffsets")] 
  df2$name = paste(df2$chr, ":", df2$start, "-", df2$end, sep="")
  colnames(df2)[4] = "circSite"
  colnames(df2)[5] = samples[i]
  cc <- merge(cc,df2,by.x=c("chr", "start", "end", "name", "strand", "circType", "geneName", "isoformName", "index", "flankIntron", "exonCount", "exonSizes", "exonOffsets"),by.y=c("chr", "start", "end", "circSite", "strand", "circType", "geneName", "isoformName", "index", "flankIntron", "exonCount", "exonSizes", "exonOffsets"),all=TRUE) 
}

head(cc)
dim(cc)

# change NA to 0 
circTable = cc
circTable[is.na(circTable)] = 0
head(circTable)
dim(circTable)
length(unique(circTable$name))
#circTable[which(circTable$name == "chr20:49279208-49280484"), ]

# CIRCexplorer2 is 0-based. Add 1 to start coordinate
class(circTable$start)
circTable$start = circTable$start+1
circTable$name = paste(circTable$chr, ":", circTable$start, "-", circTable$end, sep="")
#rownames(circTable) = circTable$name

rownames(circTable) = make.names(circTable$name, unique=TRUE) # Set unique row.names. The issue arise if a ciRNA with same coordinates are present in multiple rows with different isoformName and flankIntron, e.g. chr20:49279208-49280484 in the ENCODE tissue data set

# Divide data frames into circRNAs and ciRNAs
ciRNAs_circTable = circTable[which(circTable$circType == "ciRNA"), ]
circRNAs_circTable = circTable[which(circTable$circType == "circRNA"), ]
dim(ciRNAs_circTable)
dim(circRNAs_circTable)

ciRNAs_circTable.mat = as.matrix(ciRNAs_circTable[, 14:ncol(ciRNAs_circTable)])
circRNAs_circTable.mat = as.matrix(circRNAs_circTable[, 14:ncol(circRNAs_circTable)])
head(ciRNAs_circTable)
head(ciRNAs_circTable.mat)

# Add 1 to start coordinate of flankIntron column as well
ciRNAs_circTable$flankIntron_chr = sapply(ciRNAs_circTable$flankIntron, function(x) strsplit(x, split = ":")[[1]][1]) # column with chr, split based on ":"

ciRNAs_circTable$flankIntron_start_plus1 = sapply(ciRNAs_circTable$flankIntron, function(x) strsplit(x, split = ":")[[1]][2]) # column with start +1, split based on ":"
ciRNAs_circTable$flankIntron_start_plus1 = sapply(ciRNAs_circTable$flankIntron_start_plus1, function(x) strsplit(x, split = "-")[[1]][1]) # column with start +1, split based on "-"
ciRNAs_circTable$flankIntron_start_plus1 = as.integer(ciRNAs_circTable$flankIntron_start_plus1)+1 # column with start +1, add 1 to start position

ciRNAs_circTable$flankIntron_end = sapply(ciRNAs_circTable$flankIntron, function(x) strsplit(x, split = "-")[[1]][2]) # column with end, split based on "-"

ciRNAs_circTable$flankIntron = paste(ciRNAs_circTable$flankIntron_chr, ":", ciRNAs_circTable$flankIntron_start_plus1, "-", ciRNAs_circTable$flankIntron_end, sep="") # flankIntron coordinates with start+1

ciRNAs_circTable$flankIntron_chr = NULL # remove columns again
ciRNAs_circTable$flankIntron_start_plus1 = NULL
ciRNAs_circTable$flankIntron_end = NULL

head(ciRNAs_circTable)
ciRNAs_circTable[grep("WDR13", ciRNAs_circTable$geneName), ]


########################################################################################################################
######################################### Generate circAnno data frame #################################################
########################################################################################################################
# circAnno contains information about the circRNAs

################ circAnno for ciRNAs ################
ciRNAs_circAnno = ciRNAs_circTable[1:13] 
head(ciRNAs_circAnno)
names(ciRNAs_circAnno)[names(ciRNAs_circAnno)=="exonSizes"] = "intronSizes" # change exonSizes to intronSizes for intronic ciRNAs
table(ciRNAs_circAnno$exonCount)
table(ciRNAs_circAnno$exonOffsets)

ciRNAs_circAnno$exonCount = NULL # I remove exonOffsets and exonCounts for intronic ciRNAs, since it doesn't make any sense here
ciRNAs_circAnno$exonOffsets = NULL

ciRNAs_circAnno$sumReads = apply(ciRNAs_circTable.mat, 1, sum)
head(ciRNAs_circAnno)
x = which.max(ciRNAs_circAnno$sumReads)
ciRNAs_circAnno[x, ]
min(ciRNAs_circAnno$sumReads)
max(ciRNAs_circAnno$sumReads)

ciRNAs_circAnno$mean = apply(ciRNAs_circTable.mat, 1, mean)
ciRNAs_circAnno$max = apply(ciRNAs_circTable.mat, 1, max)
ciRNAs_circAnno$sumExp = apply(ciRNAs_circTable.mat, 1, function(x) sum(x > 0))
ciRNAs_circAnno$median = apply(ciRNAs_circTable.mat, 1, median)

head(ciRNAs_circTable.mat)
head(ciRNAs_circAnno)
ciRNAs_circAnno[grep("WDR13", ciRNAs_circAnno$geneName), ]

################ circAnno for circRNAs ################
circRNAs_circAnno = circRNAs_circTable[1:13] 
head(circRNAs_circAnno)
table(circRNAs_circAnno$exonCount)

circRNAs_circAnno$sumReads = apply(circRNAs_circTable.mat, 1, sum)
head(circRNAs_circAnno)
x = which.max(circRNAs_circAnno$sumReads)
circRNAs_circAnno[x, ]
min(circRNAs_circAnno$sumReads)
max(circRNAs_circAnno$sumReads)

circRNAs_circAnno$mean = apply(circRNAs_circTable.mat, 1, mean)
circRNAs_circAnno$max = apply(circRNAs_circTable.mat, 1, max)
circRNAs_circAnno$sumExp = apply(circRNAs_circTable.mat, 1, function(x) sum(x > 0))
circRNAs_circAnno$median = apply(circRNAs_circTable.mat, 1, median)

head(circRNAs_circTable.mat)
head(circRNAs_circAnno)

which(circRNAs_circAnno$geneName == "CDR1")
circRNAs_circAnno[which(circRNAs_circAnno$geneName == "TULP4"), ]

# as last output when I run script, get the number of ciRNAs
nrow(ciRNAs_circAnno)
ciRNAs_circAnno[grep("WDR13", ciRNAs_circAnno$geneName), ]


########################################################################################################################
######################################### Save data frames as R data ###################################################
########################################################################################################################
CIRCexplorer2 = c("CIRCexplorer2", "ciRNAs_circTable", "ciRNAs_circTable.mat", "ciRNAs_circAnno", "circRNAs_circTable", "circRNAs_circTable.mat", "circRNAs_circAnno")
save(list=CIRCexplorer2, file=paste0(path, "CIRCexplorer2.Rdata"))
