#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
###
#
# For each region in the BED files, fetches the read count in that regions from the BAM file.
#
###

##bed files

res_ctrl <- read.csv(file = "/mnt/data/DeSeq_H4/ctrl_h4ac_peaks.csv")
res_rest <- read.csv(file = "/mnt/data/DeSeq_H4/rest_h4ac_peaks.csv")

##filter for pval cut off, distance to TSS and peak size

ctrl <- res_ctrl[which(res_ctrl[,"neg10log10pVal"] >= 100),]
rest <- res_rest[which(res_rest[,"neg10log10pVal"] >= 100),]

ctrl <- ctrl[which(ctrl[,"Peak_width"]<=10000),]
rest <- rest[which(rest[,"Peak_width"]<=10000),]

ctrl <- ctrl[which(abs(ctrl[,"distancetoFeature"]) <= 10000),]
rest <- rest[which(abs(rest[,"distancetoFeature"]) <= 10000),]

ctrl <- ctrl[,3:6]
rest <- rest[,3:6]

for(i in 1:length(ctrl[,1])){
	ctrl[i,5] <- ctrl[i,"Peak_width"] / 2
	ctrl[i,6] <- (ctrl[i,"Peak_start"] + ctrl[i,5]) - 250
	ctrl[i,7] <- (ctrl[i,"Peak_end"] - ctrl[i,5]) + 250
}

for(i in 1:length(rest[,1])){
	rest[i,5] <- rest[i,"Peak_width"] / 2
	rest[i,6] <- (rest[i,"Peak_start"] + rest[i,5]) - 250
	rest[i,7] <- (rest[i,"Peak_end"] - rest[i,5]) + 250
}

ctrl[,1] <- paste("chr",ctrl[,1],sep = "")
rest[,1] <- paste("chr",rest[,1],sep = "")

##remove 1 from start as bed starts at 0

ctrl[,6] <- ctrl[,6] - 1
rest[,6] <- rest[,6] - 1

##reorder so that it works

ctrl <- ctrl[,c(1,6,7)]
rest <- rest[,c(1,6,7)]

colnames(ctrl) <- c("Chromosome", "Peak_start","Peak_end")
colnames(rest) <- c("Chromosome", "Peak_start","Peak_end")

write.table(ctrl, file = "/mnt/data/DeSeq_H4/ctrl_H4peaks_for_deseq_middlepeak.bed",row.names = F, col.names = F, quote = F)
write.table(rest, file = "/mnt/data/DeSeq_H4/rest_H4peaks_for_deseq_middlepeak.bed",row.names = F, col.names = F, quote = F)

# rest_h4ac     = 'FloResCre__C18__H4ac_CME117_s_2_export_sorted_nodups.bed'
# ctrl_h4ac     = 'FloRes_H4acs_CME118_3_export_sorted_nodups.bed'


# you will need:
# The BAM files for each of the samples
# A BED file containing the peak regions for each of the samples


#!/usr/local/bin/Rscript peak_read_counts.R /path/to/outdir N_threads /path/to/bam1 /path/to/bed1 /path/to/bam2 /path/to/bed2 ...

options(stringsAsFactors = FALSE);

args <- commandArgs(trailingOnly=TRUE)

#for testing
#args<-c("/mnt/data/",4, "/mnt/astro/BAM/IP.bam" ,"/mnt/astro/Macs/NA_peaks.bed", "/mnt/esc/chip_export_sorted_nodups.bam","/mnt/esc/macs_300_1.0e-05/EscChIPseqREST_peaks.bed")
#args <- c("/mnt/data/",4,"/mnt/TC/CMN066_s_8_export_sorted_nodups.bam", "/mnt/TC/macs_300_1.0e-05/Mash1_TC_peaks.bed", "/mnt/SC/chip_export_sorted_nodups.bam", "/mnt/SC/macs_300_1.0e-05/Mash1_SC_peaks.bed")

args<-c("/mnt/data/",4, "/mnt/data/FloRes_H4acs_CME118_3_export_sorted_nodups.bam" ,"/mnt/data/DeSeq_H4/ctrl_H4peaks_for_deseq_middlepeak.bed", "/mnt/data/FloResCre__C18__H4ac_CME117_s_2_export_sorted_nodups.bam","/mnt/data/DeSeq_H4/rest_H4peaks_for_deseq_middlepeak.bed")

resdir <- args[1]
threads <- args[2]
args <- args[-1:-2]

inds <- 1:length(args)
bam.files <- args[which(inds%%2!=0)]
bed.files <- args[which(inds%%2==0)]


library(IRanges)
library(ShortRead)
library(snow)
library(baySeq)
library(DESeq)
library(rtracklayer)
library(Rsamtools)

# For ChIPseq data we aren't dealing with that many locations, so probably
# we're good using just the number of cores on the AWS machine. 
#if(is.null(threads) || threads==1){
#  cl <- NULL
#}else{
#  cl <- makeCluster(threads,"SOCK")
#}


# Read in the BED files as RangedData.
# Retrieve the read data for those regions
# Build the counts table
# This seems to take about an hour for 4K regions.

beds <- list()
counts <- NULL
seglens <- NULL
for(i in 1:length(bed.files)){
  beds[[i]] <- import(bed.files[i])

  #get the count data for these ranges from each bam file
  bam.counts <- list()
  for(j in 1:length(bam.files)){
    what <- c("qname") 
    param <- ScanBamParam(which=beds[[i]], what=what)
    bam <- scanBam(bam.files[[j]], param=param)
    bam.counts[[j]] <- sapply(bam, function(x){length(x$qname)})
  }

  beds[[i]] <- as.data.frame(beds[[i]])
  nms <- paste(beds[[i]][,"space"], paste(beds[[i]][,"start"], beds[[i]][,"end"], sep="-"), sep=":")
  these <- do.call(cbind, bam.counts)
  col.nms <- sub("-","", gsub("/","-", bam.files))
  colnames(these) <- col.nms
  these <- data.frame(sample=bed.files[i], these)
  if(is.null(counts)){
    counts <- these
  }else{
    counts <- rbind(counts, these)
  }
  if(is.null(seglens)){
    seglens <- beds[[i]][,"end"]-beds[[i]][,"start"]+1
  }else{
    seglens <- c(seglens,beds[[i]][,"end"]-beds[[i]][,"start"]+1)
  }
}

save(counts, file=paste(resdir,"DeSeq_H4/counts_middlepeak.RData", sep="/"))



# We'll need the library sizes later, so get them from the summary files
# Note that these should be generated using the samtools flagstat program
summaries <- gsub('.bam', '.summary', bam.files )
libsizes <- sapply(summaries, function(x){
  l <- readLines(x,1)
  l <- sub("\\s+.*","",l, perl=T)
})


save(libsizes, file=paste(resdir,"DeSeq/libsizes_middlepeak.RData", sep="/"))
