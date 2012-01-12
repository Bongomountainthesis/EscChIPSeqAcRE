
options(stringsAsFactors = FALSE)

################
# Take H3K9ac changes from DeSeq and annotate back to nearest_peak_to_tss file
# Check this works
# Use ensemblID to link back to microarray dataset
# Compare amount of H3K9ac to expression
#
# Repeat for REST...
#
################


ctrl <- read.csv(file = "ctrl_h3k9ac_peaks.csv")
rest <- read.csv(file = "rest_h3k9ac_peaks.csv")

##filter for pval cut off, distance to TSS and peak size

ctrl <- ctrl[which(ctrl[,"neg10log10pVal"] >= 100),]
rest <- rest[which(rest[,"neg10log10pVal"] >= 100),]

ctrl <- ctrl[which(ctrl[,"Peak_width"]<=6000),]
rest <- rest[which(rest[,"Peak_width"]<=6000),]

ctrl <- ctrl[which(abs(ctrl[,"distancetoFeature"]) <= 3000),]
rest <- rest[which(abs(rest[,"distancetoFeature"]) <= 3000),]

ctrl <- ctrl[,2:21]
rest <- rest[,2:21]

## add in deseq data

deseq <- read.csv(file = "peak_compare.csv")

deseq <- deseq[,1:8]

##merge peaks from deseq to original annotation

res <- merge(rest, deseq, by.x = "Peak", by.y = "id", all.x = TRUE)

# tidy up abit

res <- res[,c(1,2,3,4,5,6,7,8,9,10,15,16,17,18,19,11,20,22,23,24,25,26,27)]

## pull in expression data

#############timecourse microarray data
exp1 <- read.csv(file = "day0.csv")

## remove duplicate probes

exp1.o <- exp1[order(exp1[,"adj.P.Val"], decreasing = FALSE),]
exp1.od <- exp1.o[!duplicated(exp1.o[,"EnsemblID"]),]

exp1.sig <- exp1.od[which(exp1.od[,"adj.P.Val"] >= 0.05),]

## merge with deseq based on ensemblid

peakexp1 <- merge(res, exp1.sig, by.x = "EnsemblID", by.y = "EnsemblID")

## tidy up

peakexp1 <- peakexp1[,c(1,2,16,17,3,4,5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,32,35,36)]

colnames(peakexp1)[17] <- "Cont_count"
colnames(peakexp1)[18] <- "Rest_count"
colnames(peakexp1)[12] <- "FC_peak"
colnames(peakexp1)[19] <- "FC_RESTKO"
colnames(peakexp1)[20] <- "logFC_RESTKO"
colnames(peakexp1)[23] <- "logFC_exp"

postscript(file = "exp_vs_h3k9ac_TCdata.ps", horizontal = FALSE)
plot(peakexp1[,"logFC_RESTKO"], peakexp1[,"logFC_exp"], pch = ".", xlab = "Change in H3K9ac", ylab = "Change in Gene Expression (TC data)", main = "Correlation between H3K9ac and gene expression")
dev.off()

#############seperate array data
exp2 <- read.csv(file = "limma_results.csv")

## remove duplicate probes

exp2.o <- exp2[order(exp2[,"adj.P.Val"], decreasing = FALSE),]
exp2.od <- exp2.o[!duplicated(exp2.o[,"EnsemblID"]),]

exp2.sig <- exp2.od[which(exp2.od[,"adj.P.Val"] >= 0.05),]

## merge with deseq based on ensemblid

peakexp2 <- merge(res, exp2.sig, by.x = "EnsemblID", by.y = "EnsemblID")

## tidy up

peakexp2 <- peakexp2[,c(1,2,16,17,3,4,5,6,7,8,9,10,11,12,13,14,18,19,20,21,22,23,32,35,36)]

colnames(peakexp2)[17] <- "Cont_count"
colnames(peakexp2)[18] <- "Rest_count"
colnames(peakexp2)[12] <- "FC_peak"
colnames(peakexp2)[19] <- "FC_RESTKO"
colnames(peakexp2)[20] <- "logFC_RESTKO"
colnames(peakexp2)[23] <- "logFC_exp"

postscript(file = "exp_vs_h3k9ac_otherdata.ps", horizontal = FALSE)
plot(peakexp2[,"logFC_RESTKO"], peakexp2[,"logFC_exp"], pch = ".",xlab = "Change in H3K9ac", ylab = "Change in Gene Expression (Other data)", main = "Correlation between H3K9ac and gene expression")
dev.off()

### hmm not much correlation - maybe try for genes that have a RE1 site?

### 'superimpose' the rest binding events....

restbind <- read.csv(file = "../REST_binding_sites.csv")

#tidy up...

restbind <- restbind[,c(2,3,4,5,6,7,8,9,10,11,15,16,17)]

# take big sites

restbind <- restbind[which(restbind[,"neg10log10pVal"] >= 100),]

# take close ones - 10kb apparently...

restbind <- restbind[which(restbind[,"distancetoFeature"] <= 10000),]

# merge with data

peakexprest <- merge(peakexp1, restbind, by.x = "EnsemblID", by.y = "EnsemblID", all.x = TRUE, suffixes = c("_K9Ac","_REST"))

## which have REST bound?

restbound <- peakexprest[which(!is.na(peakexprest[,"Peak_REST"])),]

postscript(file = "exp_vs_h3k9ac_timecourse_REST_bound_timecourse.ps", horizontal = FALSE)
plot(restbound[,"logFC_RESTKO"], restbound[,"logFC_exp"], pch = ".",xlab = "Change in H3K9ac", ylab = "Change in Gene Expression (Timecourse data)", main = "Correlation between H3K9ac and gene expression - RE1 sites")
dev.off()


### just how different are the microarray datasets?

# timecourse exp1
# other set exp2

set1 <- exp1[,c("EnsemblID", "logFC")]
set2 <- exp2[,c("EnsemblID","logFC")]

res <- merge(set1,set2, by.x = "EnsemblID", by.y = "EnsemblID")

#### bollocks to this - pretty sure the counts.RData is normalised for read counts?

for(i in 1:length(rest_counts_avg[,1])){
     rest_counts_avg[i,4] <- rest_counts_avg[i,3] / rest_counts_avg[i,2]
    }

##merge back into peakexp data

peakexp_raw <- merge(peakexp1, rest_counts_avg, by.x = "Peak", by.y = "row.names")













