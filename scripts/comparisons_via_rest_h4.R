
options(stringsAsFactors = FALSE)

#pull in REST binding sites from Larry's REST ESC data

restbind <- read.csv(file = "results/REST_binding_sites.csv")

restbind <- restbind[which(restbind[,"neg10log10pVal"] >=100),]
restbind <- restbind[which(abs(restbind[,"distancetoFeature"]) <= 10000),]

#get DeSeq data for H4ac peak changes

peaks <- read.csv(file = "results/DeSeq_H4/peak_compare_middlepeak.csv")

# need data about all H3K9ac peaks in ctrl and rest ko

ctrl <- read.csv(file = "results/DeSeq_H4/ctrl_h4ac_peaks.csv")
rest <- read.csv(file = "results/DeSeq_H4/rest_h4ac_peaks.csv")

##filter for pval cut off, distance to TSS and peak size

ctrl <- ctrl[which(ctrl[,"neg10log10pVal"] >= 100),]
rest <- rest[which(rest[,"neg10log10pVal"] >= 100),]


## remove peaks that are massive - in RESTKO are 22.8 and 33.4kb, rest below 10kb

ctrl <- ctrl[which(ctrl[,"Peak_width"]<=10000),]
rest <- rest[which(rest[,"Peak_width"]<=10000),]

ctrl <- ctrl[which(abs(ctrl[,"distancetoFeature"]) <= 10000),]
rest <- rest[which(abs(rest[,"distancetoFeature"]) <= 10000),]

## find unique and shared peaks, will be useful later...

shared <- intersect(ctrl[,"EnsemblID"],rest[,"EnsemblID"])

unique_ctrl <- ctrl[which(!(ctrl[,"EnsemblID"] %in% shared)),"EnsemblID"]
unique_rest <- rest[which(!(rest[,"EnsemblID"] %in% shared)),"EnsemblID"]

## draw box and whisker plot of K9Ac size/width in each sample

library(vioplot)

names <- c("ESC WT", "ESC REST-KO")

boxplot(ctrl[which(ctrl[,"EnsemblID"] %in% unique_ctrl), "FoldEnrichment"], 
        rest[which(rest[,"EnsemblID"] %in% unique_rest), "FoldEnrichment"],
	names = names,
	ylab = "H3K9ac Fold Enrichment",
	outline = FALSE
	)

## try with all peaks

boxplot(ctrl[,"FoldEnrichment"], 
        rest[,"FoldEnrichment"],
	names = names,
	ylab = "H3K9ac Fold Enrichment",
	outline = FALSE
	)

## try with those that have REST bound in control

ctrl_rest <- merge(ctrl, restbind, by.x = "EnsemblID", by.y = "EnsemblID", suffixes = c("_K9ac","_REST"))
rest_rest <- merge(rest, restbind, by.x = "EnsemblID", by.y = "EnsemblID", suffixes = c("_K9ac","_REST"))

boxplot(ctrl_rest[,"FoldEnrichment_K9ac"], 
        rest_rest[,"FoldEnrichment_K9ac"],
	names = names,
	ylab = "H3K9ac Fold Enrichment",
	outline = FALSE
	)

## try those that are shared between ctrl and restko

boxplot(ctrl[which(ctrl[,"EnsemblID"] %in% shared), "FoldEnrichment"], 
        rest[which(rest[,"EnsemblID"] %in% shared), "FoldEnrichment"],
	names = names,
	ylab = "H3K9ac Fold Enrichment",
	outline = FALSE
	)


## concatenate ctrl and restko peak info together, then pull out genes with REST peaks for plotting...

res_all <- rbind(ctrl,rest)

res_all <- res_all[,2:21]

###because some peaks will be in both ctrl and restko, need to remove duplicates - order by size of peak and keep biggest

res_all.o <- res_all[order(res_all[,"FDR"], decreasing = FALSE),]
res_all.od <- res_all.o[!duplicated(res_all.o[,"Peak"]),]

# then merge with REST binding sites

#tidy up before merging

restbind <- restbind[,c(2,3,4,5,6,7,8,9,10,11,16)]

h4_rest <- merge(restbind, res_all.od, by.x = "EnsemblID", by.y = "EnsemblID", suffixes = c("_REST","_H4"))

## and then merge with deseq data

#but first give sensible colnames

peaks <- peaks[,1:8]

colnames(peaks) <- c("Peak","baseMean","Ctrl_count","Rest_count", "FC_DeSeq", "logFC_DeSeq", "PVal_DeSeq", "PAdj_DeSeq")

# also names dont match - as now taking 500bp region in middle of peak - rehash this from k9_rest

for(i in 1:length(h4_rest[,1])){
	h4_rest[i,31] <- h4_rest[i,"Peak_width_H4"] / 2
	h4_rest[i,32] <- (h4_rest[i,"Peak_start_H4"] + h4_rest[i,31]) - 250
	h4_rest[i,33] <- (h4_rest[i,"Peak_end_H4"] - h4_rest[i,31]) + 250
	h4_rest[i,34] <- paste("chr",paste(paste(h4_rest[i,"Chromosome_H4"],h4_rest[i,32],sep = ":"),h4_rest[i,33],sep = "-"), sep ="")
}

h4_rest <- h4_rest[,c(1:30,34)]
colnames(h4_rest)[31] <- "Peak_window_DeSeq"


h4_rest_deseq <- merge(h4_rest, peaks, by.x = "Peak_window_DeSeq", by.y = "Peak")

# save

write.csv(h4_rest_deseq, file = "results/H4ac_vs_REST_binding_vs_DeSeqH4ac_middlepeak.csv")

########## plot change in H4ac vs REST and try to draw correlation

# remove REST (line 463) from dataset and any other outliers

h4vsRest <- h4_rest_deseq[which(!(h4_rest_deseq[,"Symbol"] == "Rest")),]

plot(h4vsRest[,"FoldEnrichment_REST"], h4vsRest[,"logFC_DeSeq"])

#cor(h4vsRest[,"FoldEnrichment_REST"], h4vsRest[,"logFC_DeSeq"])

fit <- lm(h4vsRest[,"logFC_DeSeq"] ~ h4vsRest[,"FoldEnrichment_REST"])

abline(fit)

## try robust regression - ignores outliers better

library(MASS)

fitrob <- rlm(h4vsRest[,"logFC_DeSeq"] ~ h4vsRest[,"FoldEnrichment_REST"] + 0)

plot(h4vsRest[,"FoldEnrichment_REST"], h4vsRest[,"logFC_DeSeq"])

abline(fitrob)

## try making this look pretty - smoothscatter like in Cao2010...

library(geneplotter)

postscript(file = "results/DeSeq_H4/REST_vs_H4ac_changes.ps", horizontal = FALSE)
smoothScatter(h4vsRest[,"FoldEnrichment_REST"], h4vsRest[,"logFC_DeSeq"], main = "Correlation between REST binding and change in H4ac", xlab = "Amount of REST bound", ylab = "Change in H4ac")
abline(fitrob, lty = 2, lwd = 3, col = "red")
dev.off()

## then pull in expression data

exp <- read.csv(file = "results/DeSeq_K9ac/day0.csv")

# remove duplicates

exp.o <- exp[order(exp[,"adj.P.Val"], decreasing = FALSE),]
exp.od <- exp.o[!duplicated(exp.o[,"EnsemblID"]),]

exp.odt <- exp.od[,c(1,2,10,13,14)]

h4_rest_deseq_exp <- merge(h4_rest_deseq, exp.odt, by.x = "EnsemblID", by.y = "EnsemblID")

##should also merge with all K9ac vs exp - remove REST genes?
## merge res_all.od with deseq and expression data

for(i in 1:length(res_all.od[,1])){
	res_all.od[i,21] <- res_all.od[i,"Peak_width"] / 2
	res_all.od[i,22] <- (res_all.od[i,"Peak_start"] + res_all.od[i,21]) - 250
	res_all.od[i,23] <- (res_all.od[i,"Peak_end"] - res_all.od[i,21]) + 250
	res_all.od[i,24] <- paste("chr",paste(paste(res_all.od[i,"Chromosome"],res_all.od[i,22],sep = ":"),res_all.od[i,23],sep = "-"), sep ="")
}

res <- res_all.od[,c(1:20,24)]

# remove REST from data

res <- res[which(!(res[,"Symbol"] == "Rest")),]

colnames(res)[21] <- "Peak_window_DeSeq"

# merge with deseq

h4_deseq <- merge(res, peaks, by.x = "Peak_window_DeSeq", by.y = "Peak")

# merge with expression data

h4_deseq_exp <- merge(h4_deseq, exp.odt, by.x = "EnsemblID", by.y = "EnsemblID")

## now can plot change in K9ac vs gene expression

plot(h4_deseq_exp[,"logFC"],h4_deseq_exp[,"logFC_DeSeq"])

## looks nicer in hexbin

library(hexbin)

bin <- hexbin(h4_deseq_exp[,"logFC"],h4_deseq_exp[,"logFC_DeSeq"],xbins = 50)

postscript(file = "results/DeSeq_H4/change_in_h4ac_versus_gene_expression.ps", horizontal = FALSE)
plot(bin, main = "Change in H4ac versus Gene Expression", xlab = "Gene Expression (log Fold Change)", ylab = "Change in H4ac")
dev.off()


## remove REST and Syt4 and plot change in K9Ac and change in gene expression between WT and REST-KO

h4_rest_deseq_exp <- h4_rest_deseq_exp[which(!(h4_rest_deseq_exp[,"Symbol"] == "Rest")

postscript(file = "results/DeSeq_H4/exp_vs_h4ac_direct_REST_targets.ps", horizontal = FALSE)
plot(h4_rest_deseq_exp[c(1:172,174:408),"logFC"], h4_rest_deseq_exp[c(1:172,174:408),"logFC_DeSeq"], ylab = "H4ac Changes", xlab = "Gene Expression Changes", main = "H4ac vs Gene Expression - direct REST targets", pch = 20)
dev.off()


###plot direct targets and see expression changes

sym <- read.csv(file = "results/DeSeq_K9ac/direct_REST_targets_forGO.csv")

sym <- sym[,2]

sym_id <- as.character(k9_rest_deseq_exp[which(k9_rest_deseq_exp[,"Symbol"] %in% sym),"ID"])


## get timecourse expression data

library(beadarray)
library(limma)
library(biomaRt)
library(illuminaMousev2BeadID.db)
library(illuminaMousev2.db)
library(gplots)
library(Mfuzz)
library(hopach)

BSData <- get(load("/space/matt/Chiara_d4_c18_timecourse_expression/results/BSData.quantile.RData"))
E <- exprs(BSData)
E <- E[,c(1:4,24:27,28:31,5:8,9:11,16:19,20:23,12:15)]

##calculate averages

D4.0 <- c(1,2,3,4)
D4.1 <- c(5,6,7,8)
D4.2 <- c(9,10,11,12)
D4.4 <- c(13,14,15,16)
C18.0 <- c(17,18,19)
C18.1 <- c(20,21,22,23)
C18.2 <- c(24,25,26,27)
C18.4 <- c(28,29,30,31)

aves <- apply(E, 1, function(x){
        D4.0av <- sum(x[D4.0])/length(D4.0)
        D4.1av <- sum(x[D4.1])/length(D4.1)
        D4.2av <- sum(x[D4.2])/length(D4.2)
        D4.4av <- sum(x[D4.4])/length(D4.4)
        C18.0av <- sum(x[C18.0])/length(C18.0)
        C18.1av <- sum(x[C18.1])/length(C18.1)
        C18.2av <- sum(x[C18.2])/length(C18.2)
        C18.4av <- sum(x[C18.4])/length(C18.4) 
            return(c(D4.0av,D4.1av,D4.2av,D4.4av,C18.0av,C18.1av,C18.2av,C18.4av))
        }
)        

aves <- t(aves)

colnames(aves) <- c("D4_0","D4_1","D4_2","D4_4","C18_0","C18_1","C18_2","C18_4")

sym_aves <- aves[sym_id,]

names <- c("Day 0", "Day 1", "Day 2", "Day 4")

plot(sym_aves[1,1:4], type = "l",ylab="Expression level (log2)", xlab = "Timecourse", xaxt = "n", ylim=c(6.5,14), main = "Generation of Neurons",col="blue")
axis(1, at=1:4, labels=names)
lines(sym_aves[1,5:8], col = "red")




## maybe pull out probe ids from day0 to 4 comparison and take top changing genes

expdif <- read.csv(file = "/space/matt/Chiara_d4_c18_timecourse_expression/results/day4vs0_D4.csv")

expdif.o <- expdif[order(expdif[,"adj.P.Val"], decreasing = FALSE),]
expdif.od <- expdif.o[!duplicated(expdif.o[,"EnsemblID"]),]

expsym <- expdif.od[which(expdif.od[,"symbol"] %in% sym),c("ID","logFC")]

expsymtop <- as.character(expsym[order(abs(expsym[,"logFC"]), decreasing = TRUE),"ID"])

sym_aves <- aves[expsymtop,]

postscript(file="results/direct_targets_timecourse.ps", horizontal=FALSE)
attach(mtcars)
par(mfrow=c(2,3))

plot(sym_aves[1,1:4], type = "l",ylab="Expression level (log2)", xlab = "Timecourse", xaxt = "n", ylim=c(6.5,14), main = "Igsf21",col="blue")
axis(1, at=1:4, labels=names)
lines(sym_aves[1,5:8], col = "red")

plot(sym_aves[2,1:4], type = "l",ylab="Expression level (log2)", xlab = "Timecourse", xaxt = "n", ylim=c(6.5,14), main = "Tcea3",col="blue")
axis(1, at=1:4, labels=names)
lines(sym_aves[2,5:8], col = "red")

plot(sym_aves[3,1:4], type = "l",ylab="Expression level (log2)", xlab = "Timecourse", xaxt = "n", ylim=c(6.5,14), main = "Trim36",col="blue")
axis(1, at=1:4, labels=names)
lines(sym_aves[3,5:8], col = "red")

plot(sym_aves[4,1:4], type = "l",ylab="Expression level (log2)", xlab = "Timecourse", xaxt = "n", ylim=c(6.5,14), main = "Dynll2",col="blue")
axis(1, at=1:4, labels=names)
lines(sym_aves[4,5:8], col = "red")

plot(sym_aves[5,1:4], type = "l",ylab="Expression level (log2)", xlab = "Timecourse", xaxt = "n", ylim=c(6.5,14), main = "Cyp51",col="blue")
axis(1, at=1:4, labels=names)
lines(sym_aves[5,5:8], col = "red")

plot(sym_aves[6,1:4], type = "l",ylab="Expression level (log2)", xlab = "Timecourse", xaxt = "n", ylim=c(6.5,14), main = "Tcea3",col="blue")
axis(1, at=1:4, labels=names)
lines(sym_aves[6,5:8], col = "red")

dev.off()

## try to plot difference in size between ctrl and restko - take unique K9ac size compared to REST height
## or just simply the size of the top say 100 K9ac peaks in each sample

# need data about all H3K9ac peaks in ctrl and rest ko

ctrl <- read.csv(file = "results/DeSeq_K9ac/ctrl_h3k9ac_peaks.csv")
rest <- read.csv(file = "results/DeSeq_K9ac/rest_h3k9ac_peaks.csv")

ctrl <- ctrl[which(ctrl[,"neg10log10pVal"] >= 100),]
rest <- rest[which(rest[,"neg10log10pVal"] >= 100),]

ctrl <- ctrl[which(ctrl[,"Peak_width"]<=6000),]
rest <- rest[which(rest[,"Peak_width"]<=6000),]

ctrl <- ctrl[which(abs(ctrl[,"distancetoFeature"]) <= 3000),]
rest <- rest[which(abs(rest[,"distancetoFeature"]) <= 3000),]



shared <- intersect(ctrl[,"EnsemblID"], rest[,"EnsemblID"])

uni_ctrl <- ctrl[which(!(ctrl[,"EnsemblID"] %in% shared)),]
uni_rest <- rest[which(!(rest[,"EnsemblID"] %in% shared)),]

## order by size and draw 100 biggest
uni_ctrl <- uni_ctrl[order(uni_ctrl[,"neg10log10pVal"], decreasing = TRUE),]
uni_rest <- uni_rest[order(uni_rest[,"neg10log10pVal"], decreasing = TRUE),]

plot(uni_ctrl[1:200,"Peak_width"], uni_rest[1:200,"Peak_width"])

shared_ctrl <- ctrl[which(ctrl[,"EnsemblID"] %in% shared),]
shared_rest <- rest[which(rest[,"EnsemblID"] %in% shared),]

shared_res <- merge(shared_ctrl,shared_rest, by.x = "EnsemblID", by.y = "EnsemblID", suffixes = c("_ctrl","_rest"))

plot(shared_res[,"Peak_width_ctrl"], shared_res[,"Peak_width_rest"])

















