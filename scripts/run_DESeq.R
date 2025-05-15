# library(DESeq2)
# library(dplyr)
# library(clusterProfiler)
# library(tidyverse)
# library(biomaRt)
# library(enrichplot)
# library(pheatmap)
# library(RColorBrewer)
# library(vsn)
# library(apeglm)
# options(connectionObserver = NULL)
# library(org.Hs.eg.db)
# library(AnnotationHub)
# library(EnhancedVolcano)
#
# readcounts <- read.table("featCounts_stoon_genes.txt",header = TRUE)
# ognames <- colnames(readcounts)
# colnames(readcounts)[7:21] <-c(
# "MSK107Li.cc6.ClbPpos_1",
# "MSK107Li.cc6.ClbPpos_2",
# "MSK107Li.cc6.ClbPpos_3",
# "MSK107Li.cc6.dClbP_1",
# "MSK107Li.cc6.dClbP_2",
# "MSK107Li.cc6.dClbP_3",
# "MSK107Li.Week0_1",
# "MSK107Li.Week0_2",
# "MSK107Li.Week0_3",
# "MSK107Li.Week0.ClbPpos_1",
# "MSK107Li.Week0.ClbPpos_2",
# "MSK107Li.Week0.ClbPpos_3",
# "MSK107Li.Week0.dClbP_1",
# "MSK107Li.Week0.dClbP_2",
# "MSK107Li.Week0.dClbP_3"
# )
#
# row.names(readcounts)<- readcounts$Geneid
# de <- readcounts[,c(7:21)]
# sample_info <-DataFrame(condition =gsub("_[0-9]+", "",
#               names(de)),row.names =names(de) )
# de <-DESeqDataSetFromMatrix(
#   countData =as.matrix(de),
#   colData = sample_info,
#   design =~condition)
# de <-estimateSizeFactors(de)
# dim(de)
# de$condition <- relevel(de$condition, ref = "MSK107Li.Week0")
# dds <- DESeq(de)
# res <- results(dds, contrast=c("condition","MSK107Li.Week0.ClbPpos","MSK107Li.Week0"))


#!/usr/local/bin/Rscript

## INPUT
library(DESeq2)
padj_cut <- 0.05 # P value cutoff per test

args <- commandArgs(trailingOnly=TRUE)
cts_path <- args[1]
cts <- as.matrix(read.csv(cts_path, sep='\t', row.names='gene'))

coldata_path <- args[2] #'/juno/work/shah/users/chois7/tickets/crcorganoid98/136PC2/rnaseq/results/featurecounts/_ensg_level/ClbPpos_vs_W0.coldata.tsv'
coldata <- read.csv(coldata_path, sep='\t', row.names=1)
coldata$condition <- factor(coldata$condition)
if (!all(rownames(coldata) == colnames(cts))) { 
  stop(sprintf('Not all rownames == colnames'))
}

dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=~condition)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "Case", "Control"))
resOrdered <- res[order(res$pvalue),]
resSig <- subset(resOrdered, padj < padj_cut)
foutName <- paste(strsplit(cts_path, ".table")[1], "deseq.csv", sep=".") # ControlID_vs_CaseID.deseq.csv
write.csv(resSig, foutName, row.names=T, quote=F)

## INIT
#foutName = paste(strsplit(cts_path, ".table")[1], "deseq.csv", sep=".") # ControlID_vs_CaseID.deseq.csv
#foutTable = paste(strsplit(cts_path, ".table")[1], "norm.csv", sep=".") # ControlID_vs_CaseID.norm.csv
#foutRawTable = paste(strsplit(cts_path, ".table")[1], "raw.csv", sep=".") # ControlID_vs_CaseID.raw.csv
#foutLog2Table = paste(strsplit(cts_path, ".table")[1], "normlog2.csv", sep=".") # ControlID_vs_CaseID.normlog2.csv
#tb = read.delim(cts_path,row.names=1)   # row.names="Symbol"
#print(cts_path)
##tb = read.delim(cts_path,row.names=NULL)   # row.names="Symbol"
#
#tb = trunc(tb)
#tb_g1 = tb[,1:controlNum]
#caseNumS = controlNum + 1
#caseNumE = controlNum + caseNum
#
#tb_g1 = tb[,1:controlNum] # 1:10 --> 10 controls
#tb_g2 = tb[,caseNumS:caseNumE] # 11:28 --> 18 cases
#rowMed1 = apply(tb_g1, 1, median)
#rowMed2 = apply(tb_g2, 1, median)
#rowTake = rep(TRUE, nrow(tb))
#for (i in 1:nrow(tb)) {
#    med1 = rowMed1[i]
#    med2 = rowMed2[i]
#    if (med1 < 10 && med2 < 10)
#        rowTake[i] = FALSE
#}
#
##rowMed = apply(tb, 1, median)
##med10_tb = tb[rowMed>=10,] # median more than 10
#med10_tb = tb[rowTake,] # median more than 10
#group = factor(c(rep(1,controlNum), rep(2,caseNum)))
#
### DESeq functions
#dds = DESeqDataSetFromMatrix(med10_tb, DataFrame(group), ~ group)
#dds = DESeq(dds)
#res = results(dds, addMLE=TRUE, alpha=0.05)
#resOrdered = res[order(res$pvalue),]
#resSig = subset(resOrdered, padj < padj_cut)
##resSig = subset(resOrdered, pvalue < padj_cut)
#summary(resSig)
#
#normalized_table = counts(dds, normalized=T)
#log2_norm_table = log2(normalized_table+1)
#
##write.table(normalized_table, foutTable, sep='\t', row.names=T, quote=F)
#write.csv(normalized_table, foutTable, row.names=T, quote=F)
#write.csv(log2_norm_table, foutLog2Table, row.names=T, quote=F)
#write.csv(med10_tb, foutRawTable, row.names=T, quote=F)
#
##write.table(resSig, foutName, sep='\t', row.names=T, quote=F)
#write.csv(resSig, foutName, row.names=T, quote=F)


