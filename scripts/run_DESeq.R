#!/usr/local/bin/Rscript

## INPUT
library(DESeq2)

padj_cut <- 0.05 # P value cutoff per test
args <- commandArgs(trailingOnly = TRUE)

cts_path <- args[1]
coldata_path <- args[2]
out_dir <- args[3]

# Read input files
cts <- as.matrix(read.csv(cts_path, sep = '\t', row.names = 'gene'))
coldata <- read.csv(coldata_path, sep = '\t', row.names = 1)
coldata$condition <- factor(coldata$condition)

# Sanity check
if (!all(rownames(coldata) == colnames(cts))) {
  stop("Not all rownames of coldata match colnames of count matrix")
}

# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~condition)
dds <- DESeq(dds)

# Output: normalized counts
normalized_table <- counts(dds, normalized = TRUE)
base_name <- sub("\\.table.*", "", basename(cts_path))  # Extract base sample name
fout_norm <- file.path(out_dir, paste0(base_name, ".norm.csv"))
write.csv(normalized_table, fout_norm, row.names = TRUE, quote = FALSE)

# Output: DEG results
res <- results(dds, contrast = c("condition", "Case", "Control"))
resOrdered <- res[order(res$pvalue), ]
resSig <- subset(resOrdered, padj < padj_cut)

fout_deg <- file.path(out_dir, paste0(base_name, ".DEG.csv"))
write.csv(resSig, fout_deg, row.names = TRUE, quote = FALSE)
