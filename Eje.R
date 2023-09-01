setwd("C:/Users/gusta/Downloads/Trinotate-Trinotate-v3.2.2/sample_data/data/DESeq2_gene")

if (! require(edgeR)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("edgeR")
  library(edgeR)
}

if (! require(DESeq2)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("DESeq2")
  library(DESeq2)
}

data = read.table("Trinity.gene.counts.matrix", header=T, row.names=1, com='')
head(data)
col_ordering = c(4,5,6,1,2,3)
rnaseqMatrix = data[,col_ordering]
head(rnaseqMatrix)
rnaseqMatrix = round(rnaseqMatrix)
head(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
conditions = data.frame(conditions=factor(c(rep("GSNO", 3), rep("wt", 3))))
rownames(conditions) = colnames(rnaseqMatrix)
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = rnaseqMatrix,
  colData = conditions,
  design = ~ conditions)
dds = DESeq(ddsFullCountTable)
contrast=c("conditions","GSNO","wt")
dds
head(dds)
res = results(dds, contrast)
baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "GSNO"])
baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "wt"])
res = cbind(baseMeanA, baseMeanB, as.data.frame(res))
res = cbind(sampleA="GSNO", sampleB="wt", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
res = as.data.frame(res[order(res$pvalue),])
write.table(res, file='Trinity.gene.counts.matrix.GSNO_vs_wt.DESeq2.DE_results', sep='	', quote=FALSE)
write.table(rnaseqMatrix, file='Trinity.gene.counts.matrix.GSNO_vs_wt.DESeq2.count_matrix', sep='	', quote=FALSE)
source("/seq/RNASEQ/tmp/trinityrnaseq-v2.10.0/Analysis/DifferentialExpression/R/rnaseq_plot_funcs.R")
pdf("Trinity.gene.counts.matrix.GSNO_vs_wt.DESeq2.DE_results.MA_n_Volcano.pdf")
plot_MA_and_Volcano(rownames(res), log2(res$baseMean+1), res$log2FoldChange, res$padj)
dev.off()
