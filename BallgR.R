library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

setwd("G:/CONTEOS")
dir <- getwd()
#run <- c( "F0_1_2023_C","F0_2_2023_C","F0_3_2023_C","F24_1_2023_C","F24_2_2023_C","F24_3_2023_C","F72_1_2023_C","F72_2_2023_C","F72_3_2023_C","FH0_1_2023_C","FH0_2_2023_C","FH0_3_2023_C","FH24_1_2023_C","FH24_2_2023_C","FH24_3_2023_C","FH72_1_2023_C","FH72_2_2023_C","FH72_3_2023_C")
run <- c( "F0_1_2023_C","F0_2_2023_C","F0_3_2023_C","FH0_1_2023_C","FH0_2_2023_C","FH0_3_2023_C")
Direct <- file.path(dir,"ballgown",run)
pheno_data = read.csv("Descripcion0T.csv")
bg_chrX = ballgown( samples = Direct, pData=pheno_data)

bg_chrX_filt = subset(bg_chrX,"rowVars(texpr(bg_chrX)) >1",genomesubset=TRUE)
results_transcripts = stattest(bg_chrX_filt, feature="transcript",covariate="Tratamiento", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg_chrX_filt, feature="gene",covariate="Tratamiento", getFC=TRUE,meas="FPKM")                      
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_chrX_filt), geneIDs=ballgown::geneIDs(bg_chrX_filt), transcript_id= ballgown::transcriptNames(bg_chrX_filt), results_transcripts, LogFC=log2(results_transcripts$fc))
results_transcripts = arrange(results_transcripts,pval)
results_genes = arrange(results_genes,pval)


write.csv(results_transcripts[results_transcripts$pval<=0.05,], "chrX_transcript_results0T.csv", row.names=FALSE)
write.csv(results_genes[results_genes$pval<=0.05,], "chrX_gene_results0T.csv",row.names=FALSE)
subset(results_transcripts,results_transcripts$qval<0.05)
subset(results_genes,results_genes$qval<0.05)
tropical= c('darkorange', 'dodgerblue','hotpink', 'limegreen', 'yellow')
palette(tropical)

run24 <- c( "F24_1_2023_C","F24_2_2023_C","F24_3_2023_C","FH24_1_2023_C","FH24_2_2023_C","FH24_3_2023_C")
Direct <- file.path(dir,"ballgown",run24)
pheno_data = read.csv("Descripcion24T.csv")
bg_chrX = ballgown( samples = Direct, pData=pheno_data)
bg_chrX_filt = subset(bg_chrX,"rowVars(texpr(bg_chrX)) >1",genomesubset=TRUE)
results_transcripts = stattest(bg_chrX_filt, feature="transcript",covariate="Tratamiento", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg_chrX_filt, feature="gene",covariate="Tratamiento", getFC=TRUE,meas="FPKM")                      
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_chrX_filt), geneIDs=ballgown::geneIDs(bg_chrX_filt), transcript_id= ballgown::transcriptNames(bg_chrX_filt), results_transcripts, LogFC=log2(results_transcripts$fc))
results_transcripts = arrange(results_transcripts,pval)
results_genes = arrange(results_genes,pval)
write.csv(results_transcripts[results_transcripts$pval<=0.05,], "chrX_transcript_results24T.csv", row.names=FALSE)
write.csv(results_genes[results_genes$pval<=0.05,], "chrX_gene_results24T.csv",row.names=FALSE)


run72 <- c( "F72_1_2023_C","F72_2_2023_C","F72_3_2023_C","FH72_1_2023_C","FH72_2_2023_C","FH72_3_2023_C")
Direct <- file.path(dir,"ballgown",run72)
pheno_data = read.csv("Descripcion72T.csv")
bg_chrX = ballgown( samples = Direct, pData=pheno_data)
bg_chrX_filt = subset(bg_chrX,"rowVars(texpr(bg_chrX)) >1",genomesubset=TRUE)
results_transcripts = stattest(bg_chrX_filt, feature="transcript",covariate="Tratamiento", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg_chrX_filt, feature="gene",covariate="Tratamiento", getFC=TRUE,meas="FPKM")                      
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_chrX_filt), geneIDs=ballgown::geneIDs(bg_chrX_filt), transcript_id= ballgown::transcriptNames(bg_chrX_filt), results_transcripts, LogFC=log2(results_transcripts$fc))
results_transcripts = arrange(results_transcripts,pval)
results_genes = arrange(results_genes,pval)
write.csv(results_transcripts[results_transcripts$pval<=0.05,], "chrX_transcript_results72T.csv", row.names=FALSE)
write.csv(results_genes[results_genes$pval<=0.05,], "chrX_gene_results72T.csv",row.names=FALSE)
