# Análisis ballgown
library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

#Asignar carpeta
setwd("/home/alejandra/Escritorio/AlinTes2021")
pheno_Data <- read.csv("DRT0.csv")
list.files(path = "ballgown" ,pattern = "F")
bg_chrX <- ballgown(dataDir = "ballgown", samplePattern = "F", pData=pheno_Data)
bg_chrX_filt = subset(bg_chrX,"rowVars(texpr(bg_chrX)) >1", genomesubset = TRUE)
results_transcripts = stattest(bg_chrX_filt, feature="transcript", covariate= "Trat", meas = "FPKM")
results_genes = stattest(bg_chrX_filt, feature= "gene", covariate="Trat", meas= "FPKM")
results_transcripts = arrange(results_transcripts, pval)
results_genes = arrange(results_genes, pval)
write.csv(results_transcripts, "chrX_transcript_results_BT0.csv", row.names= FALSE)
write.csv(results_genes, "chreX_gene_results_BT0.csv", row.names = FALSE)