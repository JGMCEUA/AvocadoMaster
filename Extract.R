library(readr)
library(tidyverse)
library(dslabs)

################################################################################
#                     RT0P
################################################################################
setwd("E:/R/NOCONTAR/CrearGTF/PruebaPy/RT0P")
geidr <- read.csv("ATHxPAH.csv")

gtf <- rtracklayer::import("PAH_ittepic.gtf")
gtf_df<-as.data.frame(gtf)
gtf_df <- subset(gtf_df, subset = gtf_df$type != "exon")


head(gtf_df)
head(geidr)

Res <- inner_join(geidr,gtf_df, by = "gene_id")
write.csv(Res,"Pruebasini.csv")


################################################################################
#                     RT1PPP
################################################################################