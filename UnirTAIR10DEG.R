library(tidyverse)
setwd("F:/ETAPA 2/R/DEG/RSEM")

Seq <- as.data.frame(read.csv("Seq.csv"))
Tab <- as.data.frame(read.csv("Genes Gustavo.csv"))

colnames(Seq)
colnames(Tab)

D <- inner_join(Tab,Seq,"transcript_id")
colnames(D)

write.csv(D,"TablaConc.csv")
