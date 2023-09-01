library(tidyverse)

setwd("F:/ETAPA 2/R/DEG/RSEM")

T0  <- as.data.frame(read.csv("res05-T0.csv"))
T24 <- as.data.frame(read.csv("res05-T24.csv"))
T72 <- as.data.frame(read.csv("res05-T72.csv"))

ANOT <- as.data.frame(read.csv("Unicos.csv"))

colnames(ANOT)

ANOT[ANOT$T0 != "NA",]
TANOT0 <- inner_join(T0,ANOT[,1],"transcript_id")
head(TANOT0)

