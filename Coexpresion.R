library(tidyverse)

setwd("F:/ETAPA 2/R/DEG/RSEM")

T0 <- as.data.frame(read.csv("res05-T0.csv"))
colnames(T0) <- c("transcript_id", "T0exp","logCPMT0","PValueT0","T0sig")
T24 <- as.data.frame(read.csv("res05-T24.csv"))
colnames(T24) <- c("transcript_id", "T24exp","logCPMT24","PValueT24","T24sig")
T72 <- as.data.frame(read.csv("res05-T72.csv"))
colnames(T72) <- c("transcript_id", "T72exp","logCPMT72","PValueT72","T72sig")

DB1 <- full_join(T0,T24, "transcript_id")
DB2 <- full_join(DB1, T72, "transcript_id")

head(DB2)
colnames(DB2)
DB3 <- DB2[,-c(3,4,7,8,11,12)]
head(DB3)
DB4 <- DB3[,c(1,2,4,6,3,5,7)]
head(DB4)
write.csv(DB4,"Coexpresion.csv")

Anot <- as.data.frame(read.csv("ReportexTiempoTOT_TrP_SP_T10.csv"))
colnames(Anot)
Anot2 <- Anot[,c(1,5,6,8,9,10,11,12,13)]
colnames(Anot2)
write.csv(Anot2,"AnotTOT.csv")

AnoTAIR <- as.data.frame(read.csv("TAIR10ANOT.csv"))
colnames(AnoTAIR)

FDB <- full_join(DB4,AnoTAIR,"transcript_id")
colnames(FDB)
FDB2 <- FDB[,c(1,8,2,3,4,5,6,7)]

colnames(FDB2)
write.csv(FDB2,"TAIRCoexp.csv")

AnotLocus <- as.data.frame(read.csv("AnotLocus.csv"))
colnames(AnotLocus)

FFV <- inner_join(FDB2,AnotLocus[,c(1,3)],"TAIR10pepcdna_BLASTX")
write.csv(FFV,"CoeTAIRF.csv")
