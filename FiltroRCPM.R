library(tximport)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(tidyverse)
edgeRUsersGuide()
#setwd("/media/alejandra/AgroBio_Tec/R/DEG/RSEM")
#setwd("/media/alejandra/AgroBio_Tec/ETAPA 2/R/DEG/RSEM")
setwd("F:/ETAPA 2/DEG/AnalisisR")
dir <- getwd()
run <- c("FH0_1_RHC.genes.results","FH0_2_RHC.genes.results","FH0_3_RHC.genes.results",
         "F0_1_RHC.genes.results","F0_2_RHC.genes.results","F0_3_RHC.genes.results",
         "FH24_1_RHC.genes.results","FH24_2_RHC.genes.results","FH24_3_RHC.genes.results",
         "F24_1_RHC.genes.results","F24_2_RHC.genes.results","F24_3_RHC.genes.results",
         "FH72_1_RHC.genes.results","FH72_2_RHC.genes.results","FH72_3_RHC.genes.results",
         "F72_1_RHC.genes.results","F72_2_RHC.genes.results","F72_3_RHC.genes.results")
         
files <- file.path(dir,run)
names(files) <- c("FC0_1","FC0_2","FC0_3",
                  "FF0_1","FF0_2","FF0_3",
                  "FC24_1","FC24_2","FC24_3",
                  "FF24_1","FF24_2","FF24_3",
                  "FC72_1","FC72_2","FC72_3",
                  "FF72_1","FF72_2","FF72_3" )

txi.rsem <- tximport(files, type= "rsem", txIn = F, txOut = F)
#write.csv(txi.rsem$counts,"Conteos-OK.csv")
class <- c( rep("FC0",3), 
            rep("FF0",3),
            rep("FC24",3), 
            rep("FF24",3),
            rep("FC72",3), 
            rep("FF72",3))

Mat <- as.data.frame(txi.rsem$counts, row.names = rownames(txi.rsem$counts))
coldata <- data.frame(row.names = colnames(Mat) ,class)
Gene <- Mat$gene_id



y <- DGEList(counts = Mat, group = class, genes = Gene)

keep = rowSums(cpm(y)>1) >= 2
y2 = y[keep, , keep.lib.sizes=FALSE]
y3 <- DGEList(counts = y2$counts, group = class, genes = rownames(y2$counts))
y3 <- calcNormFactors(y3)
y3 <- estimateCommonDisp(y3, verbose=TRUE)
y3 <- estimateTagwiseDisp(y3)

################################################################################
#                             T0
################################################################################
et = exactTest(y3, dispersion=y3$tagwise.dispersion , pair=c("FC0","FF0")) #Primero COntrol despu?s Tratay3iento.
eT0 <- et
y3T0 <- topTags(et, n=Inf)$table
yUDT0 = rownames(y3T0)[y3T0$FDR < 0.05 & abs(y3T0$logFC) >= 1]
y05T0 = rownames(y3T0)[y3T0$FDR < 0.05]
sum(y3T0$logFC>=1 & y3T0$FDR<0.05)
sum(y3T0$logFC<=-1 & y3T0$FDR<0.05)
sum(abs(y3T0$logFC)>=1 & y3T0$FDR<0.05)
summary(deT0 <- decideTestsDGE(eT0, adjust.method="BH", p=0.05))
sum(y3T0$FDR<0.05)

################################################################################
#                             T1
################################################################################
et = exactTest(y3, dispersion=y3$tagwise.dispersion , pair=c("FC1","FF1")) #Primero COntrol despu?s Tratay3iento.
eT1 <- et
y3T1 <- topTags(et, n=Inf)$table
yUDT1 = rownames(y3T1)[y3T1$FDR < 0.05 & abs(y3T1$logFC) >= 1]
y05T1 = rownames(y3T1)[y3T1$FDR < 0.05]
sum(y3T1$logFC>=1 & y3T1$FDR<0.05)
sum(y3T1$logFC<=-1 & y3T1$FDR<0.05)
sum(abs(y3T1$logFC)>=1 & y3T1$FDR<0.05)
summary(deT1 <- decideTestsDGE(eT1, adjust.method="BH", p=0.05))
sum(y3T1$FDR<0.05)

################################################################################
#                             T24
################################################################################
et = exactTest(y3, dispersion=y3$tagwise.dispersion , pair=c("FC24","FF24")) #Primero COntrol despu?s Tratay3iento.
eT24 <- et
y3T24 <- topTags(et, n=Inf)$table
yUDT24 = rownames(y3T24)[y3T24$FDR < 0.05 & abs(y3T24$logFC) >= 1]
y05T24 = rownames(y3T24)[y3T24$FDR < 0.05]
sum(y3T24$logFC>=1 & y3T24$FDR<0.05)
sum(y3T24$logFC<=-1 & y3T24$FDR<0.05)
sum(abs(y3T24$logFC)>=1 & y3T24$FDR<0.05)
summary(deT24 <- decideTestsDGE(eT24, adjust.method="BH", p=0.05))
sum(y3T24$FDR<0.05)

################################################################################
#                             T72
################################################################################
et = exactTest(y3, dispersion=y3$tagwise.dispersion , pair=c("FC72","FF72")) #Primero COntrol despu?s Tratay3iento.
eT72 <- et
y3T72 <- topTags(et, n=Inf)$table
yUDT72 = rownames(y3T72)[y3T72$FDR < 0.05 & abs(y3T72$logFC) >= 1]
y05T72 = rownames(y3T72)[y3T72$FDR < 0.05]
sum(y3T72$logFC>=1 & y3T72$FDR<0.05)
sum(y3T72$logFC<=-1 & y3T72$FDR<0.05)
sum(abs(y3T72$logFC)>=1 & y3T72$FDR<0.05)
summary(deT72 <- decideTestsDGE(eT72, adjust.method="BH", p=0.05))
sum(y3T72$FDR<0.05)

#______________________________________________________________________________#
#                       Volcano plot
#______________________________________________________________________________#
#Librer?as utilizadas
library(tidyverse)
library(ggrepel)

################################################################################
#                             T0
################################################################################
#Volcano plot
dataT0 <- y3T0 %>% 
  mutate(
    expression = case_when(logFC >= log(1) & FDR < 0.05 ~ "Reg +",
                           logFC <= -log(1) & FDR < 0.05 ~ "Reg -",
                           TRUE ~ "Sin cambio")
  )

#generar un volcano con los genes up y down
volc = ggplot(dataT0, aes(logFC, -log(FDR,10))) +
  geom_point(aes(color = expression), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_manual(values = c("red", "green", "gray50")) +
  guides(colour = guide_legend(override.aes = list(size=1.5)))+
  ggtitle("Volcano plot 0 horas")+
  geom_vline(xintercept =c(0), linetype = 1,
             col = "black")+
  geom_hline(yintercept =-log10(0.05), linetype = 1, col="black")
ggsave("VolcanoplotT0_lines_new_P.jpeg", device="jpeg") 

################################################################################
#                             T1
################################################################################
#Volcano plot
dataT1 <- y3T1 %>% 
  mutate(
    expression = case_when(logFC >= log(1) & FDR < 0.05 ~ "Reg +",
                           logFC <= -log(1) & FDR < 0.05 ~ "Reg -",
                           TRUE ~ "Sin cambio")
  )

#generar un volcano con los genes up y down
volc = ggplot(dataT1, aes(logFC, -log(FDR,10))) +
  geom_point(aes(color = expression), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_manual(values = c("red", "green", "gray50")) +
  guides(colour = guide_legend(override.aes = list(size=1.5)))+
  ggtitle("Volcano plot 0 horas")+
  geom_vline(xintercept =c(0), linetype = 1,
             col = "black")+
  geom_hline(yintercept =-log10(0.05), linetype = 1, col="black")
ggsave("VolcanoplotT1_lines_new_P.jpeg", device="jpeg") 

################################################################################
#                             T24
################################################################################
#Volcano plot
dataT24 <- y3T24 %>% 
  mutate(
    expression = case_when(logFC >= log(1) & FDR < 0.05 ~ "Reg +",
                           logFC <= -log(1) & FDR < 0.05 ~ "Reg -",
                           TRUE ~ "Sin cambio")
  )

#generar un volcano con los genes up y down
volc = ggplot(dataT24, aes(logFC, -log(FDR,10))) +
  geom_point(aes(color = expression), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_manual(values = c("red", "green", "gray50")) +
  guides(colour = guide_legend(override.aes = list(size=1.5)))+
  ggtitle("Volcano plot 0 horas")+
  geom_vline(xintercept =c(0), linetype = 1,
             col = "black")+
  geom_hline(yintercept =-log10(0.05), linetype = 1, col="black")
ggsave("VolcanoplotT24_lines_new_P.jpeg", device="jpeg") 

################################################################################
#                             T72
################################################################################
#Volcano plot
dataT72 <- y3T72 %>% 
  mutate(
    expression = case_when(logFC >= log(1) & FDR < 0.05 ~ "Reg +",
                           logFC <= -log(1) & FDR < 0.05 ~ "Reg -",
                           TRUE ~ "Sin cambio")
  )

#generar un volcano con los genes up y down
volc = ggplot(dataT72, aes(logFC, -log(FDR,10))) +
  geom_point(aes(color = expression), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_manual(values = c("red", "green", "gray50")) +
  guides(colour = guide_legend(override.aes = list(size=1.5)))+
  ggtitle("Volcano plot 0 horas")+
  geom_vline(xintercept =c(0), linetype = 1,
             col = "black")+
  geom_hline(yintercept =-log10(0.05), linetype = 1, col="black")
ggsave("VolcanoplotT72_lines_new_P.jpeg", device="jpeg") 


#______________________________________________________________________________#
#Imprimir tabla conr esultados de la expresi?n diferencial de los genes.
#______________________________________________________________________________#
write.csv(y3T0[y05T0,],"res05-T0.csv", row.names = F)
write.csv(y3T0[yUDT0,],"res05UD-T0.csv", row.names = F)

write.csv(y3T1[y05T1,],"res05-T1.csv", row.names = F)
write.csv(y3T1[yUDT1,],"res05UD-T1.csv", row.names = F)

write.csv(y3T24[y05T24,],"res05-T24.csv", row.names = F)
write.csv(y3T24[yUDT24,],"res05UD-T24.csv", row.names = F)

write.csv(y3T72[y05T72,],"res05-T72.csv", row.names = F)
write.csv(y3T72[yUDT72,],"res05UD-T72.csv", row.names = F)


####
# Tabla COnteos totales
####
RT0 <- as.data.frame(y3T0[y05T0,])
RT0$gene <- rownames(RT0)
resT0 <- RT0[,c(1,2)]

RT1 <- as.data.frame(y3T1[y05T1,])
RT1$gene <- rownames(RT1)
resT1 <- RT1[,c(1,2)]

RT24 <- as.data.frame(y3T24[y05T24,])
RT24$gene <- rownames(RT24)
resT24 <- RT24[,c(1,2)]

RT72 <- as.data.frame(y3T72[y05T72,])
RT72$gene <- rownames(RT72)
resT72 <- RT72[,c(1,2)]

####
# Tabla COnteos solo diferenciales con abs(LgFC)>=1
####
RT0 <- as.data.frame(y3T0[yUDT0,])
RT0$gene <- rownames(RT0)
resDFT0 <- RT0[,c(1,2)]

RT1 <- as.data.frame(y3T1[yUDT1,])
RT1$gene <- rownames(RT1)
resDFT1 <- RT1[,c(1,2)]

RT24 <- as.data.frame(y3T24[yUDT24,])
RT24$gene <- rownames(RT24)
resDFT24 <- RT24[,c(1,2)]

RT72 <- as.data.frame(y3T72[yUDT72,])
RT72$gene <- rownames(RT72)
resDFT72 <- RT72[,c(1,2)]


#Unificar todos los resultados de la expresi?n diferencial
TT1 <- full_join(resT0,resT24, "genes")
TT2 <- full_join(TT1, resT72, "genes")
write.csv(TT2,"Res05TOTAL.csv")

TTDF1 <- full_join(resDFT0,resDFT24, "genes")
TTDF2 <- full_join(TTDF1, resDFT72, "genes")
write.csv(TTDF2,"resDF05TOTALDF.csv")

KO <- as.data.frame(read.csv("Ko.csv"))
KoD <- inner_join(TT3,KO,"genes")
write.csv(KoD,"KOsas.csv")

#______________________________________________________________________________#
#                    Imprimir listas para Venn
#______________________________________________________________________________#

write( y3T0[y3T0$PValue < 0.05 & y3T0$logFC > 0,1], "UPT0.txt" )
write( y3T24[y3T24$PValue < 0.05 & y3T24$logFC > 0,1], "UPT24.txt" )
write( y3T72[y3T72$PValue < 0.05 & y3T72$logFC > 0,1], "UPT72.txt" )

write( y3T0[y3T0$PValue < 0.05 & y3T0$logFC < 0,1], "DNT0.txt" )
write( y3T24[y3T24$PValue < 0.05 & y3T24$logFC < 0,1], "DNT24.txt" )
write( y3T72[y3T72$PValue < 0.05 & y3T72$logFC < 0,1], "DNT72.txt" )

write( y3T0[y3T0$PValue < 0.05 & y3T0$logFC>=1,1], "UPUDT0.txt" )
write( y3T24[y3T24$PValue < 0.05 & y3T24$logFC>=1,1], "UPUDT24.txt" )
write( y3T72[y3T72$PValue < 0.05 & y3T72$logFC>=1,1], "UPUDT72.txt" )

write( y3T0[y3T0$PValue < 0.05 & y3T0$logFC<=-1,1], "DNuUDT0.txt" )
write( y3T24[y3T24$PValue < 0.05 & y3T24$logFC<=-1,1], "DNUDT24.txt" )
write( y3T72[y3T72$PValue < 0.05 & y3T72$logFC<=-1,1], "DNUDT72.txt" )

#______________________________________________________________________________#
#                    Asignacion de rutas metab?licas por PathView
#______________________________________________________________________________#
library(pathview)
KODEG <- as.data.frame(read.csv("KOxDEG.csv"))
testVect <- structure(KODEG[KODEG$T0!="NA",3], .Names =KODEG[KODEG$T0!="NA",2])
pv.out <- pathview(gene.data = testVect, pathway.id = "04075", species = "ko", out.suffix = "Pru.Fin", kegg.native = T)
str(pv.out)
head(pv.out$plot.data.gene)

