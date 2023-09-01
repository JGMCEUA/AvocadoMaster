library("ggVennDiagram")
library("ggvenn")
library("tidyverse")

#install.packages("ggVennDiagram")
setwd("G:/AgroBio_4TB/ETAPA 2/DEG/AnalisisR")
#setwd("/media/alejandra/AgroBio_Tec/R/DEG/Venn")

T0 <- as.data.frame(read.csv("res05-T0.csv"))
T24 <- as.data.frame(read.csv("res05-T24.csv"))
T72 <- as.data.frame(read.csv("res05-T72.csv"))
#RES <- list(c("T0up",T0[T0$logFC>0,1]),c("T0dn",T0[T0$logFC<0,1]), c("T24up",T24[T24$logFC>0,1]),c("T24dn",T24[T24$logFC<0,1]), c("T72up",T72[T72$logFC>0,1]),c("T72dn",T72[T72$logFC<0,1]))
#head(RES)
write(T0[T0$logFC>=1,1], "UPT0_DF.txt")
write(T0[T0$logFC<=-1,1], "DNT0_DF.txt")

write(T24[T24$logFC>=1,1], "UPT24_DF.txt")
write(T24[T24$logFC<=-1,1], "DNT24_DF.txt")

write(T72[T72$logFC>=1,1], "UPT72_DF.txt")
write(T72[T72$logFC<=-1,1], "DNT72_DF.txt")

write(T0[,1], "T0DG.txt")
write(T24[,1], "T24DG.txt")
write(T72[,1], "T72DG.txt")

RES <- as.data.frame(read.csv("vennDEG.csv"))
head(RES)
#xu<- list(T0 = subset(RES$T0up, subset = RES$T0up!="") ,T1 = subset(RES$T1up, subset = RES$T1up!=""),T24 = subset(RES$T24up, subset = RES$T24up!=""),T72 = subset(RES$T72up, subset = RES$T72up!=""))
#xd<- list(T0 = subset(RES$T0dn, subset = RES$T0dn!="") ,T1 = subset(RES$T1dn, subset = RES$T1dn!=""),T24 = subset(RES$T24dn, subset = RES$T24dn!=""),T72 = subset(RES$T72dn, subset = RES$T72dn!=""))
#x <- list(T0 = c(xu$T0,xd$T0),T1 = c(xu$T1,xd$T1),T24 = c(xu$T24,xd$T24),T72 = c(xu$T72,xd$T72))
#ggVennDiagram(x) + ggplot2::scale_fill_gradient(low="#2C9236", high = "#04D613") + ggplot2::geom_contour(color = "Black")
xu<- list("1 hpt" = subset(RES$T0up, subset = RES$T0up!="") ,"24 hpt" = subset(RES$T24up, subset = RES$T24up!=""),"72 hpt" = subset(RES$T72up, subset = RES$T72up!=""))
xd<- list("1 hpt" = subset(RES$T0dn, subset = RES$T0dn!="") ,"24 hpt" = subset(RES$T24dn, subset = RES$T24dn!=""),"72 hpt" = subset(RES$T72dn, subset = RES$T72dn!=""))
x <- list("1 hpt" = c(xu$`1 hpt`,xd$`1 hpt`),"24 hpt" = c(xu$`24 hpt`,xd$`24 hpt`),"72 hpt" = c(xu$`72 hpt`,xd$`72 hpt`))

ggvenn(x, fill_color = c("#2C9236", "#AF7F3F", "#868686FF", "#0EB905"), stroke_size = 0.5, set_name_size = 4, show_elements = F, ) + ggtitle("Transcritos totales") + theme(plot.title = element_text(hjust = 0.5)) 
  
ggvenn(xu, fill_color = c("#2C9236", "#AF7F3F", "#868686FF", "#0EB905"), stroke_size = 0.5, set_name_size = 4)+ ggtitle("Transcritos regulados positivamente") + theme(plot.title = element_text(hjust = 0.5))
ggvenn(xd, fill_color = c("#2C9236", "#AF7F3F", "#868686FF", "#0EB905"), stroke_size = 0.5, set_name_size = 4)+ ggtitle("Transcritos regulados negativamente") + theme(plot.title = element_text(hjust = 0.5))


cGEFTT <- as.data.frame(read.csv("DEGTT.csv"))
colnames(GEFTT) <- c("T0",colnames(GEFTT[,-1]))
colnames(GEFTT)

t04 = subset(GEFTT[,1], subset = GEFTT[,1]!="")
t03 = subset(GEFTT[,7], subset = GEFTT[,7]!="")
t243 = subset(GEFTT[,9], subset = GEFTT[,9]!="")
t244 = subset(GEFTT[,3], subset = GEFTT[,3]!="")
t724 = subset(GEFTT[,5], subset = GEFTT[,5]!="")
t723 = subset(GEFTT[,11], subset = GEFTT[,11]!="")

t0 <- list(T04 = t04, T03 = t03)
t24 <- list(T244=t244,T243=t243)
t72 <- list(T724=t724,T723=t723)

ggvenn(t0, fill_color = c("#2c9236", "#af7f3f", "#868686ff", "#0eb905"), stroke_size = 0.5, set_name_size = 4)
ggvenn(t24, fill_color = c("#2c9236", "#af7f3f", "#868686ff", "#0eb905"), stroke_size = 0.5, set_name_size = 4)
ggvenn(t72, fill_color = c("#2c9236", "#af7f3f", "#868686ff", "#0eb905"), stroke_size = 0.5, set_name_size = 4)

write.csv(full_join(GEFTT[c(1,2)],GEFTT[c(7,8)],"T0"),"ListT0.csv")
write.csv(full_join(GEFTT[c(3,4)],GEFTT[c(9,10)],"T24"),"ListT24.csv")
write.csv(full_join(GEFTT[c(5,6)],GEFTT[c(11,12)],"T72"),"ListT72.csv")

