counts.all<-read.table(file="/home/sandra/Documents/transcriptome_results/statistic_analysis/EdgeR/phosphite/count_all.txt",header=T, sep="\t",row.names=1)

#groups
samples<-factor(c("MaP4H","MaP8D","MaPPhi4H","MaPPhi8D","MeP4H","MeP8D","MePPhi4H","MePPhi8D"))
#Analisis completo
d.all<- DGEList(counts.all,group=samples)#toda la matrix completa
d.all <- d.all[rowSums(d.all$counts) >= 3, ]
d.all<- calcNormFactors(d.all)#Normalizacion


#diseño
design <- model.matrix(~0+ lines)
colnames(design)<- levels(lines)
contVector <- c(
"diffMePandMaPPhi_8D"= "(MeP8D-MaP8D)-(MaPPhi8D-MaP8D)", #genes que responden diferencialmente a carencia de pi y en presencia de phi 8D
"commMePandMaPPhi_8D"= "((MeP8D-MaP8D)+(MaPPhi8D-MaP8D))", #genes que responden igual a carencia de pi y en presencia de phi 8D
"commMePandMaPPhi_8D2"= "((MeP8D-MaP8D)+(MaPPhi8D-MaP8D)+(MeP8D-MePPhi8D))", 
"Resp_Phi_8D"="(MaPPhi8D-MaP8D)+(MePPhi8D-MeP8D)",
"Resp_local1"= "(MeP4H-MaP4H)+(MeP4H-MePPhi4H)-(MeP8D-MaP8D)",
"Resp_local2"="(MeP4H-MaP4H)+(MeP4H-MePPhi4H)",
"Resp_Phi_4H"="(MaPPhi4H-MaP4H)+(MePPhi4H-MeP4H)",
"RespMavsMePPhi_4H"="(MePPhi4H-MaP4H)-((MeP4H-MaP4H)+(MeP4H-MePPhi4H))",
"MePPhivsMaP_4H"="MePPhi4H-MaP4H",
"MaPvsMeP_4H"="MeP4H-MaP4H", 
"MaPvsMeP_8D"="MeP8D-MaP8D",
"MaPvsMePPhi_4H"="MePPhi4H-MaP4H",
"MaPvsMePPhi_8D"="MePPhi8D-MaP8D",
"MaPvsMaPPhi_4H"="MaPPhi4H-MaP4H",
"MaPvsMaPPhi_8D"="MaPPhi8D-MaP8D",
"MePvsMePPhi_8D"="MePPhi8D-MeP8D",
"MaPvsMaPPhi_4H"="MaPPhi4H-MaP4H")

contMatrix <- makeContrasts(contrasts=contVector, levels=design)
colnames(contMatrix) <- names(contVector)

#Haciendo el test
glmfit.all <- glmFit(d.all, design,dispersion=0.05)#para la corrida2 de mi experimento

for (cont in colnames(contMatrix)){
print(cont)
lrt.all<-glmLRT(d.all, glmfit.all, contrast=contMatrix[,cont])
print(topTags(lrt.all))
print(summary(decideTestsDGE(lrt.all, p.value=0.01)))
top<-topTags(lrt.all, n=nrow(d.all$counts))
write.table(top$table[,c(1,5)], file=paste(cont, ".txt", sep=""),sep="\t",col.names=TRUE)
top <- rownames(topTags(lrt.all)$table)
print(d.all$counts[top,order(d.all$samples$group)])
}


#"MaPvsMeP_4H"="MeP4H-MaP4H", 
#"MaPvsMeP_8D"="MeP8D-MaP8D",
#"MaPvsMePPhi_4H"="MePPhi4H-MaP4H",
#"MaPvsMePPhi_8D"="MePPhi8D-MaP8D",
#"MePvsMePPhi_4H"="MePPhi4H-MeP4H",#genes reprimidos en presencia de Phi
#"MePvsMePPhi_8D"="MePPhi8D-MeP8D",
#"MaPvsMaPPhi_4H"="MaPPhi4H-MaP4H",
#"MaPvsMaPPhi_8D"="MaPPhi8D-MaP8D",

