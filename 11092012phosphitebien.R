#metiendo datos
counts.all<-read.delim(file="library.txt",stringsAsFactors=FALSE)
d.all<- readDGE(counts.all, skip = 5, comment.char = "!")#toda la matrix completa
colnames(d.all$counts)<-counts.all[,3]
write.table(d.all$counts, file="count_all.txt",sep="\t")
d.all <- d.all[rowSums(d.all$counts) >= 3, ]
d.all<- calcNormFactors(d.all)#Normalizacion
cpmd.all<-cpm(d.all,normalized.lib.sizes=TRUE)
write.table(cpmd.all, file ="cpm_norm.txt",sep="\t")


#Analisis completo
d.all<- readDGE(counts.all, skip = 5, comment.char = "!")#toda la matrix completa
colnames(d.all$counts)<-counts.all[,3]
d.all <- d.all[rowSums(d.all$counts) >= 3, ]
d.all<- calcNormFactors(d.all)#Normalizacion
pdf(file="phosphiteMDSplot.pdf")
plotMDS(d.all, main = "phosphite", xlim = c(-3, 3))
dev.off()

#Factores
treatment<-factor(c("MaP","MaP","MaP","MaP","MeP","MeP","MeP","MeP"))
time<-factor(rep(c("4H","8D"),4))
lines<-factor(c("MaP4H","MaP8D","MaPPhi4H","MaPPhi8D","MeP4H","MeP8D","MePPhi4H","MePPhi8D"))


#diseño
design <- model.matrix(~0+ lines)
colnames(design)<- levels(lines)

contVector <- c(
"MaPvsMeP_4H"="MeP4H-MaP4H", 
"MaPvsMeP_8D"="MeP8D-MaP8D",
"MaPvsMePPhi_4H"="MePPhi4H-MaP4H",
"MaPvsMePPhi_8D"="MePPhi8D-MaP8D",
"MePvsMePPhi_4H"="MePPhi4H-MeP4H",#genes reprimidos en presencia de Phi
"MePvsMePPhi_8D"="MePPhi8D-MeP8D",
"MaPvsMaPPhi_4H"="MaPPhi4H-MaP4H",
"MaPvsMaPPhi_8D"="MaPPhi8D-MaP8D",
"4Hvs8D_MeP"="(MeP8D-MeP4H)-(MaP8D-MaP4H)",
"4Hvs8D_MePPhi"="(MePPhi8D-MePPhi4H)-(MaP8D-MaP4H)",
"8DMeP4HMeP"="(MeP8D-MaP8D)+(MePPhi4H-MeP4H)", #genes que se reprimen en Phivs -P a 4H y 8D en -P
"MaPvsMeP"= "(MeP4H-MaP4H)+(MeP8D-MaP8D)+(MePPhi4H-MaP4H)+(MePPhi8D-MaP8D)")#+(MePPhi8D-MeP8D)+(MePPhi4H-MeP4H)+(MeP8D-MePPhi8D)+(MeP4H-MePPhi4H)")#genes responsivos a -P
#"4Hvs8D_MePPhi"="(MeP4H-MePPhi4H)-(MeP8D-MePPhi8D)",#Genes diferenciales que solo reponden a phi en carencia a pi en tiempos cortos
#"MaPPhivsMePPhi_4H"="(MaP4H-MaPP

hi4H)-(MeP4H-MePPhi4H)",
#"MaPPhivsMePPhi_8D"="(MaP8D-MaPPhi8D)-(MeP8D-MePPhi8D)")
#"4Hvs8D_MeP"="(MaP4H-MeP4H)-(MaP4H-MeP8D)",#Genes diferenciales que se comparten en tiempo cortos y largos en carencia de pi
#"4Hvs8D_MePPhi"="(MaP4H-MePPhi4H)-(MaP4H-MePPhi8D)",#Genes diferenciales que solo responden a phi en tiempos cortos
contMatrix <- makeContrasts(contrasts=contVector, levels=design)
colnames(contMatrix) <- names(contVector)

#Haciendo el test
glmfit.all <- glmFit(d.all, design,dispersion=0.05)#para la corrida2 de mi experimento
lrt.all <- glmLRT(d.all, glmfit.all, contrast=contMatrix[,"8DMeP4HMeP"])
topTags(lrt.all,n=50)
summary(decideTestsDGE(lrt.all, p.value=0.05))
top <- rownames(topTags(lrt.all,n=50)$table)
print(d.all$counts[top,order(d.all$samples$group)])


for (cont in colnames(contMatrix)){
print(cont)
lrt.all<-glmLRT(d.all, glmfit.all, contrast=contMatrix[,cont])
print(topTags(lrt.all))
print(summary(decideTestsDGE(lrt.all, p.value=0.05)))
top<-topTags(lrt.all, n=nrow(d.all$counts))
#write.table(top$table[,c(1,4,5)], file=paste(cont, ".txt", sep=""),sep="\t",col.names=TRUE)
#top <- rownames(topTags(lrt.all)$table)
#print(d.all$counts[top,order(d.all$samples$group)])
}




