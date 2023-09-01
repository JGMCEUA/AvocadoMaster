#para el analisis completo 
library(edgeR)
#metiendo datos
counts.all<-read.table(file="~/Documents/transcriptome_results/statistic_analysis/EdgeR/experiment_com/analisis2/count_all.txt",header=T, sep="\t",row.names=1)
#Factores
lines<-factor(c(rep(c("ws15","ws17","ws21","lec215","lec217","lec221","lec115","lec117","lec121","ler15","ler17","ler21","abi3115","abi3117","abi3121","abi3515", "abi3517","abi3521","col015","col017","col021","fus315","fus317","fus321"),2)))
d.all<- DGEList(counts.all,group=lines)#toda la matrix completa
d.all <- d.all[rowSums(d.all$counts) >= 5,]
d.all<- calcNormFactors(d.all)#Normalizacion
cpmd.all<-cpm(d.all,normalized.lib.sizes=TRUE)
write.table(cpmd.all, file ="cpm_norm.txt",sep="\t")

pdf(file="Desiccation_MDSplot.pdf")
plotMDS(d.all, main = "Desiccation tolerance", xlim = c(-3, 3))
dev.off()


#diseño
design <- model.matrix(~0+ lines)
d.all <- estimateGLMCommonDisp(d.all, design, verbose=TRUE)
d.all<-estimateGLMTrendedDisp(d.all, design)
d.all<- estimateGLMTagwiseDisp(d.all, design)

pdf(file="BCV_Desiccation.pdf")
plotBCV(d.all, main ="BCV")
dev.off()

colnames(design)<- levels(lines)
design

#haciendo contrastes
contVector <- c(
"lec1vsws_15"="lec115-ws15",
"lec1vsws_17"="lec117-ws17",
"lec1vsws_21"="lec121-ws21",
"lec1vsws"="(lec115-ws15)+(lec117-ws17)+(lec121-ws21)",
"fus3vscol0_15"="fus315-col015",
"fus3vscol0_17"="fus317-col017",
"fus3vscol0_21"="fus321-col021",
"fus3vscol0"="(fus315-col015)+(fus317-col017)+(fus321-col021)",
"lec1vswslec2_15"="(lec115-ws15)-(lec215-ws15)",
"lec1vswslec2_17"="(lec117-ws17)-(lec217-ws17)",
"lec1vswslec2_21"="(lec121-ws21)-(lec221-ws21)",
"lec1vswslec2"="((lec115-ws15)-(lec215-ws15))+((lec117-ws17)-(lec217-ws17))+((lec121-ws21)-(lec221-ws21))",
"abi35vslerabi31_15"= "(abi3515-ler15)-(abi3115-ler15)",
"abi35vslerabi31_17"= "(abi3517-ler17)-(abi3117-ler17)",
"abi35vslerabi31_21"= "(abi3521-ler21)-(abi3121-ler21)",
"abi35vslerabi31"="((abi3515-ler15)-(abi3115-ler15))+((abi3517-ler17)-(abi3117-ler17))+((abi3521-ler21)-(abi3121-ler21))",
"intolerantvstolerant_15"="((lec115-ws15)-(lec215-ws15))+((abi3515-ler15)-(abi3115-ler15))+(fus315-col015)",
"intolerantvstolerant_17"="((lec117-ws17)-(lec217-ws17))+((abi3517-ler17)-(abi3117-ler17))+(fus317-col017)",
"intolerantvstolerant_21"="((lec121-ws21)-(lec221-ws21))+((abi3521-ler21)-(abi3121-ler21))+(fus321-col021)",
"intolerantvstolerant"="(((lec115-ws15)-(lec215-ws15))+((abi3515-ler15)-(abi3115-ler15))+(fus315-col015))+(((lec117-ws17)-(lec217-ws17))+((abi3517-ler17)-(abi3117-ler17))+(fus317-col017))+(((lec121-ws21)-(lec221-ws21))+((abi3521-ler21)-(abi3121-ler21))+(fus321-col021))",
"difftolerant_17vsdifftolerant_15"="(((lec117-ws17)-(lec217-ws17))+((abi3517-ler17)-(abi3117-ler17))+(fus317-col017))-(((lec115-ws15)-(lec215-ws15))+((abi3515-ler15)-(abi3115-ler15))+(fus315-col015))",#intolerantvstolerant_17:intolerantvstolerant_15
"difftolerant_21vsdifftolerant_15"="(((lec121-ws21)-(lec221-ws21))+((abi3521-ler21)-(abi3121-ler21))+(fus321-col021))-(((lec115-ws15)-(lec215-ws15))+((abi3515-ler15)-(abi3115-ler15))+(fus315-col015))")#intolerantvstolerant_21:intolerantvstolerant_15


#"lec1vswslec2_15"="lec115-(ws15+lec215)/2",
#"lec1vswslec2_17"="lec117-(ws17+lec217)/2",
#"lec1vswslec2_21"="lec121-(ws21+lec221)/2",
#"abi35vslerabi31_15"= "abi3515-(ler15+abi3115)/2",
#"abi35vslerabi31_17"= "abi3517-(ler17+abi3117)/2",
#"abi35vslerabi31_21"= "abi3521-(ler21+abi3121)/2",
#"intoleranvstolerant_15"="(lec115-(ws15+lec215)/2)+(abi3515-(ler15+abi3115)/2)+(fus315-col015)",
#"intoleranvstolerant_17"="(lec117-(ws17+lec217)/2)+(abi3517-(ler17+abi3117)/2)+(fus317-col017)",
#"intoleranvstolerant_21"="(lec121-(ws21+lec221)/2)+(abi3521-(ler21+abi3121)/2)+(fus321-col021)",





contMatrix <- makeContrasts(contrasts=contVector, levels=design)
colnames(contMatrix) <- names(contVector)
#Haciendo el test
glmfit.all <- glmFit(d.all, design)
lrt.all <- glmLRT(glmfit.all, coef=18:20) #Para hacer una ANOVA :-D
topTags(lrt.all)
top <- rownames(topTags(lrt.all,n=50)$table)
print(d.all$counts[top,order(d.all$samples$group)])

for (cont in colnames(contMatrix)){
print(cont)
lrt.all<-glmLRT(glmfit.all, contrast=contMatrix[,cont])
print(topTags(lrt.all))
print(summary(decideTestsDGE(lrt.all, p.value=0.05)))
top<-topTags(lrt.all, n=nrow(d.all$counts))
write.table(top$table[,c(1,2,5)], file=paste(cont, ".txt", sep=""),sep="\t",col.names=TRUE)
#top <- rownames(topTags(lrt.all)$table)
#print(d.all$counts[top,order(d.all$samples$group)])
}


