## Ejercicio de anotacion
##Georgina Hernandez Montes
## viernes 30 de agosto 2021
## ejemplo tomado de https://www.bioconductor.org/packages/release/workflows/vignettes/annotation/inst/doc/Annotation_Resources.html
## Genomic Annotation Resources

##Cargar la paqueteria

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationHub")
library(AnnotationHub)   
### generar un objeto con las bases de datos
ah <- AnnotationHub()
ah
## exploracion del objeto
unique(ah$dataprovider) ###proveedores
unique(ah$rdataclass)  ### tipo de objeto
orgs <- subset(ah, ah$rdataclass == "OrgDb")
orgs

### seleccion de la base de datos de una especie
dog <- query(orgs, "Canis familiaris")[[1]]
dog
columns(dog)
keytypes(dog)
head(keys(dog, keytype="ENTREZID"))
head(keys(dog, keytype="SYMBOL"))
## seleccion de genes de interes
keys(dog, keytype="SYMBOL", pattern="COX")
keys(dog, keytype="ENTREZID", pattern="COX", column="SYMBOL")
keys(dog, keytype="ENSEMBL", pattern="COX", column="SYMBOL")
select(dog, keys="804478", columns=c("SYMBOL","REFSEQ"), keytype="ENTREZID")
select(dog, keys="804478", columns="GO", keytype="ENTREZID")
############
keys(dog, keytype="SYMBOL", pattern="AKT")
keys(dog, keytype="ENTREZID", pattern="AKT", column="SYMBOL")
ids<-keys(dog, keytype="ENTREZID", pattern="AKT", column="SYMBOL")
select(dog, keys=ids, columns=c("GENENAME","SYMBOL", "ONTOLOGY"), keytype="ENTREZID")
tabla_anot<-select(dog, keys=ids, columns=c("GENENAME","SYMBOL", "ONTOLOGY"), keytype="ENTREZID")
tabla_anot
write.table(tabla_anot,file="C:/Users/User/Desktop/Material_R_cursos/Clase_agosto_2021/misresultados.csv", sep= "\t")

sessionInfo()
###################
# Ejercicio
# 1. Construir el objeto hs donde guarden la base de datos para Homo Sapiens
# 2. Sobre ese objeto buscar los genes que correspondan al patron ATG
# 3 Crear un objeto con los identificadores ENTREZID
# 3. Crear una tabla con ENTREZID, SYMBOL, GENENAME y PATH y exportarla a su carpeta