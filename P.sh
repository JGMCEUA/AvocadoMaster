#!/bin/sh
#Prueba con Bowtie
cd "/home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/bowtie";
var=$(date); echo "Comenzó Alineamiento LSU F f1 [ $var]">> "Resumen.txt";
bowtie2 --quiet --very-sensitive-local --phred33  -x /home/lbmv/Resp_070921_JGMC/rRNA_ref/LSU/bowtie2/LSU_DNA_Silva_1 -1 /home/lbmv/Resp_070921_JGMC/FrutoTratado/seqCorre/unfixrm_F0_1_R1_T_PO.cor.fq.gz -2 /home/lbmv/Resp_070921_JGMC/FrutoTratado/seqCorre/unfixrm_F0_1_R2_T_PO.cor.fq.gz --threads 8 --met-file F0_1_bowtie2_metrics.txt --al-conc-gz blacklist_paired_aligned_F0_1.fq.gz --un-conc-gz blacklist_paired_unaligned_F0_1.fq.gz --al-gz blacklist_unpaired_aligned_F0_1.fq.gz --un-gz blacklist_unpaired_unaligned_F0_1.fq.gz
var=$(date); echo "Terminó Alineamiento LSU F f1 [ $var]">> "Resumen.txt";
#Prueba con sortmerna
cd "/home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/sortmerna";
var=$(date); echo "Comenzó Alineamiento LSU F f1 [ $var]">> "Resumen.txt";
sortmerna --ref /home/lbmv/Resp_070921_JGMC/rRNA_ref/LSU/DNA_SILVA_138.1_LSUParc_tax_silva.fasta --reads /home/lbmv/Resp_070921_JGMC/FrutoTratado/seqCorre/unfixrm_F0_1_R1_T_PO.cor.fq.gz /home/lbmv/Resp_070921_JGMC/FrutoTratado/seqCorre/unfixrm_F0_1_R2_T_PO.cor.fq.gz >> Res
var=$(date); echo "Terminó Alineamiento LSU F f1 [ $var]">> "Resumen.txt";
#Descompresión de archivos
cd /home/lbmv/JGMC_copias/FrutoHass/filtro2/filtradosFH;
var=$(date); echo "Comenzó duplicado FH[ $var]">> "Resumen.txt";
cp * /home/lbmv/JGMC_copias/FrutoHass/filtro2/filtradosFH/ensamble
var=$(date); echo "Terminó duplicado FH[ $var]">> "Resumen.txt";
cd /home/lbmv/JGMC_copias/FrutoHass/filtro2/filtradosFH/ensamble
var=$(date); echo "Comenzó descompresión FH[ $var]">> "Resumen.txt";
gzip -d *fq.1.gz | gzip -d *fq.2.gz;
var=$(date); echo "Terminó descompresión FH[ $var]">> "Resumen.txt";
#Aqui empieza el del comando para Fruto Tratado
cd /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2;
mv F*PO_ali_f2*.fq /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF;
cd /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF;
var=$(date); echo "Comenzó duplicado F[ $var]">> "Resumen.txt";
cp * /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble;
var=$(date); echo "Terminó duplicado F[ $var]">> "Resumen.txt";
cd /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble;
var=$(date); echo "Comenzó descompresión F[ $var]">> "Resumen.txt";
gzip -d *fq.1.gz | gzip -d *fq.2.gz;
var=$(date); echo "Terminó descompresión F[ $var]">> "Resumen.txt";
