#!/bin/bash
cd /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2;
hisat2 --very-sensitive --phred33 --dta -p 8 -x /home/lbmv/Resp_070921_JGMC/rRNA_ref/SSU/hisat2/SSU_DNA_Silva --reorder --un-gz F1_1_UP_fail_f2.fq.gz --al-gz F1_1_UP_ali_f2.fq.gz --un-conc-gz F1_1_PO_fail_f2.fq.gz --al-conc-gz F1_1_PO_ali_f2.fq.gz --summary-file F1_1_res_f2.txt --met-file F1_1_Hisat_metrics_f2.txt --fr -1 /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro1/F1_1_PO_fail_f1.fq.1.gz -2 /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro1/F1_1_PO_fail_f1.fq.2.gz;
###Mover los resultados a la carpeta original
mv F1_1_PO_fail_f2.fq.2.gz F1_1_PO_fail_f2.fq.1.gz /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF;
cd /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF;
###Copiar carpetas para descompromirlos
cp F1_1_PO_fail_f2.fq.1.gz F1_1_PO_fail_f2.fq.2.gz F24_1_PO_fail_f2.fq.1.gz F24_1_PO_fail_f2.fq.2.gz F24_3_PO_fail_f2.fq.1.gz F24_3_PO_fail_f2.fq.2.gz /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble/readsOKF;
cd /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble/readsOKF;
### F0_1_PO_fail_f2.fq.1.gz -> F0_1_PO_fail_f2.R1.fq; se descomprimen y cambian los nombres para el análisis
gzip -d *;
mv  F1_1_PO_fail_f2.fq.1 F1_1_PO_fail_f2.R1.fq;
mv  F1_1_PO_fail_f2.fq.2 F1_1_PO_fail_f2.R1.fq;
mv  F24_1_PO_fail_f2.fq.1 F24_1_PO_fail_f2.R1.fq;
mv  F24_1_PO_fail_f2.fq.2 F24_1_PO_fail_f2.R2.fq;
mv  F24_3_PO_fail_f2.fq.1 F24_3_PO_fail_f2.R1.fq;
mv  F24_3_PO_fail_f2.fq.2 F24_3_PO_fail_f2.R2.fq;
### Concatenar datos que faltan c(72) ->  F72_PO_trinity_R1.fq
cat F1_[123]_PO_fail_f2.R1.fq >> F1_PO_trinity_R1.fq | cat F1_[123]_PO_fail_f2.R2.fq >> F1_PO_trinity_R2.fq;
cat F24_[123]_PO_fail_f2.R1.fq >> F24_PO_trinity_R1.fq | cat F24_[123]_PO_fail_f2.R2.fq >> F24_PO_trinity_R2.fq;
mv F1_PO_trinity_R1.fq F24_PO_trinity_R1.fq F1_PO_trinity_R2.fq F24_PO_trinity_R2.fq /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble;
cd /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble;
###Inicio de ensamble
var=$(date); echo "Comenzó ensamble de F1c[ $var]">> "ResumenE.txt";
$TRINITY_HOME/Trinity --seqType fq --max_memory 120G --left F1_PO_trinity_R1.fq --right F1_PO_trinity_R2.fq --SS_lib_type FR --CPU 10 --output /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble/F1Ctrinity --normalize_max_read_cov 100 &>>F1C_TR.txt;
var=$(date); echo "Teminó ensamble de F1c[ $var]">> "ResumenE.txt";
var=$(date); echo "Comenzó ensamble de F24c[ $var]">> "ResumenE.txt";
$TRINITY_HOME/Trinity --seqType fq --max_memory 120G --left F24_PO_trinity_R1.fq --right F24_PO_trinity_R2.fq --SS_lib_type FR --CPU 10 --output /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble/F24Ctrinity --normalize_max_read_cov 100 &>>F24C_TR.txt;
var=$(date); echo "Teminó ensamble de F24c[ $var]">> "ResumenE.txt";
mv F1_PO_trinity_R1.fq F24_PO_trinity_R1.fq F1_PO_trinity_R2.fq F24_PO_trinity_R2.fq /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble/readsOKF;
cd /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble/readsOKF;
###Ensamble de F
var=$(date); echo "Comenzó concatenación de Fc Total[ $var]">> "Resumen.txt";
cat F0_PO_trinity_R1.fq F1_PO_trinity_R1.fq F24_PO_trinity_R1.fq F72_PO_trinity_R1.fq >> Fc_PO_trinity_R1.fq;
cat F0_PO_trinity_R2.fq F1_PO_trinity_R2.fq F24_PO_trinity_R2.fq F72_PO_trinity_R2.fq >> Fc_PO_trinity_R2.fq;
mv Fc_PO_trinity_R1.fq Fc_PO_trinity_R2.fq /home/lbmv/EnsambleF;
var=$(date); echo "Terminó concatenación de Fc Total[ $var]">> "Resumen.txt";
###Concatenación Final suprema
cd /home/lbmv/EnsambleF;
cat FH_PO_trinity_R1.fq Fc_PO_trinity_R1.fq >> Persea_americana_Hass_POR1c.fq;
cat FH_PO_trinity_R2.fq Fc_PO_trinity_R2.fq >> Persea_americana_Hass_POR2c.fq;
###Ensamblado supremo
var=$(date); echo "Comenzó ensamble de P.A.Hassc[ $var]">> "ResumenE.txt";
$TRINITY_HOME/Trinity --seqType fq --max_memory 120G --left Persea_americana_Hass_POR1c.fq --right Persea_americana_Hass_POR2c.fq --SS_lib_type FR --CPU 10 --output /home/lbmv/EnsambleF/Persea_americana_HassCtrinity --normalize_max_read_cov 50 &>> P_A_H_TRC.txt;
var=$(date); echo "Terminó ensamble de P.A.Hassc[ $var]">> "ResumenE.txt";
#Ensamblado independiente F
cd /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble/readsOKF;
var=$(date); echo "Comenzó ensamble de F0_ coma[ $var]">> "ResumenE.txt";
$TRINITY_HOME/Trinity --seqType fq --max_memory 120G --left F0_1_PO_fail_f2.R1.fq,F0_2_PO_fail_f2.R1.fq,F0_3_PO_fail_f2.R1.fq --right F0_1_PO_fail_f2.R2.fq,F0_2_PO_fail_f2.R2.fq,F0_3_PO_fail_f2.R2.fq --SS_lib_type FR --CPU 10 --output /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble/F0comatrinity --normalize_max_read_cov 100 &>>F0_TR_coma.txt;
var=$(date); echo "Teminó ensamble de F0_coma[ $var]">> "ResumenE.txt";
cd /home/lbmv;
echo "Por fin, todo acabó, pero a que precio" >> mensaje.txt;
