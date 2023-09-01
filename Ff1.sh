#!/bin/sh
cd /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble/readsOKF;
cat F1_[123]_PO_fail_f2.R1.fq >> F1_PO_trinity_R1.fq | cat F1_[123]_PO_fail_f2.R2.fq >> F1_PO_trinity_R2.fq;
mv F1_PO_trinity_R1.fq F1_PO_trinity_R2.fq /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble;
#t ###Inicio de ensamble
#t var=$(date); echo "Comenzó ensamble de F1c[ $var]">> "ResumenE.txt";
#t $TRINITY_HOME/Trinity --seqType fq --max_memory 120G --left F1_PO_trinity_R1.fq --right F1_PO_trinity_R2.fq --SS_lib_type FR --CPU 10 --output /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble/F1Ctrinity --normalize_max_read_cov 100 &>>F1C_TR.txt;
#t var=$(date); echo "Teminó ensamble de F1c[ $var]">> "ResumenE.txt";
mv F1_PO_trinity_R1.fq F1_PO_trinity_R2.fq /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble/readsOKF;
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
#t ###Ensamblado supremo
#t var=$(date); echo "Comenzó ensamble de P.A.Hassc[ $var]">> "ResumenE.txt";
#t $TRINITY_HOME/Trinity --seqType fq --max_memory 120G --left Persea_americana_Hass_POR1c.fq --right Persea_americana_Hass_POR2c.fq --SS_lib_type FR --CPU 10 --output /home/lbmv/EnsambleF/Persea_americana_HassCtrinity --normalize_max_read_cov 50 &>> P_A_H_TRC.txt;
#t echo "Por fin, todo acabó, pero a que precio" >> mensaje.txt;
