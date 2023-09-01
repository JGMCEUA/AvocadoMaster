#!/bin/sh
#Concatenación de librerias Tratadas
#t cd /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble;
#t var=$(date); echo "Comenzó concatenación de F[ $var]">> "Resumen.txt";
#t cat F0_[123]_PO_fail_f2.R1.fq >> F0_PO_trinity_R1.fq | cat F0_[123]_PO_fail_f2.R2.fq >> F0_PO_trinity_R2.fq;
#t cat F1_[123]_PO_fail_f2.R1.fq >> F1_PO_trinity_R1.fq | cat F1_[123]_PO_fail_f2.R2.fq >> F1_PO_trinity_R2.fq;
#t cat F24_[123]_PO_fail_f2.R1.fq >> F24_PO_trinity_R1.fq | cat F24_[123]_PO_fail_f2.R2.fq >> F24_PO_trinity_R2.fq;
#t cat F72_[123]_PO_fail_f2.R1.fq >> F72_PO_trinity_R1.fq | cat F72_[123]_PO_fail_f2.R2.fq >> F72_PO_trinity_R2.fq;
#t var=$(date); echo "Terminó concatenación de F[ $var]">> "Resumen.txt";
#t #Concatenación de librerias Control
#t cd /home/lbmv/JGMC_copias/FrutoHass/filtro2/filtradosFH/ensamble;
#t var=$(date); echo "Comenzó concatenación de FH[ $var]">> "Resumen.txt";
#t cat FH0_[123]_PO_fail_f2.R1.fq >> FH0_PO_trinity_R1.fq | cat FH0_[123]_PO_fail_f2.R2.fq >> FH0_PO_trinity_R2.fq;
#t cat FH1_[123]_PO_fail_f2.R1.fq >> FH1_PO_trinity_R1.fq | cat FH1_[123]_PO_fail_f2.R2.fq >> FH1_PO_trinity_R2.fq;
#t cat FH24_[123]_PO_fail_f2.R1.fq >> FH24_PO_trinity_R1.fq | cat FH24_[123]_PO_fail_f2.R2.fq >> FH24_PO_trinity_R2.fq;
#t cat FH72_[123]_PO_fail_f2.R1.fq >> FH72_PO_trinity_R1.fq | cat FH72_[123]_PO_fail_f2.R2.fq >> FH72_PO_trinity_R2.fq;
#t var=$(date); echo "Terminó concatenación de FH[ $var]">> "Resumen.txt";
#Ensamblado independiente F
#t cd /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble;
#t var=$(date); echo "Comenzó ensamble de F0[ $var]">> "ResumenE.txt";
#t $TRINITY_HOME/Trinity --seqType fq --max_memory 120G --left F0_PO_trinity_R1.fq --right F0_PO_trinity_R2.fq --SS_lib_type FR --CPU 10 --output /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble/F0trinity --normalize_max_read_cov 100 &>>F0_TR.txt;
#t var=$(date); echo "Teminó ensamble de F0[ $var]">> "ResumenE.txt";
#t var=$(date); echo "Comenzó ensamble de F1[ $var]">> "ResumenE.txt";
#t $TRINITY_HOME/Trinity --seqType fq --max_memory 120G --left F1_PO_trinity_R1.fq --right F1_PO_trinity_R2.fq --SS_lib_type FR --CPU 10 --output /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble/F1trinity --normalize_max_read_cov 100 &>>F1_TR.txt;
#t var=$(date); echo "Teminó ensamble de F1[ $var]">> "ResumenE.txt";
#t var=$(date); echo "Comenzó ensamble de F24[ $var]">> "ResumenE.txt";
#t $TRINITY_HOME/Trinity --seqType fq --max_memory 120G --left F24_PO_trinity_R1.fq --right F24_PO_trinity_R2.fq --SS_lib_type FR --CPU 10 --output /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble/F24trinity --normalize_max_read_cov 100 &>>F24_TR.txt;
#t var=$(date); echo "Teminó ensamble de F24[ $var]">> "ResumenE.txt";
#t var=$(date); echo "Comenzó ensamble de F72[ $var]">> "ResumenE.txt";
#t $TRINITY_HOME/Trinity --seqType fq --max_memory 120G --left F72_PO_trinity_R1.fq --right F72_PO_trinity_R2.fq --SS_lib_type FR --CPU 10 --output /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble/F72trinity --normalize_max_read_cov 100 &>>F72_TR.txt;
#t var=$(date); echo "Teminó ensamble de F72[ $var]">> "ResumenE.txt";
#t #Ensamblado independiente FH
#t cd /home/lbmv/JGMC_copias/FrutoHass/filtro2/filtradosFH/ensamble;
#t var=$(date); echo "Comenzó ensamble de FH0[ $var]">> "ResumenE.txt";
#t $TRINITY_HOME/Trinity --seqType fq --max_memory 120G --left FH0_PO_trinity_R1.fq --right FH0_PO_trinity_R2.fq --SS_lib_type FR --CPU 10 --output /home/lbmv/JGMC_copias/FrutoHass/filtro2/filtradosFH/ensamble/FH0trinity --normalize_max_read_cov 100 &>> FH0_TR.txt;
#t var=$(date); echo "Terminó ensamble de FH0[ $var]">> "ResumenE.txt";
#t var=$(date); echo "Comenzó ensamble de FH1[ $var]">> "ResumenE.txt";
#t $TRINITY_HOME/Trinity --seqType fq --max_memory 120G --left FH1_PO_trinity_R1.fq --right FH1_PO_trinity_R2.fq --SS_lib_type FR --CPU 10 --output /home/lbmv/JGMC_copias/FrutoHass/filtro2/filtradosFH/ensamble/FH1trinity --normalize_max_read_cov 100 &>> FH1_TR.txt;
#t var=$(date); echo "Terminó ensamble de FH1[ $var]">> "ResumenE.txt";
#t var=$(date); echo "Comenzó ensamble de FH24[ $var]">> "ResumenE.txt";
#t $TRINITY_HOME/Trinity --seqType fq --max_memory 120G --left FH24_PO_trinity_R1.fq --right FH24_PO_trinity_R2.fq --SS_lib_type FR --CPU 10 --output /home/lbmv/JGMC_copias/FrutoHass/filtro2/filtradosFH/ensamble/FH24trinity --normalize_max_read_cov 100 &>> FH24_TR.txt;
#t var=$(date); echo "Terminó ensamble de FH24[ $var]">> "ResumenE.txt";
#t var=$(date); echo "Comenzó ensamble de FH72[ $var]">> "ResumenE.txt";
#t $TRINITY_HOME/Trinity --seqType fq --max_memory 120G --left FH72_PO_trinity_R1.fq --right FH72_PO_trinity_R2.fq --SS_lib_type FR --CPU 10 --output /home/lbmv/JGMC_copias/FrutoHass/filtro2/filtradosFH/ensamble/FH72trinity --normalize_max_read_cov 100 &>> FH72_TR.txt;
#t var=$(date); echo "Terminó ensamble de FH72[ $var]">> "ResumenE.txt";
#t ###Emsamblado total de las lecturas F
#t cd /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble;
#t var=$(date); echo "Comenzó concatenación de F Total[ $var]">> "Resumen.txt";
#t cat F0_PO_trinity_R1.fq F1_PO_trinity_R1.fq F24_PO_trinity_R1.fq F72_PO_trinity_R1.fq >> F_PO_trinity_R1.fq;
#t cat F0_PO_trinity_R2.fq F1_PO_trinity_R2.fq F24_PO_trinity_R2.fq F72_PO_trinity_R2.fq >> F_PO_trinity_R2.fq;
#t mv F_PO_trinity_R1.fq F_PO_trinity_R2.fq /home/lbmv/EnsambleF;
#t var=$(date); echo "Terminó concatenación de F Total[ $var]">> "Resumen.txt";
#t ###Emsamblado total de las lecturas FH
#t cd /home/lbmv/JGMC_copias/FrutoHass/filtro2/filtradosFH/ensamble;
#t var=$(date); echo "Comenzó concatenación de F Total[ $var]">> "Resumen.txt";
#t cat FH0_PO_trinity_R1.fq FH1_PO_trinity_R1.fq FH24_PO_trinity_R1.fq FH72_PO_trinity_R1.fq >> FH_PO_trinity_R1.fq;
#t cat FH0_PO_trinity_R2.fq FH1_PO_trinity_R2.fq FH24_PO_trinity_R2.fq FH72_PO_trinity_R2.fq >> FH_PO_trinity_R2.fq;
#t mv FH_PO_trinity_R1.fq FH_PO_trinity_R2.fq /home/lbmv/EnsambleF;
#t var=$(date); echo "Terminó concatenación de F Total[ $var]">> "Resumen.txt";
#t ###Concatenación Final suprema
cd /home/lbmv/EnsambleF;
#t cat FH_PO_trinity_R1.fq F_PO_trinity_R1.fq >> Persea_americana_Hass_POR1.fq;
#t cat FH_PO_trinity_R2.fq F_PO_trinity_R2.fq >> Persea_americana_Hass_POR2.fq;
###Ensamblado supremo
var=$(date); echo "Comenzó ensamble de P.A.Hass[ $var]">> "ResumenE.txt";
$TRINITY_HOME/Trinity --seqType fq --max_memory 120G --left Persea_americana_Hass_POR1.fq --right Persea_americana_Hass_POR2.fq --SS_lib_type FR --CPU 10 --output /home/lbmv/EnsambleF/Persea_americana_Hasstrinity --normalize_max_read_cov 50 &>> P_A_H_TR.txt;
var=$(date); echo "Terminó ensamble de P.A.Hass[ $var]">> "ResumenE.txt";
