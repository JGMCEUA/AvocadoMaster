#!/bin/sh
###Frutos Hass Comprimir 
var=$(date); echo "Comando de Compresión empieza en F [ $var]"var=$(date);
gzip -9 unfixrm_F0_1_R1_T_PO.cor.fq |gzip -9 unfixrm_F0_1_R2_T_PO.cor.fq;
gzip -9 unfixrm_F0_2_R1_T_PO.cor.fq |gzip -9 unfixrm_F0_2_R2_T_PO.cor.fq;
gzip -9 unfixrm_F0_3_R1_T_PO.cor.fq |gzip -9 unfixrm_F0_3_R2_T_PO.cor.fq;
gzip -9 unfixrm_F1_1_R1_T_PO.cor.fq |gzip -9 unfixrm_F1_1_R2_T_PO.cor.fq;
gzip -9 unfixrm_F1_2_R1_T_PO.cor.fq |gzip -9 unfixrm_F1_2_R2_T_PO.cor.fq;
gzip -9 unfixrm_F1_3_R1_T_PO.cor.fq |gzip -9 unfixrm_F1_3_R2_T_PO.cor.fq;
gzip -9 unfixrm_F24_1_R1_T_PO.cor.fq |gzip -9 unfixrm_F24_1_R2_T_PO.cor.fq;
gzip -9 unfixrm_F24_2_R1_T_PO.cor.fq |gzip -9 unfixrm_F24_2_R2_T_PO.cor.fq;
gzip -9 unfixrm_F24_3_R1_T_PO.cor.fq |gzip -9 unfixrm_F24_3_R2_T_PO.cor.fq;
gzip -9 unfixrm_F72_1_R1_T_PO.cor.fq |gzip -9 unfixrm_F72_1_R2_T_PO.cor.fq;
gzip -9 unfixrm_F72_2_R1_T_PO.cor.fq |gzip -9 unfixrm_F72_2_R2_T_PO.cor.fq;
gzip -9 unfixrm_F72_3_R1_T_PO.cor.fq |gzip -9 unfixrm_F72_3_R2_T_PO.cor.fq;
var=$(date); echo "Comando de Compresión terminado en F [ $var]"
cd /home/lbmv/JGMC_copias/FrutoHass
nohup bash "FH_Correciones.sh";
