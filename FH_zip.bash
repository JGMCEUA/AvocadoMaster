#!/bin/sh
###Frutos Hass Comprimir
var=$(date); echo "Comando de Compresión inicia en FH [ $var]" >> "Aviso.txt";
gzip -9 unfixrm_FH0_1_R1_T_PO.cor.fq |gzip -9 unfixrm_FH0_1_R2_T_PO.cor.fq;
gzip -9 unfixrm_FH0_2_R1_T_PO.cor.fq |gzip -9 unfixrm_FH0_2_R2_T_PO.cor.fq;
gzip -9 unfixrm_FH0_3_R1_T_PO.cor.fq |gzip -9 unfixrm_FH0_3_R2_T_PO.cor.fq;
gzip -9 unfixrm_FH1_1_R1_T_PO.cor.fq |gzip -9 unfixrm_FH1_1_R2_T_PO.cor.fq;
gzip -9 unfixrm_FH1_2_R1_T_PO.cor.fq |gzip -9 unfixrm_FH1_2_R2_T_PO.cor.fq;
gzip -9 unfixrm_FH1_3_R1_T_PO.cor.fq |gzip -9 unfixrm_FH1_3_R2_T_PO.cor.fq;
gzip -9 unfixrm_FH24_1_R1_T_PO.cor.fq |gzip -9 unfixrm_FH24_1_R2_T_PO.cor.fq;
gzip -9 unfixrm_FH24_2_R1_T_PO.cor.fq |gzip -9 unfixrm_FH24_2_R2_T_PO.cor.fq;
gzip -9 unfixrm_FH24_3_R1_T_PO.cor.fq |gzip -9 unfixrm_FH24_3_R2_T_PO.cor.fq;
gzip -9 unfixrm_FH72_1_R1_T_PO.cor.fq |gzip -9 unfixrm_FH72_1_R2_T_PO.cor.fq;
gzip -9 unfixrm_FH72_2_R1_T_PO.cor.fq |gzip -9 unfixrm_FH72_2_R2_T_PO.cor.fq;
gzip -9 unfixrm_FH72_3_R1_T_PO.cor.fq |gzip -9 unfixrm_FH72_3_R2_T_PO.cor.fq;
var=$(date); echo "Comando de Compresión terminado en FH [ $var]" >> "Aviso.txt"
#cd /home/lbmv/JGMC_copias/FrutoTratado/PO
#nohup bash F_zip.sh
