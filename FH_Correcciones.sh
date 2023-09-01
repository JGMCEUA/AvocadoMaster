#!/bin/sh
###Frutos Hass Correcciones
var=$(date); echo "Comando de Correción empieza en FH [ $var]">> 'Aviso.txt'$'\r'
python2 /home/lbmv/softwares/git/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 FH0_1_R1_T_PO.cor.fq.gz -2 FH0_1_R2_T_PO.cor.fq.gz -s FH0_1;
python2 /home/lbmv/softwares/git/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 FH0_2_R1_T_PO.cor.fq.gz -2 FH0_2_R2_T_PO.cor.fq.gz -s FH0_2;
python2 /home/lbmv/softwares/git/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 FH0_3_R1_T_PO.cor.fq.gz -2 FH0_3_R2_T_PO.cor.fq.gz -s FH0_3;
python2 /home/lbmv/softwares/git/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 FH1_1_R1_T_PO.cor.fq.gz -2 FH1_1_R2_T_PO.cor.fq.gz -s FH1_1;
python2 /home/lbmv/softwares/git/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 FH1_2_R1_T_PO.cor.fq.gz -2 FH1_2_R2_T_PO.cor.fq.gz -s FH1_2;
python2 /home/lbmv/softwares/git/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 FH1_3_R1_T_PO.cor.fq.gz -2 FH1_3_R2_T_PO.cor.fq.gz -s FH1_3;
python2 /home/lbmv/softwares/git/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 FH24_1_R1_T_PO.cor.fq.gz -2 FH24_1_R2_T_PO.cor.fq.gz -s FH24_1;
python2 /home/lbmv/softwares/git/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 FH24_2_R1_T_PO.cor.fq.gz -2 FH24_2_R2_T_PO.cor.fq.gz -s FH24_2;
python2 /home/lbmv/softwares/git/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 FH24_3_R1_T_PO.cor.fq.gz -2 FH24_3_R2_T_PO.cor.fq.gz -s FH24_3;
python2 /home/lbmv/softwares/git/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 FH72_1_R1_T_PO.cor.fq.gz -2 FH72_1_R2_T_PO.cor.fq.gz -s FH72_1;
python2 /home/lbmv/softwares/git/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 FH72_2_R1_T_PO.cor.fq.gz -2 FH72_2_R2_T_PO.cor.fq.gz -s FH72_2;
python2 /home/lbmv/softwares/git/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 FH72_3_R1_T_PO.cor.fq.gz -2 FH72_3_R2_T_PO.cor.fq.gz -s FH72_3;
var=$(date); echo "Comando de Correción terminado en FH [ $var]">> 'Aviso.txt'$'\r'
#cd /home/lbmv/JGMC_copias/FrutoTratado/PO
#bash F_Correcciones
nohup bash "FH_zip.bash";
