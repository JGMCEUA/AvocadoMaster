#!/bin/sh
###Frutos Tratado Correcciones
#python2 /home/lbmv/softwares/git/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 F0_1_R1_T_PO.cor.fq.gz -2 F0_1_R2_T_PO.cor.fq.gz -s F0_1;
#python2 /home/lbmv/softwares/git/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 F0_2_R1_T_PO.cor.fq.gz -2 F0_2_R2_T_PO.cor.fq.gz -s F0_2;
python2 /home/lbmv/softwares/git/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 F0_3_R1_T_PO.cor.fq.gz -2 F0_3_R2_T_PO.cor.fq.gz -s F0_3;
#python2 /home/lbmv/softwares/git/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 F1_1_R1_T_PO.cor.fq.gz -2 F1_1_R2_T_PO.cor.fq.gz -s F1_1;
#python2 /home/lbmv/softwares/git/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 F1_2_R1_T_PO.cor.fq.gz -2 F1_2_R2_T_PO.cor.fq.gz -s F1_2;
python2 /home/lbmv/softwares/git/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 F1_3_R1_T_PO.cor.fq.gz -2 F1_3_R2_T_PO.cor.fq.gz -s F1_3;
#python2 /home/lbmv/softwares/git/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 F24_1_R1_T_PO.cor.fq.gz -2 F24_1_R2_T_PO.cor.fq.gz -s F24_1;
#python2 /home/lbmv/softwares/git/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 F24_2_R1_T_PO.cor.fq.gz -2 F24_2_R2_T_PO.cor.fq.gz -s F24_2;
#python2 /home/lbmv/softwares/git/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 F24_3_R1_T_PO.cor.fq.gz -2 F24_3_R2_T_PO.cor.fq.gz -s F24_3;
#python2 /home/lbmv/softwares/git/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 F72_1_R1_T_PO.cor.fq.gz -2 F72_1_R2_T_PO.cor.fq.gz -s F72_1;
python2 /home/lbmv/softwares/git/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 F72_2_R1_T_PO.cor.fq.gz -2 F72_2_R2_T_PO.cor.fq.gz -s F72_2;
python2 /home/lbmv/softwares/git/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 F72_3_R1_T_PO.cor.fq.gz -2 F72_3_R2_T_PO.cor.fq.gz -s F72_3;
var=$(date); echo "Comando de Correción terminado en FH [ $var]">> "Aviso.txt"
nohup bash "F_zip.bash";
