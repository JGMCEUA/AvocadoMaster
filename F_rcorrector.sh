#!/bin/sh
###Frutos Tratado
#perl /home/lbmv/softwares/rcorrector/run_rcorrector.pl -1 F0_1_R1_T_PO.fastq.gz -2 F0_1_R2_T_PO.fastq.gz -t 12 &> "F0_1_rcorrrector.log"
#perl /home/lbmv/softwares/rcorrector/run_rcorrector.pl -1 F0_2_R1_T_PO.fastq.gz -2 F0_2_R2_T_PO.fastq.gz -t 12 &> "F0_2_rcorrrector.log"
perl /home/lbmv/softwares/rcorrector/run_rcorrector.pl -t 12 -1 F0_3_R1_T_PO.fastq.gz -2 F0_3_R2_T_PO.fastq.gz &> "F0_3_rcorrrector.log";
#perl /home/lbmv/softwares/rcorrector/run_rcorrector.pl -1 F1_1_R1_T_PO.fastq.gz -2 F1_1_R2_T_PO.fastq.gz -t 12 &> "F1_1_rcorrrector.log"
#perl /home/lbmv/softwares/rcorrector/run_rcorrector.pl -1 F1_2_R1_T_PO.fastq.gz -2 F1_2_R2_T_PO.fastq.gz -t 12 &> "F1_2_rcorrrector.log"
perl /home/lbmv/softwares/rcorrector/run_rcorrector.pl -t 12 -1 F1_3_R1_T_PO.fastq.gz -2 F1_3_R2_T_PO.fastq.gz &> "F1_3_rcorrrector.log";
#perl /home/lbmv/softwares/rcorrector/run_rcorrector.pl -1 F24_1_R1_T_PO.fastq.gz -2 F24_1_R2_T_PO.fastq.gz -t 12 &> "F24_1_rcorrrector.log"
#perl /home/lbmv/softwares/rcorrector/run_rcorrector.pl -1 F24_2_R1_T_PO.fastq.gz -2 F24_2_R2_T_PO.fastq.gz -t 12 &> "F24_2_rcorrrector.log"
#perl /home/lbmv/softwares/rcorrector/run_rcorrector.pl -1 F24_3_R1_T_PO.fastq.gz -2 F24_3_R2_T_PO.fastq.gz -t 12 &> "F24_3_rcorrrector.log"
#perl /home/lbmv/softwares/rcorrector/run_rcorrector.pl -1 F72_1_R1_T_PO.fastq.gz -2 F72_1_R2_T_PO.fastq.gz -t 12 &> "F72_1_rcorrrector.log"
perl /home/lbmv/softwares/rcorrector/run_rcorrector.pl -t 12 -1 F72_2_R1_T_PO.fastq.gz -2 F72_2_R2_T_PO.fastq.gz &> "F72_2_rcorrrector.log";
#perl /home/lbmv/softwares/rcorrector/run_rcorrector.pl -1 F72_3_R1_T_PO.fastq.gz -2 F72_3_R2_T_PO.fastq.gz -t 12 &> "F72_3_rcorrrector.log"
var=$(date); echo "Hola hoy es [ $var], respeta al peaton" >> "Aviso.txt";
nohup bash "F_Correcciones.sh";
