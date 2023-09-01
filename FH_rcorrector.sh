#!/bin/sh
###Frutos Hass
perl /home/lbmv/softwares/rcorrector/run_rcorrector.pl -1 FH0_1_R1_T_PO.fastq.gz -2 FH0_1_R2_T_PO.fastq.gz -t 12 &> "FH0_1_rcorrrector.log"
perl /home/lbmv/softwares/rcorrector/run_rcorrector.pl -1 FH0_2_R1_T_PO.fastq.gz -2 FH0_2_R2_T_PO.fastq.gz -t 12 &> "FH0_2_rcorrrector.log"
perl /home/lbmv/softwares/rcorrector/run_rcorrector.pl -1 FH0_3_R1_T_PO.fastq.gz -2 FH0_3_R2_T_PO.fastq.gz -t 12 &> "FH0_3_rcorrrector.log"
perl /home/lbmv/softwares/rcorrector/run_rcorrector.pl -1 FH1_1_R1_T_PO.fastq.gz -2 FH1_1_R2_T_PO.fastq.gz -t 12 &> "FH1_1_rcorrrector.log"
perl /home/lbmv/softwares/rcorrector/run_rcorrector.pl -1 FH1_2_R1_T_PO.fastq.gz -2 FH1_2_R2_T_PO.fastq.gz -t 12 &> "FH1_2_rcorrrector.log"
perl /home/lbmv/softwares/rcorrector/run_rcorrector.pl -1 FH1_3_R1_T_PO.fastq.gz -2 FH1_3_R2_T_PO.fastq.gz -t 12 &> "FH1_3_rcorrrector.log"
perl /home/lbmv/softwares/rcorrector/run_rcorrector.pl -1 FH24_1_R1_T_PO.fastq.gz -2 FH24_1_R2_T_PO.fastq.gz -t 12 &> "FH24_1_rcorrrector.log"
perl /home/lbmv/softwares/rcorrector/run_rcorrector.pl -1 FH24_2_R1_T_PO.fastq.gz -2 FH24_2_R2_T_PO.fastq.gz -t 12 &> "FH24_2_rcorrrector.log"
perl /home/lbmv/softwares/rcorrector/run_rcorrector.pl -1 FH24_3_R1_T_PO.fastq.gz -2 FH24_3_R2_T_PO.fastq.gz -t 12 &> "FH24_3_rcorrrector.log"
perl /home/lbmv/softwares/rcorrector/run_rcorrector.pl -1 FH72_1_R1_T_PO.fastq.gz -2 FH72_1_R2_T_PO.fastq.gz -t 12 &> "FH72_1_rcorrrector.log"
perl /home/lbmv/softwares/rcorrector/run_rcorrector.pl 30 -1 FH72_2_R1_T_PO.fastq.gz -2 FH72_2_R2_T_PO.fastq.gz -t 12 &> "FH72_2_rcorrrector.log"
perl /home/lbmv/softwares/rcorrector/run_rcorrector.pl 30 -1 FH72_3_R1_T_PO.fastq.gz -2 FH72_3_R2_T_PO.fastq.gz -t 12 &> "FH72_3_rcorrrector.log"
echo "Acabé por aquí" >> "Aviso.txt";
#cd /home/lbmv/JGMC_copias/FrutoTratado/PO
#bash F_rcorrector.sh
nohup bash "FH_Correcciones.sh";
