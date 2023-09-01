#!/bin/sh
var=$(date); echo "Inició [ $var]" >>resumen.txt;
hisat2 --very-sensitive --phred33 --dta -p 8 -x /home/lbmv/JGMC_copias/rRNA_ref/LSU/LSU_DNA_Silva --reorder --un-gz FH0_1_UP_fail.fq.gz --al FH0_1_UP_ali.fq.gz --un-conc-gz FH0_1_PO_fail.fq.gz --al-conc-gz FH0_1_PO_ali.fq.gz --summary-file FH0_1_res.txt --met-file FH0_1_Hisat_metrics.txt --fr -1 unfixrm_F0_1_R1_T_PO.cor.fq.gz -2 unfixrm_F0_1_R2_T_PO.cor.fq.gz &>> resumen.txt
var=$(date); echo "Acabó [ $var]" >> resumen.txt;
