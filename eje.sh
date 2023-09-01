#!/bin/sh
var=$(date); echo "Inició [ $var]" >> resumen.txt; 
hisat2-build /home/lbmv/JGMC_copias/rRNA_ref/LSU/DNA_SILVA_138.1_LSUParc_tax_silva.fasta LSU_DNA_Silva &> resumen.txt; 
var=$(date); echo "Acabó [ $var]" >> resumen.txt
