#!/bin/sh
sudo apt-get update
sudo apt-get upgrade
###F0_1
stringtie F0_1_PO_QL_Hass.bam -o ./F0_1_PO_QL_Hass_Stringtie/F0_1_PO_QL_Hass_Stringtie.gtf -p 30 --rf -l F_QL_Hass_JGMC -A gene_abudant1.txt 
###F0_2
stringtie F0_2_PO_QL_Hass.bam -o ./F0_2_PO_QL_Hass_Stringtie/F0_2_PO_QL_Hass_Stringtie.gtf -p 30 --rf -l F_QL_Hass_JGMC -A gene_abudant2.txt 
###F0_3
stringtie F0_3_PO_QL_Hass.bam -o ./F0_3_PO_QL_Hass_Stringtie/F0_3_PO_QL_Hass_Stringtie.gtf -p 30 --rf -l F_QL_Hass_JGMC -A gene_abudant3.txt 

###F1_1
stringtie F1_1_PO_QL_Hass.bam -o ./F1_1_PO_QL_Hass_Stringtie/F1_1_PO_QL_Hass_Stringtie.gtf -p 30 --rf -l F_QL_Hass_JGMC -A gene_abudant4.txt 
###F1_2
stringtie F1_2_PO_QL_Hass.bam -o ./F1_2_PO_QL_Hass_Stringtie/F1_2_PO_QL_Hass_Stringtie.gtf -p 30 --rf -l F_QL_Hass_JGMC -A gene_abudant5.txt 
###F1_3
stringtie F1_3_PO_QL_Hass.bam -o ./F1_3_PO_QL_Hass_Stringtie/F1_3_PO_QL_Hass_Stringtie.gtf -p 30 --rf -l F_QL_Hass_JGMC -A gene_abudant6.txt 

###F24_1
stringtie F24_1_PO_QL_Hass.bam -o ./F24_1_PO_QL_Hass_Stringtie/F24_1_PO_QL_Hass_Stringtie.gtf -p 30 --rf -l F_QL_Hass_JGMC -A gene_abudant7.txt 
###F24_2
stringtie F24_2_PO_QL_Hass.bam -o ./F24_2_PO_QL_Hass_Stringtie/F24_2_PO_QL_Hass_Stringtie.gtf -p 30 --rf -l F_QL_Hass_JGMC -A gene_abudant8.txt 
###F24_3
stringtie F24_3_PO_QL_Hass.bam -o ./F24_3_PO_QL_Hass_Stringtie/F24_3_PO_QL_Hass_Stringtie.gtf -p 30 --rf -l F_QL_Hass_JGMC -A gene_abudant9.txt 

###F72_1
stringtie F72_1_PO_QL_Hass.bam -o ./F72_1_PO_QL_Hass_Stringtie/F72_1_PO_QL_Hass_Stringtie.gtf -p 30 --rf -l F_QL_Hass_JGMC -A gene_abudanta1.txt 
###F72_2
stringtie F72_2_PO_QL_Hass.bam -o ./F72_2_PO_QL_Hass_Stringtie/F72_2_PO_QL_Hass_Stringtie.gtf -p 30 --rf -l F_QL_Hass_JGMC -A gene_abudanta2.txt 
###F72_3
stringtie F72_3_PO_QL_Hass.bam -o ./F72_3_PO_QL_Hass_Stringtie/F72_3_PO_QL_Hass_Stringtie.gtf -p 30 --rf -l F_QL_Hass_JGMC -A gene_abudanta3.txt

sudo shutdown now
