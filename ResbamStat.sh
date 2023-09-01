#!/bin/sh
###Concatenar todos los archivos LOG que nos dió el alinemaiento en uno solo llamado ResBamStat_QL_Hass.txt
###	Jose GUstavo Main Contreras
###	02 08 2021

####
#Codigo
####

#Crear archivo  ResBamStat_QL_Hass.txt
#F0
echo "###	F0_1_PO_QL_Hass" >> ResBamStat_QL_Hass.txt | cat F0_1_PO_QL_Hass_Rseqc.txt >> ResBamStat_QL_Hass.txt
echo "###	F0_2_PO_QL_Hass" >> ResBamStat_QL_Hass.txt | cat F0_2_PO_QL_Hass_Rseqc.txt >> ResBamStat_QL_Hass.txt
echo "###	F0_3_PO_QL_Hass" >> ResBamStat_QL_Hass.txt | cat F0_3_PO_QL_Hass_Rseqc.txt >> ResBamStat_QL_Hass.txt

#F1
echo "###	F1_1_PO_QL_Hass" >> ResBamStat_QL_Hass.txt | cat F1_1_PO_QL_Hass_Rseqc.txt >> ResBamStat_QL_Hass.txt
echo "###	F1_2_PO_QL_Hass" >> ResBamStat_QL_Hass.txt | cat F1_2_PO_QL_Hass_Rseqc.txt >> ResBamStat_QL_Hass.txt
echo "###	F1_3_PO_QL_Hass" >> ResBamStat_QL_Hass.txt | cat F1_3_PO_QL_Hass_Rseqc.txt >> ResBamStat_QL_Hass.txt

#F24
echo "###	F24_1_PO_QL_Hass" >> ResBamStat_QL_Hass.txt | cat F24_1_PO_QL_Hass_Rseqc.txt >> ResBamStat_QL_Hass.txt
echo "###	F24_2_PO_QL_Hass" >> ResBamStat_QL_Hass.txt | cat F24_2_PO_QL_Hass_Rseqc.txt >> ResBamStat_QL_Hass.txt
echo "###	F24_3_PO_QL_Hass" >> ResBamStat_QL_Hass.txt | cat F24_3_PO_QL_Hass_Rseqc.txt >> ResBamStat_QL_Hass.txt

#F72
echo "###	F72_1_PO_QL_Hass" >> ResBamStat_QL_Hass.txt | cat F72_1_PO_QL_Hass_Rseqc.txt >> ResBamStat_QL_Hass.txt
echo "###	F72_2_PO_QL_Hass" >> ResBamStat_QL_Hass.txt | cat F72_2_PO_QL_Hass_Rseqc.txt >> ResBamStat_QL_Hass.txt
echo "###	F72_3_PO_QL_Hass" >> ResBamStat_QL_Hass.txt | cat F72_3_PO_QL_Hass_Rseqc.txt >> ResBamStat_QL_Hass.txt
