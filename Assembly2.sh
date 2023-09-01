#!/bin/sh
###Realizar el ensamblado de los archivos bam del Fruto Control
###F0_1
#cufflinks -p 30 -b /home/alejandra/Escritorio/estancia/referenceG/NCBI_Queensland_Hass/GCA_018408905.1_UQ_Hass_1.0_assembly_report.faa --library-type fr-firststrand -N --total-hits-norm -L F_QL_Hass_JGMC  F0_1_PO_QL_Hass.bam -o F0_1_PO_QL_Hass_LT
###F0_2
cufflinks -p 30 -b /home/alejandra/Escritorio/estancia/referenceG/NCBI_Queensland_Hass/GCA_018408905.1_UQ_Hass_1.0_assembly_report.faa --library-type fr-firststrand -N --total-hits-norm -L F_QL_Hass_JGMC  F0_2_PO_QL_Hass.bam -o F0_2_PO_QL_Hass_LT
###F0_3
cufflinks -p 30 -b /home/alejandra/Escritorio/estancia/referenceG/NCBI_Queensland_Hass/GCA_018408905.1_UQ_Hass_1.0_assembly_report.faa --library-type fr-firststrand -N --total-hits-norm -L F_QL_Hass_JGMC  F0_3_PO_QL_Hass.bam -o F0_3_PO_QL_Hass_LT

###F1_1
cufflinks -p 30 -b /home/alejandra/Escritorio/estancia/referenceG/NCBI_Queensland_Hass/GCA_018408905.1_UQ_Hass_1.0_assembly_report.faa --library-type fr-firststrand -N --total-hits-norm -L F_QL_Hass_JGMC  F1_1_PO_QL_Hass.bam -o F1_1_PO_QL_Hass_LT
###F1_2
cufflinks -p 30 -b /home/alejandra/Escritorio/estancia/referenceG/NCBI_Queensland_Hass/GCA_018408905.1_UQ_Hass_1.0_assembly_report.faa --library-type fr-firststrand -N --total-hits-norm -L F_QL_Hass_JGMC  F1_2_PO_QL_Hass.bam -o F1_2_PO_QL_Hass_LT
###F1_3
cufflinks -p 30 -b /home/alejandra/Escritorio/estancia/referenceG/NCBI_Queensland_Hass/GCA_018408905.1_UQ_Hass_1.0_assembly_report.faa --library-type fr-firststrand -N --total-hits-norm -L F_QL_Hass_JGMC  F1_3_PO_QL_Hass.bam -o F1_3_PO_QL_Hass_LT

###F24_1
cufflinks -p 30 -b /home/alejandra/Escritorio/estancia/referenceG/NCBI_Queensland_Hass/GCA_018408905.1_UQ_Hass_1.0_assembly_report.faa --library-type fr-firststrand -N --total-hits-norm -L F_QL_Hass_JGMC  F24_1_PO_QL_Hass.bam -o F24_1_PO_QL_Hass_LT
###F24_2
cufflinks -p 30 -b /home/alejandra/Escritorio/estancia/referenceG/NCBI_Queensland_Hass/GCA_018408905.1_UQ_Hass_1.0_assembly_report.faa --library-type fr-firststrand -N --total-hits-norm -L F_QL_Hass_JGMC  F24_2_PO_QL_Hass.bam -o F24_2_PO_QL_Hass_LT
###F24_3
cufflinks -p 30 -b /home/alejandra/Escritorio/estancia/referenceG/NCBI_Queensland_Hass/GCA_018408905.1_UQ_Hass_1.0_assembly_report.faa --library-type fr-firststrand -N --total-hits-norm -L F_QL_Hass_JGMC  F24_3_PO_QL_Hass.bam -o F24_3_PO_QL_Hass_LT

###F72_1
cufflinks -p 30 -b /home/alejandra/Escritorio/estancia/referenceG/NCBI_Queensland_Hass/GCA_018408905.1_UQ_Hass_1.0_assembly_report.faa --library-type fr-firststrand -N --total-hits-norm -L F_QL_Hass_JGMC  F72_1_PO_QL_Hass.bam -o F72_1_PO_QL_Hass_LT
###F72_2
cufflinks -p 30 -b /home/alejandra/Escritorio/estancia/referenceG/NCBI_Queensland_Hass/GCA_018408905.1_UQ_Hass_1.0_assembly_report.faa --library-type fr-firststrand -N --total-hits-norm -L F_QL_Hass_JGMC  F72_2_PO_QL_Hass.bam -o F72_2_PO_QL_Hass_LT
###F72_3
cufflinks -p 30 -b /home/alejandra/Escritorio/estancia/referenceG/NCBI_Queensland_Hass/GCA_018408905.1_UQ_Hass_1.0_assembly_report.faa --library-type fr-firststrand -N --total-hits-norm -L F_QL_Hass_JGMC  F72_3_PO_QL_Hass.bam -o F72_3_PO_QL_Hass_LT

