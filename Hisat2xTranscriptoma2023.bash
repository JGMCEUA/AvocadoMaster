### SCRIPT PARA EL MAPEO DE LAS SECUENCIAS CON EL GENOMA ENSAMBLADO DE FAAGP Y LAS SECUENCIAS DE XOCA
##23/03/2023
#José Gustavo Marín Contreras
##Alieamiento de las secuencias del transcriptoma de Angelica con el genoma de referencia
#####AGUACATE+FAAGP
#### Aguacate+FAAGP a 0 hpt (F0)
hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1  '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoTratado/PO/F0_1_R1_T_PO.fastq.gz' -2 '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoTratado/PO/F0_1_R2_T_PO.fastq.gz' -S F0_1_2023.sam &> F0_1_2023.txt;
samtools sort -@ 30 -o F0_1_2023.bam F0_1_2023.sam;
rm F0_1_2023.sam;
hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1  '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoTratado/PO/F0_2_R1_T_PO.fastq.gz' -2 '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoTratado/PO/F0_2_R2_T_PO.fastq.gz' -S F0_2_2023.sam &> F0_2_2023.txt;
samtools sort -@ 30 -o F0_2_2023.bam F0_2_2023.sam;
rm F0_2_2023.sam;
hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1  '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoTratado/PO/F0_3_R1_T_PO.fastq.gz' -2 '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoTratado/PO/F0_3_R2_T_PO.fastq.gz' -S F0_3_2023.sam &> F0_3_2023.txt;
samtools sort -@ 30 -o F0_3_2023.bam F0_3_2023.sam;
rm F0_3_2023.sam;
#### Aguacate+FAAGP a 24 hpt (F24)
hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1  '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoTratado/PO/F24_1_R1_T_PO.fastq.gz' -2 '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoTratado/PO/F24_1_R2_T_PO.fastq.gz' -S F24_1_2023.sam &> F24_1_2023.txt;
samtools sort -@ 30 -o F24_1_2023.bam F24_1_2023.sam;
rm F24_1_2023.sam;
hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1  '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoTratado/PO/F24_2_R1_T_PO.fastq.gz' -2 '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoTratado/PO/F24_2_R2_T_PO.fastq.gz' -S F24_2_2023.sam &> F24_2_2023.txt;
samtools sort -@ 30 -o F24_2_2023.bam F24_2_2023.sam;
rm F24_2_2023.sam;
hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1  '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoTratado/PO/F24_3_R1_T_PO.fastq.gz' -2 '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoTratado/PO/F24_3_R2_T_PO.fastq.gz' -S F24_3_2023.sam &> F24_3_2023.txt;
samtools sort -@ 30 -o F24_3_2023.bam F24_3_2023.sam;
rm F24_3_2023.sam;
#### Aguacate+FAAGP a 72 hpt (F72)
hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1  '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoTratado/PO/F72_1_R1_T_PO.fastq.gz' -2 '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoTratado/PO/F72_1_R2_T_PO.fastq.gz' -S F72_1_2023.sam &> F72_1_2023.txt;
samtools sort -@ 30 -o F72_1_2023.bam F72_1_2023.sam;
rm F72_1_2023.sam;
hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1  '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoTratado/PO/F72_2_R1_T_PO.fastq.gz' -2 '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoTratado/PO/F72_2_R2_T_PO.fastq.gz' -S F72_2_2023.sam &> F72_2_2023.txt;
samtools sort -@ 30 -o F72_2_2023.bam F72_2_2023.sam;
rm F72_2_2023.sam;
hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1  '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoTratado/PO/F72_3_R1_T_PO.fastq.gz' -2 '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoTratado/PO/F72_3_R2_T_PO.fastq.gz' -S F72_3_2023.sam &> F72_3_2023.txt;
samtools sort -@ 30 -o F72_3_2023.bam F72_3_2023.sam;
rm F72_3_2023.sam;
#####AGUACATE CONTROL
#### Aguacate CONTROL a 0 hpt (FH0)
hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1  '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoHass/PO/FH0_1_R1_T_PO.fastq.gz' -2 '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoHass/PO/FH0_1_R2_T_PO.fastq.gz' -S FH0_1_2023.sam &> FH0_1_2023.txt;
samtools sort -@ 30 -o FH0_1_2023.bam FH0_1_2023.sam;
rm FH0_1_2023.sam;
hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1  '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoHass/PO/FH0_2_R1_T_PO.fastq.gz' -2 '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoHass/PO/FH0_2_R2_T_PO.fastq.gz' -S FH0_2_2023.sam &> FH0_2_2023.txt;
samtools sort -@ 30 -o FH0_2_2023.bam FH0_2_2023.sam;
rm FH0_2_2023.sam;
hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1  '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoHass/PO/FH0_3_R1_T_PO.fastq.gz' -2 '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoHass/PO/FH0_3_R2_T_PO.fastq.gz' -S FH0_3_2023.sam &> FH0_3_2023.txt;
samtools sort -@ 30 -o FH0_3_2023.bam FH0_3_2023.sam;
rm FH0_3_2023.sam;
#### Aguacate CONTROL a 0 hpt (FH24)
hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1  '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoHass/PO/FH24_1_R1_T_PO.fastq.gz' -2 '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoHass/PO/FH24_1_R2_T_PO.fastq.gz' -S FH24_1_2023.sam &> FH24_1_2023.txt;
samtools sort -@ 30 -o FH24_1_2023.bam FH24_1_2023.sam;
rm FH24_1_2023.sam;
hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1  '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoHass/PO/FH24_2_R1_T_PO.fastq.gz' -2 '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoHass/PO/FH24_2_R2_T_PO.fastq.gz' -S FH24_2_2023.sam &> FH24_2_2023.txt;
samtools sort -@ 30 -o FH24_2_2023.bam FH24_2_2023.sam;
rm FH24_2_2023.sam;
hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1  '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoHass/PO/FH24_3_R1_T_PO.fastq.gz' -2 '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoHass/PO/FH24_3_R2_T_PO.fastq.gz' -S FH24_3_2023.sam &> FH24_3_2023.txt;
samtools sort -@ 30 -o FH24_3_2023.bam FH24_3_2023.sam;
rm FH24_3_2023.sam;
#### Aguacate CONTROL a 0 hpt (FH72)
hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1  '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoHass/PO/FH72_1_R1_T_PO.fastq.gz' -2 '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoHass/PO/FH72_1_R2_T_PO.fastq.gz' -S FH72_1_2023.sam &> FH72_1_2023.txt;
samtools sort -@ 30 -o FH72_1_2023.bam FH72_1_2023.sam;
rm FH72_1_2023.sam;
hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1  '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoHass/PO/FH72_2_R1_T_PO.fastq.gz' -2 '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoHass/PO/FH72_2_R2_T_PO.fastq.gz' -S FH72_2_2023.sam &> FH72_2_2023.txt;
samtools sort -@ 30 -o FH72_2_2023.bam FH72_2_2023.sam;
rm FH72_2_2023.sam;
hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1  '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoHass/PO/FH72_3_R1_T_PO.fastq.gz' -2 '/media/alejandra/Guzz_Infor/AgroBio_2TB/Proyecto_Gus_Angie/ETAPA1/RecorteCalidad/MuestrasTrimeadas/FrutoHass/PO/FH72_3_R2_T_PO.fastq.gz' -S FH72_3_2023.sam &> FH72_3_2023.txt;
samtools sort -@ 30 -o FH72_3_2023.bam FH72_3_2023.sam;
rm FH72_3_2023.sam;