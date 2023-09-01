#!/bin/sh

#t  ### Mover carpeta F0_1 a carpeta de trabaj
#t #t cp /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F0_1_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F0_1_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
#t  ###Ejecutar alineamiento F0_1
#t  hisat2 -p 8 --dta -x PAHCT -1 F0_1_R1_T_PO.fastq.gz -2 F0_1_R2_T_PO.fastq.gz -S F0_1_PAHCT.sam &> F0_1_PAHCT.txt;
#t  ###Convertir a archivo Bam F0_1
#t  samtools view -bS F0_1_PAHCT.sam > F0_1_PAHCT.bam;
#t  samtools sort F0_1_PAHCT.bam > F0_1_PAHCT.sorted.bam;
#t  samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta F0_1_PAHCT.sorted.bam | bcftools view -bvcg - > F0_1_PAHCT.raw.bcf;
#t  ###Borrar los archivos inecesarios
#t  cp F0_1_PAHCT.sorted.bam F0_1_PAHCT.raw.bcf F0_1_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
#t  rm F0_[123]_R[12]_T_PO.fastq.gz F0_1_R2_T_PO.fastq.gz F0_1_PAHCT.sam F0_1_PAHCT.sorted.bam F0_1_PAHCT.raw.bcf F0_1_PAHCT.txt
#t  
 ### Mover carpeta F24_3 a carpeta de trabaj
#t cp /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F24_3_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F24_3_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
###Ejecutar alineamiento F24_3
hisat2 -p 8 --dta -x PAHCT -1 F24_3_R1_T_PO.fastq.gz -2 F24_3_R2_T_PO.fastq.gz -S F24_3_PAHCT.sam &> F24_3_PAHCT.txt;
###Convertir a archivo Bam F24_3
samtools view -bS F24_3_PAHCT.sam > F24_3_PAHCT.bam;
samtools sort F24_3_PAHCT.bam > F24_3_PAHCT.sorted.bam;
samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta F24_3_PAHCT.sorted.bam | bcftools view -bvcg - > F24_3_PAHCT.raw.bcf;
###Borrar los archivos inecesarios
cp F24_3_PAHCT.sorted.bam F24_3_PAHCT.raw.bcf F24_3_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
rm F24_3_R[12]_T_PO.fastq.gz F24_3_PAHCT.sam F24_3_PAHCT.sorted.bam F24_3_PAHCT.raw.bcf F24_3_PAHCT.txt

#t  ### Mover carpeta FH1_1 a carpeta de trabaj
#t  #cp /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH1_1_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH1_1_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
 ###Ejecutar alineamiento FH1_1
 hisat2 -p 8 --dta -x PAHCT -1 FH1_1_R1_T_PO.fastq.gz -2 FH1_1_R2_T_PO.fastq.gz -S FH1_1_PAHCT.sam &> FH1_1_PAHCT.txt;
 ###Convertir a archivo Bam FH1_1
 samtools view -bS FH1_1_PAHCT.sam > FH1_1_PAHCT.bam;
 samtools sort FH1_1_PAHCT.bam > FH1_1_PAHCT.sorted.bam;
samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta FH1_1_PAHCT.sorted.bam | bcftools view -bvcg - > FH1_1_PAHCT.raw.bcf;
###Borrar los archivos inecesarios
cp FH1_1_PAHCT.sorted.bam FH1_1_PAHCT.raw.bcf FH1_1_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
rm FH1_1_R[12]_T_PO.fastq.gz FH1_1_PAHCT.sam FH1_1_PAHCT.sorted.bam FH1_1_PAHCT.raw.bcf FH1_1_PAHCT.txt

#t #t  ### Mover carpeta FH0_1 a carpeta de trabaj
#t #t  #cp /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH0_1_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH0_1_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
 ###Ejecutar alineamiento FH0_1
 hisat2 -p 8 --dta -x PAHCT -1 FH0_1_R1_T_PO.fastq.gz -2 FH0_1_R2_T_PO.fastq.gz -S FH0_1_PAHCT.sam &> FH0_1_PAHCT.txt;
 ###Convertir a archivo Bam FH0_1
 samtools view -bS FH0_1_PAHCT.sam > FH0_1_PAHCT.bam;
 samtools sort FH0_1_PAHCT.bam > FH0_1_PAHCT.sorted.bam;
 samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta FH0_1_PAHCT.sorted.bam | bcftools view -bvcg - > FH0_1_PAHCT.raw.bcf;
 ###Borrar los archivos inecesarios
 cp FH0_1_PAHCT.sorted.bam FH0_1_PAHCT.raw.bcf FH0_1_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
 rm FH0_1_R[12]_T_PO.fastq.gz FH0_1_PAHCT.sam FH0_1_PAHCT.sorted.bam FH0_1_PAHCT.raw.bcf FH0_1_PAHCT.txt

#t  ### Mover carpeta FH0_2 a carpeta de trabaj
#t  #cp /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH0_2_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH0_2_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
 ###Ejecutar alineamiento FH0_2
 hisat2 -p 8 --dta -x PAHCT -1 FH0_2_R1_T_PO.fastq.gz -2 FH0_2_R2_T_PO.fastq.gz -S FH0_2_PAHCT.sam &> FH0_2_PAHCT.txt;
 ###Convertir a archivo Bam FH0_2
 samtools view -bS FH0_2_PAHCT.sam > FH0_2_PAHCT.bam;
 samtools sort FH0_2_PAHCT.bam > FH0_2_PAHCT.sorted.bam;
samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta FH0_2_PAHCT.sorted.bam | bcftools view -bvcg - > FH0_2_PAHCT.raw.bcf;
###Borrar los archivos inecesarios
cp FH0_2_PAHCT.sorted.bam FH0_2_PAHCT.raw.bcf FH0_2_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
rm FH0_2_R[12]_T_PO.fastq.gz FH0_2_PAHCT.sam FH0_2_PAHCT.sorted.bam FH0_2_PAHCT.raw.bcf FH0_2_PAHCT.txt

#t  ### Mover carpeta FH0_3 a carpeta de trabaj
#t  #cp /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH0_3_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH0_3_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
 ###Ejecutar alineamiento FH0_3
 hisat2 -p 8 --dta -x PAHCT -1 FH0_3_R1_T_PO.fastq.gz -2 FH0_3_R2_T_PO.fastq.gz -S FH0_3_PAHCT.sam &> FH0_3_PAHCT.txt;
 ###Convertir a archivo Bam FH0_3
 samtools view -bS FH0_3_PAHCT.sam > FH0_3_PAHCT.bam;
 samtools sort FH0_3_PAHCT.bam > FH0_3_PAHCT.sorted.bam;
samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta FH0_3_PAHCT.sorted.bam | bcftools view -bvcg - > FH0_3_PAHCT.raw.bcf;
###Borrar los archivos inecesarios
cp FH0_3_PAHCT.sorted.bam FH0_3_PAHCT.raw.bcf FH0_3_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
rm FH0_3_R[12]_T_PO.fastq.gz FH0_3_PAHCT.sam FH0_3_PAHCT.sorted.bam FH0_3_PAHCT.raw.bcf FH0_3_PAHCT.txt

