#!/bin/sh


 ### Mover carpeta FH0_1 a carpeta de trabaj
 cp /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH0_1_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH0_1_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
 ###Ejecutar alineamiento FH0_1
 hisat2 -p 8 --dta -x PAHCT -1 FH0_1_R1_T_PO.fastq.gz -2 FH0_1_R2_T_PO.fastq.gz -S FH0_1_PAHCT.sam &> FH0_1_PAHCT.txt;
 ###Convertir a archivo Bam FH0_1
 samtools view -bS FH0_1_PAHCT.sam > FH0_1_PAHCT.bam;
 samtools sort FH0_1_PAHCT.bam > FH0_1_PAHCT.sorted.bam;
 samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta FH0_1_PAHCT.sorted.bam | bcftools view -bvcg - > FH0_1_PAHCT.raw.bcf;
 ###Borrar los archivos inecesarios
 cp FH0_1_PAHCT.sorted.bam FH0_1_PAHCT.raw.bcf FH0_1_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
 rm FH*_[123]_R[12]_T_PO.fastq.gz FH0_1_R2_T_PO.fastq.gz FH0_1_PAHCT.sam FH0_1_PAHCT.sorted.bam FH0_1_PAHCT.raw.bcf FH0_1_PAHCT.txt

 ### Mover carpeta FH0_2 a carpeta de trabaj
 cp /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH0_2_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH0_2_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
 ###Ejecutar alineamiento FH0_2
 hisat2 -p 8 --dta -x PAHCT -1 FH0_2_R1_T_PO.fastq.gz -2 FH0_2_R2_T_PO.fastq.gz -S FH0_2_PAHCT.sam &> FH0_2_PAHCT.txt;
 ###Convertir a archivo Bam FH0_2
 samtools view -bS FH0_2_PAHCT.sam > FH0_2_PAHCT.bam;
 samtools sort FH0_2_PAHCT.bam > FH0_2_PAHCT.sorted.bam;
 samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta FH0_2_PAHCT.sorted.bam | bcftools view -bvcg - > FH0_2_PAHCT.raw.bcf;
 ###Borrar los archivos inecesarios
 cp FH0_2_PAHCT.sorted.bam FH0_2_PAHCT.raw.bcf FH0_2_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
 rm FH*_[123]_R[12]_T_PO.fastq.gz FH0_2_R2_T_PO.fastq.gz FH0_2_PAHCT.sam FH0_2_PAHCT.sorted.bam FH0_2_PAHCT.raw.bcf FH0_2_PAHCT.txt

 ### Mover carpeta FH0_3 a carpeta de trabaj
 cp /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH0_3_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH0_3_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
 ###Ejecutar alineamiento FH0_3
 hisat2 -p 8 --dta -x PAHCT -1 FH0_3_R1_T_PO.fastq.gz -2 FH0_3_R2_T_PO.fastq.gz -S FH0_3_PAHCT.sam &> FH0_3_PAHCT.txt;
 ###Convertir a archivo Bam FH0_3
 samtools view -bS FH0_3_PAHCT.sam > FH0_3_PAHCT.bam;
 samtools sort FH0_3_PAHCT.bam > FH0_3_PAHCT.sorted.bam;
 samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta FH0_3_PAHCT.sorted.bam | bcftools view -bvcg - > FH0_3_PAHCT.raw.bcf;
 ###Borrar los archivos inecesarios
 cp FH0_3_PAHCT.sorted.bam FH0_3_PAHCT.raw.bcf FH0_3_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
 rm FH*_[123]_R[12]_T_PO.fastq.gz FH0_3_R2_T_PO.fastq.gz FH0_3_PAHCT.sam FH0_3_PAHCT.sorted.bam FH0_3_PAHCT.raw.bcf FH0_3_PAHCT.txt


 ### Mover carpeta FH1_1 a carpeta de trabaj
 cp /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH1_1_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH1_1_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
 ###Ejecutar alineamiento FH1_1
 hisat2 -p 8 --dta -x PAHCT -1 FH1_1_R1_T_PO.fastq.gz -2 FH1_1_R2_T_PO.fastq.gz -S FH1_1_PAHCT.sam &> FH1_1_PAHCT.txt;
 ###Convertir a archivo Bam FH1_1
 samtools view -bS FH1_1_PAHCT.sam > FH1_1_PAHCT.bam;
 samtools sort FH1_1_PAHCT.bam > FH1_1_PAHCT.sorted.bam;
samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta FH1_1_PAHCT.sorted.bam | bcftools view -bvcg - > FH1_1_PAHCT.raw.bcf;
###Borrar los archivos inecesarios
cp FH1_1_PAHCT.sorted.bam FH1_1_PAHCT.raw.bcf FH1_1_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
rm FH*_[123]_R[12]_T_PO.fastq.gz FH1_1_R2_T_PO.fastq.gz FH1_1_PAHCT.sam FH1_1_PAHCT.sorted.bam FH1_1_PAHCT.raw.bcf FH1_1_PAHCT.txt

### Mover carpeta FH1_2 a carpeta de trabaj
cp /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH1_2_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH1_2_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
###Ejecutar alineamiento FH1_2
hisat2 -p 8 --dta -x PAHCT -1 FH1_2_R1_T_PO.fastq.gz -2 FH1_2_R2_T_PO.fastq.gz -S FH1_2_PAHCT.sam &> FH1_2_PAHCT.txt;
###Convertir a archivo Bam FH1_2
samtools view -bS FH1_2_PAHCT.sam > FH1_2_PAHCT.bam;
samtools sort FH1_2_PAHCT.bam > FH1_2_PAHCT.sorted.bam;
samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta FH1_2_PAHCT.sorted.bam | bcftools view -bvcg - > FH1_2_PAHCT.raw.bcf;
###Borrar los archivos inecesarios
cp FH1_2_PAHCT.sorted.bam FH1_2_PAHCT.raw.bcf FH1_2_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
rm FH*_[123]_R[12]_T_PO.fastq.gz FH1_2_R2_T_PO.fastq.gz FH1_2_PAHCT.sam FH1_2_PAHCT.sorted.bam FH1_2_PAHCT.raw.bcf FH1_2_PAHCT.txt

### Mover carpeta FH1_3 a carpeta de trabaj
cp /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH1_3_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH1_3_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
###Ejecutar alineamiento FH1_3
hisat2 -p 8 --dta -x PAHCT -1 FH1_3_R1_T_PO.fastq.gz -2 FH1_3_R2_T_PO.fastq.gz -S FH1_3_PAHCT.sam &> FH1_3_PAHCT.txt;
###Convertir a archivo Bam FH1_3
samtools view -bS FH1_3_PAHCT.sam > FH1_3_PAHCT.bam;
samtools sort FH1_3_PAHCT.bam > FH1_3_PAHCT.sorted.bam;
samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta FH1_3_PAHCT.sorted.bam | bcftools view -bvcg - > FH1_3_PAHCT.raw.bcf;
###Borrar los archivos inecesarios
cp FH1_3_PAHCT.sorted.bam FH1_3_PAHCT.raw.bcf FH1_3_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
rm FH*_[123]_R[12]_T_PO.fastq.gz FH1_3_R2_T_PO.fastq.gz FH1_3_PAHCT.sam FH1_3_PAHCT.sorted.bam FH1_3_PAHCT.raw.bcf FH1_3_PAHCT.txt


### Mover carpeta FH24_1 a carpeta de trabaj
cp /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH24_1_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH24_1_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
###Ejecutar alineamiento FH24_1
hisat2 -p 8 --dta -x PAHCT -1 FH24_1_R1_T_PO.fastq.gz -2 FH24_1_R2_T_PO.fastq.gz -S FH24_1_PAHCT.sam &> FH24_1_PAHCT.txt;
###Convertir a archivo Bam FH24_1
samtools view -bS FH24_1_PAHCT.sam > FH24_1_PAHCT.bam;
samtools sort FH24_1_PAHCT.bam > FH24_1_PAHCT.sorted.bam;
samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta FH24_1_PAHCT.sorted.bam | bcftools view -bvcg - > FH24_1_PAHCT.raw.bcf;
###Borrar los archivos inecesarios
cp FH24_1_PAHCT.sorted.bam FH24_1_PAHCT.raw.bcf FH24_1_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
rm FH*_[123]_R[12]_T_PO.fastq.gz FH24_1_R2_T_PO.fastq.gz FH24_1_PAHCT.sam FH24_1_PAHCT.sorted.bam FH24_1_PAHCT.raw.bcf FH24_1_PAHCT.txt

### Mover carpeta FH24_2 a carpeta de trabaj
cp /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH24_2_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH24_2_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
###Ejecutar alineamiento FH24_2
hisat2 -p 8 --dta -x PAHCT -1 FH24_2_R1_T_PO.fastq.gz -2 FH24_2_R2_T_PO.fastq.gz -S FH24_2_PAHCT.sam &> FH24_2_PAHCT.txt;
###Convertir a archivo Bam FH24_2
samtools view -bS FH24_2_PAHCT.sam > FH24_2_PAHCT.bam;
samtools sort FH24_2_PAHCT.bam > FH24_2_PAHCT.sorted.bam;
samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta FH24_2_PAHCT.sorted.bam | bcftools view -bvcg - > FH24_2_PAHCT.raw.bcf;
###Borrar los archivos inecesarios
cp FH24_2_PAHCT.sorted.bam FH24_2_PAHCT.raw.bcf FH24_2_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
rm FH*_[123]_R[12]_T_PO.fastq.gz FH24_2_R2_T_PO.fastq.gz FH24_2_PAHCT.sam FH24_2_PAHCT.sorted.bam FH24_2_PAHCT.raw.bcf FH24_2_PAHCT.txt

### Mover carpeta FH24_3 a carpeta de trabaj
cp /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH24_3_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH24_3_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
###Ejecutar alineamiento FH24_3
hisat2 -p 8 --dta -x PAHCT -1 FH24_3_R1_T_PO.fastq.gz -2 FH24_3_R2_T_PO.fastq.gz -S FH24_3_PAHCT.sam &> FH24_3_PAHCT.txt;
###Convertir a archivo Bam FH24_3
samtools view -bS FH24_3_PAHCT.sam > FH24_3_PAHCT.bam;
samtools sort FH24_3_PAHCT.bam > FH24_3_PAHCT.sorted.bam;
samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta FH24_3_PAHCT.sorted.bam | bcftools view -bvcg - > FH24_3_PAHCT.raw.bcf;
###Borrar los archivos inecesarios
cp FH24_3_PAHCT.sorted.bam FH24_3_PAHCT.raw.bcf FH24_3_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
rm FH*_[123]_R[12]_T_PO.fastq.gz FH24_3_R2_T_PO.fastq.gz FH24_3_PAHCT.sam FH24_3_PAHCT.sorted.bam FH24_3_PAHCT.raw.bcf FH24_3_PAHCT.txt



### Mover carpeta FH72_1 a carpeta de trabaj
cp /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH72_1_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH72_1_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
###Ejecutar alineamiento FH72_1
hisat2 -p 8 --dta -x PAHCT -1 FH72_1_R1_T_PO.fastq.gz -2 FH72_1_R2_T_PO.fastq.gz -S FH72_1_PAHCT.sam &> FH72_1_PAHCT.txt;
###Convertir a archivo Bam FH72_1
samtools view -bS FH72_1_PAHCT.sam > FH72_1_PAHCT.bam;
samtools sort FH72_1_PAHCT.bam > FH72_1_PAHCT.sorted.bam;
samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta FH72_1_PAHCT.sorted.bam | bcftools view -bvcg - > FH72_1_PAHCT.raw.bcf;
###Borrar los archivos inecesarios
cp FH72_1_PAHCT.sorted.bam FH72_1_PAHCT.raw.bcf FH72_1_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
rm FH*_[123]_R[12]_T_PO.fastq.gz FH72_1_R2_T_PO.fastq.gz FH72_1_PAHCT.sam FH72_1_PAHCT.sorted.bam FH72_1_PAHCT.raw.bcf FH72_1_PAHCT.txt

### Mover carpeta FH72_2 a carpeta de trabaj
cp /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH72_2_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH72_2_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
###Ejecutar alineamiento FH72_2
hisat2 -p 8 --dta -x PAHCT -1 FH72_2_R1_T_PO.fastq.gz -2 FH72_2_R2_T_PO.fastq.gz -S FH72_2_PAHCT.sam &> FH72_2_PAHCT.txt;
###Convertir a archivo Bam FH72_2
samtools view -bS FH72_2_PAHCT.sam > FH72_2_PAHCT.bam;
samtools sort FH72_2_PAHCT.bam > FH72_2_PAHCT.sorted.bam;
samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta FH72_2_PAHCT.sorted.bam | bcftools view -bvcg - > FH72_2_PAHCT.raw.bcf;
###Borrar los archivos inecesarios
cp FH72_2_PAHCT.sorted.bam FH72_2_PAHCT.raw.bcf FH72_2_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
rm FH*_[123]_R[12]_T_PO.fastq.gz FH72_2_R2_T_PO.fastq.gz FH72_2_PAHCT.sam FH72_2_PAHCT.sorted.bam FH72_2_PAHCT.raw.bcf FH72_2_PAHCT.txt

### Mover carpeta FH72_3 a carpeta de trabaj
cp /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH72_3_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoHass/PO/FH72_3_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
###Ejecutar alineamiento FH72_3
hisat2 -p 8 --dta -x PAHCT -1 FH72_3_R1_T_PO.fastq.gz -2 FH72_3_R2_T_PO.fastq.gz -S FH72_3_PAHCT.sam &> FH72_3_PAHCT.txt;
###Convertir a archivo Bam FH72_3
samtools view -bS FH72_3_PAHCT.sam > FH72_3_PAHCT.bam;
samtools sort FH72_3_PAHCT.bam > FH72_3_PAHCT.sorted.bam;
samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta FH72_3_PAHCT.sorted.bam | bcftools view -bvcg - > FH72_3_PAHCT.raw.bcf;
###Borrar los archivos inecesarios
cp FH72_3_PAHCT.sorted.bam FH72_3_PAHCT.raw.bcf FH72_3_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
rm FH*_[123]_R[12]_T_PO.fastq.gz FH72_3_R2_T_PO.fastq.gz FH72_3_PAHCT.sam FH72_3_PAHCT.sorted.bam FH72_3_PAHCT.raw.bcf FH72_3_PAHCT.txt
