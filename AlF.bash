#!/bin/sh


#t  ### Mover carpeta F0_1 a carpeta de trabaj
#t  cp /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F0_1_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F0_1_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
#t  ###Ejecutar alineamiento F0_1
#t  hisat2 -p 8 --dta -x PAHCT -1 F0_1_R1_T_PO.fastq.gz -2 F0_1_R2_T_PO.fastq.gz -S F0_1_PAHCT.sam &> F0_1_PAHCT.txt;
#t  ###Convertir a archivo Bam F0_1
#t  samtools view -bS F0_1_PAHCT.sam > F0_1_PAHCT.bam;
#t  samtools sort F0_1_PAHCT.bam > F0_1_PAHCT.sorted.bam;
#t  samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta F0_1_PAHCT.sorted.bam | bcftools view -bvcg - > F0_1_PAHCT.raw.bcf;
#t  ###Borrar los archivos inecesarios
#t  cp F0_1_PAHCT.sorted.bam F0_1_PAHCT.raw.bcf F0_1_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
#t  rm F*_[123]_R[12]_T_PO.fastq.gz F0_1_R2_T_PO.fastq.gz F0_1_PAHCT.sam F0_1_PAHCT.sorted.bam F0_1_PAHCT.raw.bcf F0_1_PAHCT.txt

 ### Mover carpeta F0_2 a carpeta de trabaj
 cp /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F0_2_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F0_2_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
 ###Ejecutar alineamiento F0_2
 hisat2 -p 8 --dta -x PAHCT -1 F0_2_R1_T_PO.fastq.gz -2 F0_2_R2_T_PO.fastq.gz -S F0_2_PAHCT.sam &> F0_2_PAHCT.txt;
 ###Convertir a archivo Bam F0_2
 samtools view -bS F0_2_PAHCT.sam > F0_2_PAHCT.bam;
 samtools sort F0_2_PAHCT.bam > F0_2_PAHCT.sorted.bam;
 samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta F0_2_PAHCT.sorted.bam | bcftools view -bvcg - > F0_2_PAHCT.raw.bcf;
 ###Borrar los archivos inecesarios
 cp F0_2_PAHCT.sorted.bam F0_2_PAHCT.raw.bcf F0_2_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
 rm F*_[123]_R[12]_T_PO.fastq.gz F0_2_R2_T_PO.fastq.gz F0_2_PAHCT.sam F0_2_PAHCT.sorted.bam F0_2_PAHCT.raw.bcf F0_2_PAHCT.txt

 ### Mover carpeta F0_3 a carpeta de trabaj
 cp /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F0_3_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F0_3_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
 ###Ejecutar alineamiento F0_3
 hisat2 -p 8 --dta -x PAHCT -1 F0_3_R1_T_PO.fastq.gz -2 F0_3_R2_T_PO.fastq.gz -S F0_3_PAHCT.sam &> F0_3_PAHCT.txt;
 ###Convertir a archivo Bam F0_3
 samtools view -bS F0_3_PAHCT.sam > F0_3_PAHCT.bam;
 samtools sort F0_3_PAHCT.bam > F0_3_PAHCT.sorted.bam;
 samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta F0_3_PAHCT.sorted.bam | bcftools view -bvcg - > F0_3_PAHCT.raw.bcf;
 ###Borrar los archivos inecesarios
 cp F0_3_PAHCT.sorted.bam F0_3_PAHCT.raw.bcf F0_3_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
 rm F*_[123]_R[12]_T_PO.fastq.gz F0_3_R2_T_PO.fastq.gz F0_3_PAHCT.sam F0_3_PAHCT.sorted.bam F0_3_PAHCT.raw.bcf F0_3_PAHCT.txt


 ### Mover carpeta F1_1 a carpeta de trabaj
 cp /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F1_1_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F1_1_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
 ###Ejecutar alineamiento F1_1
 hisat2 -p 8 --dta -x PAHCT -1 F1_1_R1_T_PO.fastq.gz -2 F1_1_R2_T_PO.fastq.gz -S F1_1_PAHCT.sam &> F1_1_PAHCT.txt;
 ###Convertir a archivo Bam F1_1
 samtools view -bS F1_1_PAHCT.sam > F1_1_PAHCT.bam;
 samtools sort F1_1_PAHCT.bam > F1_1_PAHCT.sorted.bam;
samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta F1_1_PAHCT.sorted.bam | bcftools view -bvcg - > F1_1_PAHCT.raw.bcf;
###Borrar los archivos inecesarios
cp F1_1_PAHCT.sorted.bam F1_1_PAHCT.raw.bcf F1_1_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
rm F*_[123]_R[12]_T_PO.fastq.gz F1_1_R2_T_PO.fastq.gz F1_1_PAHCT.sam F1_1_PAHCT.sorted.bam F1_1_PAHCT.raw.bcf F1_1_PAHCT.txt

### Mover carpeta F1_2 a carpeta de trabaj
cp /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F1_2_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F1_2_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
###Ejecutar alineamiento F1_2
hisat2 -p 8 --dta -x PAHCT -1 F1_2_R1_T_PO.fastq.gz -2 F1_2_R2_T_PO.fastq.gz -S F1_2_PAHCT.sam &> F1_2_PAHCT.txt;
###Convertir a archivo Bam F1_2
samtools view -bS F1_2_PAHCT.sam > F1_2_PAHCT.bam;
samtools sort F1_2_PAHCT.bam > F1_2_PAHCT.sorted.bam;
samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta F1_2_PAHCT.sorted.bam | bcftools view -bvcg - > F1_2_PAHCT.raw.bcf;
###Borrar los archivos inecesarios
cp F1_2_PAHCT.sorted.bam F1_2_PAHCT.raw.bcf F1_2_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
rm F*_[123]_R[12]_T_PO.fastq.gz F1_2_R2_T_PO.fastq.gz F1_2_PAHCT.sam F1_2_PAHCT.sorted.bam F1_2_PAHCT.raw.bcf F1_2_PAHCT.txt

### Mover carpeta F1_3 a carpeta de trabaj
cp /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F1_3_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F1_3_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
###Ejecutar alineamiento F1_3
hisat2 -p 8 --dta -x PAHCT -1 F1_3_R1_T_PO.fastq.gz -2 F1_3_R2_T_PO.fastq.gz -S F1_3_PAHCT.sam &> F1_3_PAHCT.txt;
###Convertir a archivo Bam F1_3
samtools view -bS F1_3_PAHCT.sam > F1_3_PAHCT.bam;
samtools sort F1_3_PAHCT.bam > F1_3_PAHCT.sorted.bam;
samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta F1_3_PAHCT.sorted.bam | bcftools view -bvcg - > F1_3_PAHCT.raw.bcf;
###Borrar los archivos inecesarios
cp F1_3_PAHCT.sorted.bam F1_3_PAHCT.raw.bcf F1_3_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
rm F*_[123]_R[12]_T_PO.fastq.gz F1_3_R2_T_PO.fastq.gz F1_3_PAHCT.sam F1_3_PAHCT.sorted.bam F1_3_PAHCT.raw.bcf F1_3_PAHCT.txt


### Mover carpeta F24_1 a carpeta de trabaj
cp /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F24_1_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F24_1_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
###Ejecutar alineamiento F24_1
hisat2 -p 8 --dta -x PAHCT -1 F24_1_R1_T_PO.fastq.gz -2 F24_1_R2_T_PO.fastq.gz -S F24_1_PAHCT.sam &> F24_1_PAHCT.txt;
###Convertir a archivo Bam F24_1
samtools view -bS F24_1_PAHCT.sam > F24_1_PAHCT.bam;
samtools sort F24_1_PAHCT.bam > F24_1_PAHCT.sorted.bam;
samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta F24_1_PAHCT.sorted.bam | bcftools view -bvcg - > F24_1_PAHCT.raw.bcf;
###Borrar los archivos inecesarios
cp F24_1_PAHCT.sorted.bam F24_1_PAHCT.raw.bcf F24_1_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
rm F*_[123]_R[12]_T_PO.fastq.gz F24_1_R2_T_PO.fastq.gz F24_1_PAHCT.sam F24_1_PAHCT.sorted.bam F24_1_PAHCT.raw.bcf F24_1_PAHCT.txt

### Mover carpeta F24_2 a carpeta de trabaj
cp /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F24_2_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F24_2_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
###Ejecutar alineamiento F24_2
hisat2 -p 8 --dta -x PAHCT -1 F24_2_R1_T_PO.fastq.gz -2 F24_2_R2_T_PO.fastq.gz -S F24_2_PAHCT.sam &> F24_2_PAHCT.txt;
###Convertir a archivo Bam F24_2
samtools view -bS F24_2_PAHCT.sam > F24_2_PAHCT.bam;
samtools sort F24_2_PAHCT.bam > F24_2_PAHCT.sorted.bam;
samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta F24_2_PAHCT.sorted.bam | bcftools view -bvcg - > F24_2_PAHCT.raw.bcf;
###Borrar los archivos inecesarios
cp F24_2_PAHCT.sorted.bam F24_2_PAHCT.raw.bcf F24_2_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
rm F*_[123]_R[12]_T_PO.fastq.gz F24_2_R2_T_PO.fastq.gz F24_2_PAHCT.sam F24_2_PAHCT.sorted.bam F24_2_PAHCT.raw.bcf F24_2_PAHCT.txt

#t ### Mover carpeta F24_3 a carpeta de trabaj
#t cp /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F24_3_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F24_3_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
#t ###Ejecutar alineamiento F24_3
#t hisat2 -p 8 --dta -x PAHCT -1 F24_3_R1_T_PO.fastq.gz -2 F24_3_R2_T_PO.fastq.gz -S F24_3_PAHCT.sam &> F24_3_PAHCT.txt;
#t ###Convertir a archivo Bam F24_3
#t samtools view -bS F24_3_PAHCT.sam > F24_3_PAHCT.bam;
#t samtools sort F24_3_PAHCT.bam > F24_3_PAHCT.sorted.bam;
#t samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta F24_3_PAHCT.sorted.bam | bcftools view -bvcg - > F24_3_PAHCT.raw.bcf;
#t ###Borrar los archivos inecesarios
#t cp F24_3_PAHCT.sorted.bam F24_3_PAHCT.raw.bcf F24_3_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
#t rm F*_[123]_R[12]_T_PO.fastq.gz F24_3_R2_T_PO.fastq.gz F24_3_PAHCT.sam F24_3_PAHCT.sorted.bam F24_3_PAHCT.raw.bcf F24_3_PAHCT.txt



### Mover carpeta F72_1 a carpeta de trabaj
cp /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F72_1_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F72_1_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
###Ejecutar alineamiento F72_1
hisat2 -p 8 --dta -x PAHCT -1 F72_1_R1_T_PO.fastq.gz -2 F72_1_R2_T_PO.fastq.gz -S F72_1_PAHCT.sam &> F72_1_PAHCT.txt;
###Convertir a archivo Bam F72_1
samtools view -bS F72_1_PAHCT.sam > F72_1_PAHCT.bam;
samtools sort F72_1_PAHCT.bam > F72_1_PAHCT.sorted.bam;
samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta F72_1_PAHCT.sorted.bam | bcftools view -bvcg - > F72_1_PAHCT.raw.bcf;
###Borrar los archivos inecesarios
cp F72_1_PAHCT.sorted.bam F72_1_PAHCT.raw.bcf F72_1_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
rm F*_[123]_R[12]_T_PO.fastq.gz F72_1_R2_T_PO.fastq.gz F72_1_PAHCT.sam F72_1_PAHCT.sorted.bam F72_1_PAHCT.raw.bcf F72_1_PAHCT.txt

### Mover carpeta F72_2 a carpeta de trabaj
cp /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F72_2_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F72_2_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
###Ejecutar alineamiento F72_2
hisat2 -p 8 --dta -x PAHCT -1 F72_2_R1_T_PO.fastq.gz -2 F72_2_R2_T_PO.fastq.gz -S F72_2_PAHCT.sam &> F72_2_PAHCT.txt;
###Convertir a archivo Bam F72_2
samtools view -bS F72_2_PAHCT.sam > F72_2_PAHCT.bam;
samtools sort F72_2_PAHCT.bam > F72_2_PAHCT.sorted.bam;
samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta F72_2_PAHCT.sorted.bam | bcftools view -bvcg - > F72_2_PAHCT.raw.bcf;
###Borrar los archivos inecesarios
cp F72_2_PAHCT.sorted.bam F72_2_PAHCT.raw.bcf F72_2_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
rm F*_[123]_R[12]_T_PO.fastq.gz F72_2_R2_T_PO.fastq.gz F72_2_PAHCT.sam F72_2_PAHCT.sorted.bam F72_2_PAHCT.raw.bcf F72_2_PAHCT.txt

### Mover carpeta F72_3 a carpeta de trabaj
cp /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F72_3_R1_T_PO.fastq.gz /media/alejandra/HV620S/MuestrasT/FrutoTratado/PO/F72_3_R2_T_PO.fastq.gz /home/alejandra/Escritorio/AlinTes2021;
###Ejecutar alineamiento F72_3
hisat2 -p 8 --dta -x PAHCT -1 F72_3_R1_T_PO.fastq.gz -2 F72_3_R2_T_PO.fastq.gz -S F72_3_PAHCT.sam &> F72_3_PAHCT.txt;
###Convertir a archivo Bam F72_3
samtools view -bS F72_3_PAHCT.sam > F72_3_PAHCT.bam;
samtools sort F72_3_PAHCT.bam > F72_3_PAHCT.sorted.bam;
samtools mpileup -uf Persea_americana_HassCtrinity.Trinity.fasta F72_3_PAHCT.sorted.bam | bcftools view -bvcg - > F72_3_PAHCT.raw.bcf;
###Borrar los archivos inecesarios
cp F72_3_PAHCT.sorted.bam F72_3_PAHCT.raw.bcf F72_3_PAHCT.txt /media/alejandra/AgroBio_Tec/AlineamientoTrinity/;
rm F*_[123]_R[12]_T_PO.fastq.gz F72_3_R2_T_PO.fastq.gz F72_3_PAHCT.sam F72_3_PAHCT.sorted.bam F72_3_PAHCT.raw.bcf F72_3_PAHCT.txt
