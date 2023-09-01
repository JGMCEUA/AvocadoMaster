### SCRIPT PARA RECORTE DE CALIDAD DE LAS SECUENCIAS DE XOCA CON TRIMMOMATIC
##21/03/2023
#José Gustavo Marín Contreras

#_______________________________________________________________________________
#				CARPETA QUITOSANO+PATOGENO
#_______________________________________________________________________________

#####SECUENCIAS QIC_0
trimmomatic PE -phred33 /home/alejandra/Escritorio/transcrip_avocado_sequence/quitosano+patogeno/QIC_0_R1.fastq /home/alejandra/Escritorio/transcrip_avocado_sequence/quitosano+patogeno/QIC_0_R1.fastq /home/alejandra/Escritorio/transcrip_avocado_sequence/quitosano+patogeno/Resultados/QIC_0_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/quitosano+patogeno/Resultados/QIC_0_R1_TUP.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/quitosano+patogeno/Resultados/QIC_0_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/quitosano+patogeno/Resultados/QIC_0_R2_TUP.fastq.gz ILLUMINACLIP:/home/alejandra/Escritorio/transcrip_avocado_sequence/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 &> TrimmomaticQIC_0.txt;

#####SECUENCIAS QIC_1
trimmomatic PE -phred33 /home/alejandra/Escritorio/transcrip_avocado_sequence/quitosano+patogeno/QIC_1_R1.fastq /home/alejandra/Escritorio/transcrip_avocado_sequence/quitosano+patogeno/QIC_1_R1.fastq /home/alejandra/Escritorio/transcrip_avocado_sequence/quitosano+patogeno/Resultados/QIC_1_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/quitosano+patogeno/Resultados/QIC_1_R1_TUP.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/quitosano+patogeno/Resultados/QIC_1_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/quitosano+patogeno/Resultados/QIC_1_R2_TUP.fastq.gz ILLUMINACLIP:/home/alejandra/Escritorio/transcrip_avocado_sequence/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 &> TrimmomaticQIC_1.txt;

#####SECUENCIAS QIC_24
trimmomatic PE -phred33 /home/alejandra/Escritorio/transcrip_avocado_sequence/quitosano+patogeno/QIC_24_R1.fastq /home/alejandra/Escritorio/transcrip_avocado_sequence/quitosano+patogeno/QIC_24_R1.fastq /home/alejandra/Escritorio/transcrip_avocado_sequence/quitosano+patogeno/Resultados/QIC_24_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/quitosano+patogeno/Resultados/QIC_24_R1_TUP.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/quitosano+patogeno/Resultados/QIC_24_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/quitosano+patogeno/Resultados/QIC_24_R2_TUP.fastq.gz ILLUMINACLIP:/home/alejandra/Escritorio/transcrip_avocado_sequence/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 &> TrimmomaticQIC_24.txt;

#####SECUENCIAS QIC_69
trimmomatic PE -phred33 /home/alejandra/Escritorio/transcrip_avocado_sequence/quitosano+patogeno/QIC_69_R1.fastq /home/alejandra/Escritorio/transcrip_avocado_sequence/quitosano+patogeno/QIC_69_R1.fastq /home/alejandra/Escritorio/transcrip_avocado_sequence/quitosano+patogeno/Resultados/QIC_69_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/quitosano+patogeno/Resultados/QIC_69_R1_TUP.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/quitosano+patogeno/Resultados/QIC_69_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/quitosano+patogeno/Resultados/QIC_69_R2_TUP.fastq.gz ILLUMINACLIP:/home/alejandra/Escritorio/transcrip_avocado_sequence/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 &> TrimmomaticQIC_69.txt;
