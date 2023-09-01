### SCRIPT PARA RECORTE DE CALIDAD DE LAS SECUENCIAS DE XOCA CON TRIMMOMATIC
##21/03/2023
#José Gustavo Marín Contreras

#_______________________________________________________________________________
#				CARPETA CONTROL
#_______________________________________________________________________________

#####SECUENCIAS IA_0
#t	trimmomatic PE -phred33 /home/alejandra/Escritorio/transcrip_avocado_sequence/control/IA_0_R1.fastq /home/alejandra/Escritorio/transcrip_avocado_sequence/control/IA_0_R1.fastq /home/alejandra/Escritorio/transcrip_avocado_sequence/control/Resultados/IA_0_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/control/Resultados/IA_0_R1_TUP.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/control/Resultados/IA_0_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/control/Resultados/IA_0_R2_TUP.fastq.gz ILLUMINACLIP:/home/alejandra/Escritorio/transcrip_avocado_sequence/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 &> TrimmomaticIA_0.txt;

#####SECUENCIAS IA_1
trimmomatic PE -threads 30 -phred33 /home/alejandra/Escritorio/transcrip_avocado_sequence/control/IA_1_R1.fastq /home/alejandra/Escritorio/transcrip_avocado_sequence/control/IA_1_R1.fastq /home/alejandra/Escritorio/transcrip_avocado_sequence/control/Resultados/IA_1_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/control/Resultados/IA_1_R1_TUP.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/control/Resultados/IA_1_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/control/Resultados/IA_1_R2_TUP.fastq.gz ILLUMINACLIP:/home/alejandra/Escritorio/transcrip_avocado_sequence/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 &> TrimmomaticIA_1.txt;

#####SECUENCIAS IA_24
#t	trimmomatic PE -phred33 /home/alejandra/Escritorio/transcrip_avocado_sequence/control/IA_24_R1.fastq /home/alejandra/Escritorio/transcrip_avocado_sequence/control/IA_24_R1.fastq /home/alejandra/Escritorio/transcrip_avocado_sequence/control/Resultados/IA_24_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/control/Resultados/IA_24_R1_TUP.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/control/Resultados/IA_24_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/control/Resultados/IA_24_R2_TUP.fastq.gz ILLUMINACLIP:/home/alejandra/Escritorio/transcrip_avocado_sequence/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 &> TrimmomaticIA_24.txt;

#####SECUENCIAS IA_69
trimmomatic PE -threads 30 -phred33 /home/alejandra/Escritorio/transcrip_avocado_sequence/control/IA_69_R1.fastq /home/alejandra/Escritorio/transcrip_avocado_sequence/control/IA_69_R1.fastq /home/alejandra/Escritorio/transcrip_avocado_sequence/control/Resultados/IA_69_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/control/Resultados/IA_69_R1_TUP.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/control/Resultados/IA_69_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/control/Resultados/IA_69_R2_TUP.fastq.gz ILLUMINACLIP:/home/alejandra/Escritorio/transcrip_avocado_sequence/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 &> TrimmomaticIA_69.txt;

