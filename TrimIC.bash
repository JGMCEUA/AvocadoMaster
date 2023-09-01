### SCRIPT PARA RECORTE DE CALIDAD DE LAS SECUENCIAS DE XOCA CON TRIMMOMATIC
##21/03/2023
#José Gustavo Marín Contreras

#_______________________________________________________________________________
#				CARPETA patogeno
#_______________________________________________________________________________

#####SECUENCIAS IC_0
trimmomatic PE -phred33 /home/alejandra/Escritorio/transcrip_avocado_sequence/patogeno/IC_0_R1.fastq /home/alejandra/Escritorio/transcrip_avocado_sequence/patogeno/IC_0_R1.fastq /home/alejandra/Escritorio/transcrip_avocado_sequence/patogeno/Resultados/IC_0_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/patogeno/Resultados/IC_0_R1_TUP.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/patogeno/Resultados/IC_0_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/patogeno/Resultados/IC_0_R2_TUP.fastq.gz ILLUMINACLIP:/home/alejandra/Escritorio/transcrip_avocado_sequence/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 &> TrimmomaticIC_0.txt;

#####SECUENCIAS IC_1
trimmomatic PE -phred33 /home/alejandra/Escritorio/transcrip_avocado_sequence/patogeno/IC_1_R1.fastq /home/alejandra/Escritorio/transcrip_avocado_sequence/patogeno/IC_1_R1.fastq /home/alejandra/Escritorio/transcrip_avocado_sequence/patogeno/Resultados/IC_1_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/patogeno/Resultados/IC_1_R1_TUP.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/patogeno/Resultados/IC_1_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/patogeno/Resultados/IC_1_R2_TUP.fastq.gz ILLUMINACLIP:/home/alejandra/Escritorio/transcrip_avocado_sequence/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 &> TrimmomaticIC_1.txt;

#####SECUENCIAS IC_24
trimmomatic PE -phred33 /home/alejandra/Escritorio/transcrip_avocado_sequence/patogeno/IC_24_R1.fastq /home/alejandra/Escritorio/transcrip_avocado_sequence/patogeno/IC_24_R1.fastq /home/alejandra/Escritorio/transcrip_avocado_sequence/patogeno/Resultados/IC_24_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/patogeno/Resultados/IC_24_R1_TUP.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/patogeno/Resultados/IC_24_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/patogeno/Resultados/IC_24_R2_TUP.fastq.gz ILLUMINACLIP:/home/alejandra/Escritorio/transcrip_avocado_sequence/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 &> TrimmomaticIC_24.txt;

#####SECUENCIAS IC_69
trimmomatic PE -phred33 /home/alejandra/Escritorio/transcrip_avocado_sequence/patogeno/IC_69_R1.fastq /home/alejandra/Escritorio/transcrip_avocado_sequence/patogeno/IC_69_R1.fastq /home/alejandra/Escritorio/transcrip_avocado_sequence/patogeno/Resultados/IC_69_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/patogeno/Resultados/IC_69_R1_TUP.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/patogeno/Resultados/IC_69_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/patogeno/Resultados/IC_69_R2_TUP.fastq.gz ILLUMINACLIP:/home/alejandra/Escritorio/transcrip_avocado_sequence/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 &> TrimmomaticIC_69.txt;
