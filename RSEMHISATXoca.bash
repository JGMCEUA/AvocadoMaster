### SCRIPT PARA EL MAPEO DE LAS SECUENCIAS CON EL TRASNCRIPTOMA ENSAMBLADO DE FAAGP Y LAS SECUENCIAS DE XOCA
##21/03/2023
#José Gustavo Marín Contreras
##Duración Aproximada 8 horas por librería
#### Archivos control
#t	rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IA_0_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IA_0_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHC IA_0_RxHC &>  IA_0_RxHC.txt;
#t	rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IA_1_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IA_1_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHC IA_1_RxHC &>  IA_1_RxHC.txt;
#t	rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IA_24_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IA_24_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHC IA_24_RxHC &>  IA_24_RxHC.txt;
#t	rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IA_69_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IA_69_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHC IA_69_RxHC &>  IA_69_RxHC.txt;
#### Archivos patogeno
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IC_0_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IC_0_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHC IC_0_RxHC &>  IC_0_RxHC.txt;
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IC_1_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IC_1_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHC IC_1_RxHC &>  IC_1_RxHC.txt;
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IC_24_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IC_24_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHC IC_24_RxHC &>  IC_24_RxHC.txt;
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IC_69_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IC_69_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHC IC_69_RxHC &>  IC_69_RxHC.txt;
#### Archivos quitosano
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIA_0_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIA_0_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHC QIA_0_RxHC &>  QIA_0_RxHC.txt;
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIA_1_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIA_1_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHC QIA_1_RxHC &>  QIA_1_RxHC.txt;
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIA_24_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIA_24_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHC QIA_24_RxHC &>  QIA_24_RxHC.txt;
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIA_69_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIA_69_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHC QIA_69_RxHC &>  QIA_69_RxHC.txt;
#### Archivos quitosano+patogeno
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIC_0_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIC_0_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHC QIC_0_RxHC &>  QIC_0_RxHC.txt;
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIC_1_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIC_1_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHC QIC_1_RxHC &>  QIC_1_RxHC.txt;
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIC_24_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIC_24_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHC QIC_24_RxHC &>  QIC_24_RxHC.txt;
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIC_69_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIC_69_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHC QIC_69_RxHC &>  QIC_69_RxHC.txt;
