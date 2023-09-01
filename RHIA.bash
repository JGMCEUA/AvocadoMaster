### SCRIPT PARA EL MAPEO DE LAS SECUENCIAS CON EL TRASNCRIPTOMA ENSAMBLADO DE FAAGP Y LAS SECUENCIAS DE XOCA
##21/03/2023
#José Gustavo Marín Contreras
##Duración Aproximada 8 horas por librería
#### Archivos control
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IA_0_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IA_0_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHC IA_0_RxHC &>  IA_0_RxHC.txt;
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IA_1_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IA_1_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHC IA_1_RxHC &>  IA_1_RxHC.txt;
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IA_24_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IA_24_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHC IA_24_RxHC &>  IA_24_RxHC.txt;
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IA_69_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IA_69_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHC IA_69_RxHC &>  IA_69_RxHC.txt;
