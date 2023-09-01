### SCRIPT PARA EL MAPEO DE LAS SECUENCIAS CON EL TRASNCRIPTOMA ENSAMBLADO DE FAAGP Y LAS SECUENCIAS DE XOCA
##21/03/2023
#José Gustavo Marín Contreras
##Duración Aproximada 8 horas por librería
#### Archivos quitosano
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIA_0_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIA_0_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHCRSEM QIA_0_RxHC &>  QIA_0_RxHC.txt;
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIA_1_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIA_1_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHCRSEM QIA_1_RxHC &>  QIA_1_RxHC.txt;
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIA_24_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIA_24_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHCRSEM QIA_24_RxHC &>  QIA_24_RxHC.txt;
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIA_69_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIA_69_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHCRSEM QIA_69_RxHC &>  QIA_69_RxHC.txt;
