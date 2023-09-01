### SCRIPT PARA EL MAPEO DE LAS SECUENCIAS CON EL TRASNCRIPTOMA ENSAMBLADO DE FAAGP Y LAS SECUENCIAS DE XOCA
##21/03/2023
#José Gustavo Marín Contreras
##Duración Aproximada 8 horas por librería
#### Archivos quitosano+patogeno
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIC_0_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIC_0_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHCRSEM QIC_0_RxHC &>  QIC_0_RxHC.txt;
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIC_1_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIC_1_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHCRSEM QIC_1_RxHC &>  QIC_1_RxHC.txt;
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIC_24_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIC_24_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHCRSEM QIC_24_RxHC &>  QIC_24_RxHC.txt;
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIC_69_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIC_69_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHCRSEM QIC_69_RxHC &>  QIC_69_RxHC.txt;
