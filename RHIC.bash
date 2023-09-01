### SCRIPT PARA EL MAPEO DE LAS SECUENCIAS CON EL TRASNCRIPTOMA ENSAMBLADO DE FAAGP Y LAS SECUENCIAS DE XOCA
##21/03/2023
#José Gustavo Marín Contreras
##Duración Aproximada 8 horas por librería
#### Archivos control
#### Archivos patogeno
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IC_0_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IC_0_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHCRSEM IC_0_RxHC &>  IC_0_RxHC.txt;
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IC_1_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IC_1_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHCRSEM IC_1_RxHC &>  IC_1_RxHC.txt;
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IC_24_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IC_24_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHCRSEM IC_24_RxHC &>  IC_24_RxHC.txt;
rsem-calculate-expression --paired-end --hisat2-hca --hisat2-path /usr/local/bin --paired-end -p 30 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IC_69_R1_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IC_69_R2_TPO.fastq.gz /home/alejandra/Escritorio/transcrip_avocado_sequence/Index/PAHCRSEM IC_69_RxHC &>  IC_69_RxHC.txt;
