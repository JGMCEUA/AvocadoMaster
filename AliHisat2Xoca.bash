### SCRIPT PARA EL MAPEO DE LAS SECUENCIAS CON EL GENOMA ENSAMBLADO DE FAAGP Y LAS SECUENCIAS DE XOCA
##23/03/2023
#José Gustavo Marín Contreras
##Alieamiento de las secuencias con el genoma de referencia
#### Archivos control
#t hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IA_0_R1_TPO.fastq.gz -2 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IA_0_R2_TPO.fastq.gz -S IA_0.sam &> IA_0_Ali.txt;
#t samtools sort -@ 30 -o IA_0_NCBI.bam  IA_0.sam;
#t rm IA_0.sam;
#t hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IA_1_R1_TPO.fastq.gz -2 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IA_1_R2_TPO.fastq.gz -S IA_1.sam &> IA_1_Ali.txt;
samtools sort -@ 30 -o IA_1_NCBI.bam  IA_1.sam;
rm IA_1.sam;
#t hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IA_24_R1_TPO.fastq.gz -2 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IA_24_R2_TPO.fastq.gz -S IA_24.sam &> IA_24_Ali.txt;
samtools sort -@ 30 -o IA_24_NCBI.bam  IA_24.sam;
rm IA_24.sam;
#t hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IA_69_R1_TPO.fastq.gz -2 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IA_69_R2_TPO.fastq.gz -S IA_69.sam &> IA_69_Ali.txt;
samtools sort -@ 30 -o IA_69_NCBI.bam  IA_69.sam;
rm IA_69.sam;
#### Archivos patógeno
#t hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IC_0_R1_TPO.fastq.gz -2 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IC_0_R2_TPO.fastq.gz -S IC_0.sam &> IC_0_Ali.txt;
samtools sort -@ 30 -o IC_0_NCBI.bam  IC_0.sam;
rm IC_0.sam;
#t hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IC_1_R1_TPO.fastq.gz -2 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IC_1_R2_TPO.fastq.gz -S IC_1.sam &> IC_1_Ali.txt;
samtools sort -@ 30 -o IC_1_NCBI.bam  IC_1.sam;
rm IC_1.sam;
#t hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IC_24_R1_TPO.fastq.gz -2 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IC_24_R2_TPO.fastq.gz -S IC_24.sam &> IC_24_Ali.txt;
samtools sort -@ 30 -o IC_24_NCBI.bam  IC_24.sam;
rm IC_24.sam;
#t hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IC_69_R1_TPO.fastq.gz -2 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/IC_69_R2_TPO.fastq.gz -S IC_69.sam &> IC_69_Ali.txt;
samtools sort -@ 30 -o IC_69_NCBI.bam  IC_69.sam;
rm IC_69.sam;
#### Archivos quitosano
#t hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIA_0_R1_TPO.fastq.gz -2 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIA_0_R2_TPO.fastq.gz -S QIA_0.sam &> QIA_0_Ali.txt;
samtools sort -@ 30 -o QIA_0_NCBI.bam  QIA_0.sam;
rm QIA_0.sam;
#t hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIA_1_R1_TPO.fastq.gz -2 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIA_1_R2_TPO.fastq.gz -S QIA_1.sam &> QIA_1_Ali.txt;
samtools sort -@ 30 -o QIA_1_NCBI.bam  QIA_1.sam;
rm QIA_1.sam;
#t hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIA_24_R1_TPO.fastq.gz -2 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIA_24_R2_TPO.fastq.gz -S QIA_24.sam &> QIA_24_Ali.txt;
samtools sort -@ 30 -o QIA_24_NCBI.bam  QIA_24.sam;
rm QIA_24.sam;
#t hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIA_69_R1_TPO.fastq.gz -2 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIA_69_R2_TPO.fastq.gz -S QIA_69.sam &> QIA_69_Ali.txt;
samtools sort -@ 30 -o QIA_69_NCBI.bam  QIA_69.sam;
rm QIA_69.sam;
#### Archivos quitosano+patógeno
#t hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIC_0_R1_TPO.fastq.gz -2 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIC_0_R2_TPO.fastq.gz -S QIC_0.sam &> QIC_0_Ali.txt;
samtools sort -@ 30 -o QIC_0_NCBI.bam  QIC_0.sam;
rm QIC_0.sam;
#t hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIC_1_R1_TPO.fastq.gz -2 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIC_1_R2_TPO.fastq.gz -S QIC_1.sam &> QIC_1_Ali.txt;
samtools sort -@ 30 -o QIC_1_NCBI.bam  QIC_1.sam;
rm QIC_1.sam;
#t hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIC_24_R1_TPO.fastq.gz -2 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIC_24_R2_TPO.fastq.gz -S QIC_24.sam &> QIC_24_Ali.txt;
samtools sort -@ 30 -o QIC_24_NCBI.bam  QIC_24.sam;
rm QIC_24.sam;
#t hisat2 -p 30 --dta -x /home/alejandra/Escritorio/transcrip_avocado_sequence/GCA_023638045.1/GCA_023638045.1_Gwen_Avocado_Scaffolds_genomic -1 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIC_69_R1_TPO.fastq.gz -2 /home/alejandra/Escritorio/transcrip_avocado_sequence/ResTrim/QIC_69_R2_TPO.fastq.gz -S QIC_69.sam &> QIC_69_Ali.txt;
samtools sort -@ 30 -o QIC_69_NCBI.bam  QIC_69.sam;
rm QIC_69.sam;