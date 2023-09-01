#Ensamblado independiente F separado por comas, es la preba solicitada por el doctor
cd /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble/readsOKF;
var=$(date); echo "Comenzó ensamble de F0_ coma[ $var]">> "ResumenE.txt";
$TRINITY_HOME/Trinity --seqType fq --max_memory 120G --left F0_1_PO_fail_f2.R1.fq,F0_2_PO_fail_f2.R1.fq,F0_3_PO_fail_f2.R1.fq --right F0_1_PO_fail_f2.R2.fq,F0_2_PO_fail_f2.R2.fq,F0_3_PO_fail_f2.R2.fq --SS_lib_type FR --CPU 10 --output /home/lbmv/JGMC_copias/FrutoTratado/pruebas_al/hisat2/filtro2/filtradosF/ensamble/F0comatrinity --normalize_max_read_cov 100 &>>F0_TR_coma.txt;
var=$(date); echo "Teminó ensamble de F0_coma[ $var]">> "ResumenE.txt";
