#!/bin/bash

databases=$1
experiment=$2


#haciendo tablas final

mysql -h localhost -u root -p123 -e "use $databases;create table experiment_final( 
id_gen varchar(20), 
$experiment float (8,3), 
pv$experiment decimal(20,8), 
index gen(id_gen));
insert into experiment_final select IDE,Fold$experiment,FDR$experiment from $experiment;"

#rm $experiment.txt

#Para los merge
for library in `ls *.txt | sed 's/.txt$//g'`

do 
echo "haciendo merge fold"
mysql -h localhost -u root -p123 -e "use $databases;alter table experiment_final add column $library float (8,3);
update experiment_final, $library set experiment_final.$library=$library.Fold$library where experiment_final.id_gen=$library.IDE;"

echo "haciendo merge pvalue"
mysql -h localhost -u root -p123 -e "use $databases;alter table experiment_final add column pv$library float (20,8);
update experiment_final,$library set experiment_final.pv$library=$library.FDR$library where experiment_final.id_gen=$library.IDE;"
echo "fin $library"

done

#Para meter tablas de identificacion(archivo con la descripcion de los genes)
#mysql -h localhost -u root -p123 -e "use $databases;create table description( code varchar(30),locus varchar(25), alias text, descript text, index id(locus));
#load data local infile '/home/sandra/Documents/BGI/pepper_data/Pepper.gene.symbol' into table description;"

#mysql -h localhost -u root -p123 -e "use $databases;alter table experiment_final add column descrip text after id_gen;
#update experiment_final,description set experiment_final.descrip=description.descript where experiment_final.id_gen=description.locus;"


#mysql -h localhost -u root -p123 -e "use $databases; alter table experiment_final add column symbol text after descrip;
#update experiment_final,description set experiment_final.symbol=description.alias where experiment_final.id_gen=description.locus;"


#Para sacar la tabla en archivo csv (desde la tereminal)
#mysql -u xocanol -e "use avocado; SELECT * FROM  experiment_final" > /home/xocanol/samples/expe.csv



