#!/bin/bash
#scrip para meter todos lo datos a una base de datos
# para correr sh ./intomysqlfor_counts.sh muestra
database=$1

# determinar la ruta y nombre del script
script=$(cd ${0%/*} && echo $PWD/${0##*/})
dirscript=`dirname "$script"`

for library in `ls *.txt | sed 's/.txt$//g'`
do
#crear tablas
echo "creando tablas"
#mysql -u xocanol -p123 -e
mysql -h localhost -u root -p123 -e "use $database ; create table $library(
IDE varchar(20),
count_$library decimal(20),
index gen(IDE));"

#haciendo inputs
#mysql -u xocanol -p123 -e
#mysql -h localhost -u root -p123 -e
#mysql -h localhost -u root -p123 -e --local-infile
mysql -h localhost -u root -p123 -e --local-infile "use $database; load data local infile '$dirscript/$library.txt' into table $library;"
echo "tablas indexadas $library"

done

#revisar el error, no  deja cargar la tablas, posiblemente sea problemas de permisos del mysql




