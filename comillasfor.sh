#!/bin/bash
#scrip para borrar comillas de archivos provenientes de R
# para correr sh ./comillasfor.sh
# copiarlo en la carpeta donde estan los archivos txt a modificar

for archivo in  `ls *.txt`
do
sed 's/"//g' $archivo > b$archivo
mv b$archivo $archivo
done




