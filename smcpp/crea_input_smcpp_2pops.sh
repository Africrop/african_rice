#!/bin/sh

# Before running, ensure it is correctly formatted by using dos2unix script_name.sh
# Make it executable by chmod 755 script_name.sh

# ecrit les erreurs dans le fichier de sortie standard 
#$ -j y 

# shell que l'on veut utiliser 
#$ -S /bin/bash 

# indiquer son email pour suivre l'execution : 
#$ -M philippe.cubry@ird.fr

# obtenir un message au demarrage (b) , a la fin (e), en cas d'abandon (a) 
#$ -m bea 

# la queue que l'on veut utiliser : 
#$ -q bioinfo.q

#$ -N crea_input_smcpp

#Exporte les variables d environnement
#$ -V

# Arguments du scrip : chr à considerer
module load system/python/3.6.0a3

for i in `cat /data2/projects/africrop_models/Rice_smcpp/list_of_pseudo_diploids.txt`;
do
smc++ vcf2smc -m /data2/projects/africrop_models/Rice_smcpp/mask_final.bed.gz -d $i $i /data2/projects/africrop_models/Rice_smcpp/Rice_pseudodiplos_AllChr.vcf.gz /data2/projects/africrop_models/Rice_smcpp/data/chr$1.All.ObOg.smc.masked.$i.gz Chr$1 `cat /data2/projects/africrop_models/Rice_smcpp/Ob_list_smcpp.txt` `cat /data2/projects/africrop_models/Rice_smcpp/Og_list_smcpp.txt`
smc++ vcf2smc -m /data2/projects/africrop_models/Rice_smcpp/mask_final.bed.gz -d $i $i /data2/projects/africrop_models/Rice_smcpp/Rice_pseudodiplos_AllChr.vcf.gz /data2/projects/africrop_models/Rice_smcpp/data/chr$1.All.OgOb.smc.masked.$i.gz Chr$1 `cat /data2/projects/africrop_models/Rice_smcpp/Og_list_smcpp.txt` `cat /data2/projects/africrop_models/Rice_smcpp/Ob_list_smcpp.txt` 
done
