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

for i in {1..12};
do for j in {1..82};
do smc++ vcf2smc -d IND$j IND$j --length 20000000 /data2/projects/africrop_models/Rice_smcpp/simuls/bottleneck/simulated.bottleneck.vcf.gz /data2/projects/africrop_models/Rice_smcpp/simuls/bottleneck/simulated.bottleneck.chr$i.smc.$j.gz Chr$i `cat /data2/projects/africrop_models/Rice_smcpp/list_simulated.smcpp.txt`;
done;
done
