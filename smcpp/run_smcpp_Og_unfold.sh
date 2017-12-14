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

#$ -N smcpp

#Exporte les variables d environnement
#$ -V
module load system/python/3.6.0a3

smc++ estimate -o /data2/projects/africrop_models/Rice_smcpp/output_og_unfold/ --unfold --tK 100000 6.5e-9 /data2/projects/africrop_models/Rice_smcpp/data/*og.smc.masked*
smc++ estimate -o /data2/projects/africrop_models/Rice_smcpp/output_og_unfold_regpen6/ --unfold --tK 100000 --regularization-penalty 6 6.5e-9 /data2/projects/africrop_models/Rice_smcpp/data/*og.smc.masked*