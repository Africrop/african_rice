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

smc++ estimate -o /data2/projects/africrop_models/Rice_smcpp/new_simuls/rec_bott/output/ 6.5e-9 /data2/projects/africrop_models/Rice_smcpp/new_simuls/rec_bott/*.smc.*