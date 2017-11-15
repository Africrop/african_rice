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

#$ -N rice_multisfs

###### Definition des variables de chemin

path_to_dir="/data2/projects/africrop_models/Riz_fastsimcoal";
path_to_tmp="/scratch/cubry-$JOB_ID-$SGE_TASK_ID" 

###### Creation du repertoire temporaire sur noeud

mkdir $path_to_tmp

scp -rp nas:/$path_to_dir/allChr.recode.HAPLO.NEW.vcf $path_to_tmp/
scp -rp nas:/$path_to_dir/List_barth_all.txt $path_to_tmp/
scp -rp nas:/$path_to_dir/List_glab_all.txt $path_to_tmp/
scp -rp nas:/$path_to_dir/functions_daf.R $path_to_tmp/
scp -rp nas:/$path_to_dir/creation_sfs2d_rice_real_data_withmonomorphic.R $path_to_tmp/


echo "tranfert donnees master -> noeud";
cd $path_to_tmp

###### Execution du programme
Rscript $path_to_tmp/creation_sfs2d_rice_real_data_withmonomorphic.R

##### Transfert des données du noeud vers master
rcp -rp $path_to_tmp/* nas:/$path_to_dir/

echo "Transfert donnees node -> master";


#### Suppression du repertoire tmp noeud

rm -rf $path_to_tmp

echo "Suppression des donnees sur le noeud";
