###################################################################################################################################
#
# Copyright 2017 IRD and Grenoble-Alpes University
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.
#
# You should have received a copy of the CeCILL-C license with this program.
#If not see <http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt>
#
# Intellectual property belongs to  IRD and Grenoble-Alpes University
#
# Written by Olivier François and Philippe Cubry
#
###################################################################################################################################

###################################################################################################################################
# This script is used to launch Splatche simulation platform
# It works with the script 'splatche_input_files_creation_final_newpriorsB_symbolicLink.R
# A working version of Splatche is needed alongside these scripts
###################################################################################################################################

# Usage : splatche_launcher.R <run number> <number of source populations>

#### getting Job_Id ####
job <- paste(Sys.getenv("JOB_ID"),"-",Sys.getenv("SGE_TASK_ID"),sep="")

#### getting call arguments ####
args <- commandArgs(TRUE)

#### source useful file ####
source("./splatche_input_files_creation_final_newpriorsB.R")

#### Definition du nombre de runs et de populations sources si specifie dans l'appel du script (defaut a 1 et 1) ####
if (is.na(args[1])==TRUE){
  nrun <- 1 #nombre de runs
  norigin <- 1 #nombre de populations sources
} else {nrun<-as.numeric(args[1]); norigin<-as.numeric(args[2])} #valeurs specifiees dans l'appel du script

#### Boucle sur le nombre de runs ####
for (it in 1:nrun){
  print(paste(it,"over",nrun,sep=" ")) #compteur de run

  ##### tirage des coordonnees des populations sources #####
  coordo <- NULL #Initialisation d'une variable de stockage des coordonnees des populations sources

  for (ori in 1:norigin){
    ox = runif(1, -16, 40) #tirage de la longitude
    oy = runif(1, 0, 30) #tirage de la latitude

    while(!is.african(cbind(ox,oy))){ #test des coordonnees pour verifier qu'on est pas dans l'eau, sinon on retire
      ox = runif(1, -16, 40)
      oy = runif(1, 0, 30)
    }
    coordo = rbind(coordo,c(ox,oy)) #stockage des coordonnees tirees
  }

  ##### tirage des variables #####
  ng <- sample(gen_min:gen_max,1) #tirage du nombre de generations a simuler
  rate <- runif(1,acc_min,acc_max) #tirage du taux d'accroissement des populations
  mig <- runif(1,mig_min,mig_max) #tirage du taux de migration
  res <- sample(res_min:res_max,1) #tirage de la taille avant expansion
  Tres <- sample(Tres_min:Tres_max,1) #tirage de la taille avant expansion
  resB <- sample(resB_min:resB_max,1) #tirage de la taille avant expansion
  MutRate<- sample(MutRate_prior,1)
  K <- sample(K_min:K_max,1)

  TrecBott <- NULL
  TrecBott <- sample(TrecBott_min:ng,1)
  while(!TrecBott<= ng){
  TrecBott <- sample(TrecBott_min:TrecBott_max,1) #tirage de la generation de collapse du Ne avec condition
  }

  splatche(input = "settings.txt", ng, rate, mig, norigin,res,Tres,resB,TrecBott, MutRate, K, coord.o = coordo) #appel de la fonction Splatche
  system(paste("splatche2-01 settings",ng,"_",rate,"_",mig,"_",res,"_",
               Tres,"_",resB,"_",TrecBott,"_",MutRate,"_",K, ".txt", sep="")) #lancement de splatche avec les parametres tires

  ##### stockage des parametres uniquement si Splatche a genere des fichiers de sortie #####
  if(file.exists(paste("./datasets_1layer/GeneticsOutput/settings",ng,"_",rate,"_",mig,"_",res,"_",
                       Tres,"_",resB,"_",TrecBott,"_",MutRate,"_",K,"_GeneSamples_2.arp" ,sep=""))){

    coo1 <- NULL ; noms <- NULL
    for (npop in 1:norigin){
      coo1<-cbind(coo1,coordo[npop,1],coordo[npop,2])
      noms <-cbind(noms,paste("Long",npop,sep=""),paste("Lat",npop,sep=""))
    }
    param<- rbind(param,c(coo1,ng,rate,mig,res,Tres,resB,TrecBott,MutRate,K)) #stockage des parametres dans un fichier
    colnames(param) <- c(noms,"generations","accroissement","migration","SizeBeforeExpansion",
                         "TimeOfBottleneck","AncestralSize","TimeForRecentBott","MutationRate","MainCarryingCapacity")
  }
  ##### Sauvegarde des parametres des runs realises dans un fichier #####
  write.table(file = paste("param_",job,".txt",sep=""), param, row.names = F, quote=F)
}
