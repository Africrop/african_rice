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

#### Definie run numbers and number of source populations ####
if (is.na(args[1])==TRUE){
  nrun <- 1 #nombre de runs
  norigin <- 1 #nombre de populations sources
} else {nrun<-as.numeric(args[1]); norigin<-as.numeric(args[2])} #valeurs specifiees dans l'appel du script

#### loop over runs ####
for (it in 1:nrun){
  print(paste(it,"over",nrun,sep=" ")) #compteur de run

  ##### sampling coordinates #####
  coordo <- NULL

  for (ori in 1:norigin){
    ox = runif(1, -16, 40) # longitude
    oy = runif(1, 0, 30) # latitude

    while(!is.african(cbind(ox,oy))){
      ox = runif(1, -16, 40)
      oy = runif(1, 0, 30)
    }
    coordo = rbind(coordo,c(ox,oy)) # store coordinates
  }

  ##### sampling variables #####
  ng <- sample(gen_min:gen_max,1) # generations number
  rate <- runif(1,acc_min,acc_max) # growth rate
  mig <- runif(1,mig_min,mig_max) # migration rate
  res <- sample(res_min:res_max,1) # size just before expansion
  Tres <- sample(Tres_min:Tres_max,1) # length of pre-expansion bottleneck
  resB <- sample(resB_min:resB_max,1) # size before the pre-expansion bottleneck
  MutRate<- sample(MutRate_prior,1) # mutation rate
  K <- sample(K_min:K_max,1) # carrying capacity

  TrecBott <- NULL # time for recent collapse of effective size
  TrecBott <- sample(TrecBott_min:ng,1)
  while(!TrecBott<= ng){
  TrecBott <- sample(TrecBott_min:TrecBott_max,1)
  }

  splatche(input = "settings.txt", ng, rate, mig, norigin,res,Tres,resB,TrecBott, MutRate, K, coord.o = coordo) # call Splatche function
  system(paste("splatche2-01 settings",ng,"_",rate,"_",mig,"_",res,"_",
               Tres,"_",resB,"_",TrecBott,"_",MutRate,"_",K, ".txt", sep="")) #launch Splatche with retained parameters

  ##### store parameters only if Splatche generated output files #####
  if(file.exists(paste("./datasets_1layer/GeneticsOutput/settings",ng,"_",rate,"_",mig,"_",res,"_",
                       Tres,"_",resB,"_",TrecBott,"_",MutRate,"_",K,"_GeneSamples_2.arp" ,sep=""))){

    coo1 <- NULL ; noms <- NULL
    for (npop in 1:norigin){
      coo1<-cbind(coo1,coordo[npop,1],coordo[npop,2])
      noms <-cbind(noms,paste("Long",npop,sep=""),paste("Lat",npop,sep=""))
    }
    param<- rbind(param,c(coo1,ng,rate,mig,res,Tres,resB,TrecBott,MutRate,K))
    colnames(param) <- c(noms,"generations","accroissement","migration","SizeBeforeExpansion",
                         "TimeOfBottleneck","AncestralSize","TimeForRecentBott","MutationRate","MainCarryingCapacity")
  }
  ##### Save parameters #####
  write.table(file = paste("param_",job,".txt",sep=""), param, row.names = F, quote=F)
}
