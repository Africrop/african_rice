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
# R script for Splatche outputs ABC analysis on cluster 
# will produce summarized statistics of the simulations 
# including binned SFS and rare variants per group      
###################################################################################################################################

## Caution: only useful for cultivated pearl millet data with 146 individuals,
## to apply to other datasets modify SFS binning and maf functions
## requires file group.txt containing groups definition
## Usage : run using qsub from a launching script providing project name (e.g. mil_1pop) and output file name


## Version implemented for dealing with initial variation in pop size (bottleneck)
## Deleted options for Tau and Anc sizes that were not used

#### Gather arguments from script call ####
args <- commandArgs(TRUE)

#### Define useful functions to calculate summarized statistics ####
Nbsample = 112

##### Function to construct SFS #####
myprocess1 = function(
  input= string,
  Nbsample = 112,
  Nbloci = 100,
  Nbindiv = rep(2, 112))
{
  ## process DNA file
  X <- scan(file = input, what = character(), sep = "\t", quiet = TRUE, skip = 0, nlines = 0, comment.char = "#")
  Totindiv <- sum(Nbindiv)
  lx <- length(X)
  X <- X[30:lx]
  M <- data.frame(matrix(NA, nrow = Totindiv, ncol = (2 + Nbloci)))
  for(k in 1:Nbsample){
    shift <- 1
    if (k == 1) k.s <- 0 else k.s <- sum(Nbindiv[1:(k-1)])
    for (i in 1:(Nbindiv[k])){
      X[shift + (i-1)*(2+Nbloci)] <- k
      M[k.s + i, ] <- X[ ( shift + (i-1)*(2+Nbloci) ):(shift -1 +i*(2+Nbloci) )]
      shift <- shift + 1
    }
    lx <- length(X)
    X<-X[(10 + shift+i*(2+Nbloci)):lx]
  }
  return(M[,-2])
}

##### Function to calculate minimal allele frequency #####
maf = function(x){ f = sum(x) ; min(f , length(x) - f ) }

##### Function to bin SFS #####
histo.bin = function(x){
  n = length(x)
  y = sapply(1:56, FUN = function(i) sum(x == i) )
  c(y[1], y[2], y[3], sum(y[4:5]), sum(y[6:10]), sum(y[11:20]), sum(y[20:56]) )/n
}

##### Function to compute mean per group #####
z.groupm = function(z){
  sapply(1:12, FUN = function(i) mean(z[group == i]) )
}

#### Load group definition ####
group = as.numeric(as.vector(read.table("./og_groups_final.txt")[,2]))

#### Variables definition ####
sum.stat = NULL
param = NULL
spect.stat=NULL

#### Loading simulations parameters ####
print("loading parameters")
job_id <- gsub(".tar.gz","",gsub(paste(args[1],"_",sep=""),"",list.files(pattern = args[1])))
for (i in 1:length(job_id)){
  param.simu <- read.table(file = paste("param_",job_id[i],".txt",sep=""),sep="",header=T,colClasses="character")
  param.simu$job_id <- job_id[i]
  param = rbind(param,param.simu)
}

#### Loop for creating summarized statistics from simulated datasets ####
for (n in 1:nrow(param)){
  print(paste("loop",n, "of",nrow(param)))
    ##### creating a snp-containing object #####
  print("extracting files from archive")
  untar(tarfile=paste(args[1],"_",param["job_id"][n,],".tar.gz",sep=""),
        files = c(paste(args[1],"_",param["job_id"][n,],"/GeneticsOutput/settings",param$generations[n],
                      "_",param$accroissement[n],"_",param$migration[n],"_",param$SizeBeforeExpansion[n],
                      "_",param$TimeOfBottleneck[n],"_",param$AncestralSize[n],"_",param$TimeForRecentBott[n],
                      "_",param$MutationRate[n],"_",param$MainCarryingCapacity[n],"_GeneSamples_",seq(1,10,1),".arp" ,sep=""),
                  paste(args[1],"_",param["job_id"][n,],"/Arrival_cell_",param$generations[n],
                              "_",param$accroissement[n],"_",param$migration[n],"_output.txt" ,sep="")),
        list = FALSE, exdir = ".",compressed = "gzip", verbose = FALSE)


  string.snp = paste("./",args[1],"_",param["job_id"][n,],"/GeneticsOutput/settings",
                     param$generations[n],"_",param$accroissement[n],"_",param$migration[n],"_",
                     param$SizeBeforeExpansion[n],"_",param$TimeOfBottleneck[n],"_",
                     param$AncestralSize[n],"_",param$TimeForRecentBott[n],
                     "_",param$MutationRate[n],"_",param$MainCarryingCapacity[n],"_GeneSamples_1.arp" ,sep="")
  genotype = myprocess1(string.snp, Nbsample = 112, Nbloci = 500, Nbindiv = rep(2, 112) )[,-1]

  arrival = (read.table(paste(args[1],"_",param["job_id"][n,],"/Arrival_cell_",param$generations[n],
                  "_",param$accroissement[n],"_",param$migration[n],"_output.txt" ,sep=""),skip=1, sep=":"))
  row.names(arrival) = arrival[,1]; arrival[,1]<-NULL; arrival <- t(arrival)

  for (r in 2:10){
    string1 = paste("./",args[1],"_",param["job_id"][n,],"/GeneticsOutput/settings",param$generations[n],
                    "_",param$accroissement[n],"_",param$migration[n],"_",param$SizeBeforeExpansion[n],
                    "_",param$TimeOfBottleneck[n],"_",param$AncestralSize[n],"_",param$TimeForRecentBott[n],
                    "_",param$MutationRate[n],"_",param$MainCarryingCapacity[n],"_GeneSamples_",r,".arp" ,sep="")
    genotype = cbind(genotype, myprocess1(string1, Nbsample = 112, Nbloci = 500, Nbindiv = rep(2, 112) )[,-1])
  }
  genotype1 = genotype[ (1:224)%%2 == 1, ]
  genotype2 = genotype[ (1:224)%%2 == 0, ]
  genotype = cbind(genotype1, genotype2)

  print("removing temp archive directory")
  unlink(x = paste("./",args[1],"_",param["job_id"][n,],sep=""),recursive = TRUE)

  print("computing summary statistics")


  # Filter MR genotype
  # Load group definition from previous analysis
  group = (read.table("./og_groups_final.txt", na.strings="-9"))
  genotype = genotype[(which(group$V1 != "ND" & group$V1 != "IK"& group$V1 != "EK"& group$V1 != "KC"& group$V1 != "MR")),]
  group = group[(which(group$V1 != "ND" & group$V1 != "IK"& group$V1 != "EK"& group$V1 != "KC"& group$V1 != "MR")),]
  group = as.numeric(as.vector(group[,2]))


  spect = apply(genotype, MARGIN = 2, FUN = function(x) {maf(as.numeric(as.factor(x)) - 1)}   ) #SFS creation

  lst1 = which(spect == 1) # list singletons
  geno1 = as.data.frame(genotype[,lst1])
  ind1 = apply(geno1, MARGIN = 2, FUN = function(x) {
    y = as.numeric(as.factor(x))
    i = which( y != y[1] )
    if (length(i) > 1) 1 else i
  })

  z = sapply(1:107, FUN = function(x) sum(ind1 == x) )

  stat.sim <- c((histo.bin(spect)/sum(histo.bin(spect))), z.groupm(z)/sum(z.groupm(z))) # get summary stats
  names(stat.sim) <- c(paste("SFS",seq(1,7,1), sep=""),paste("RareVariants",seq(1,12,1), sep=""))

  stat.temp <- as.data.frame(c(param[n,],stat.sim,as.data.frame(arrival)))
  spect.temp <- as.data.frame(c(param[n,],spect))

  sum.stat <- rbind(sum.stat, stat.temp )
  spect.stat <- rbind(spect.stat,spect.temp)
  cond <- c(rep("SFS",7),rep("RareVariants",12))

  print("saving computed statistics to file")
}
print("saving computed statistics to file")
write.table(file = paste(args[2],".moinsEst_MR_sum.stat.txt",sep=""), sum.stat, row.names = F, quote=F)
write.table(file = paste(args[2],".moinsEst_MR_sum.spect.stat.txt",sep=""), spect.stat, row.names = F, quote=F)
