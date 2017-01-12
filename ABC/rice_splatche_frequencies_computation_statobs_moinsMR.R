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
# Intellectual property belongs to IRD and Grenoble-Alpes University
#
# Written by Philippe Cubry, Yves Vigouroux and Olivier François
#
###################################################################################################################################

######################################################
# R script to compute set of frequencies and summary #
# statistics on real datasets for Splatche analysis  #
######################################################

#### Useful functions ####
# Function to compute minimum allele frequency
maf3 = function(x){ 
  x = x[x != 9]
  f = sum(x) ; (min(f , length(x) - f ))}

# Function to compute missing value frequency
miss.freq = function(x){ 
  n = length(x)
  x = x[x == 9]
  f = length(x)/n}


# Function for binning SFS, CAUTION : written for african rice cultivated genotypes data
histo.bin = function(x){
  n = length(x)
  y = sapply(1:66, FUN = function(i) sum(x == i) )
  c(y[1], y[2], y[3], sum(y[4:5]), sum(y[6:10]), sum(y[11:20]), sum(y[20:66]) )/n
}

#
z.group = function(z){
  sapply(1:13, FUN = function(i) sum(z[group == i]) )
}

#
z.groupm = function(z){
  sapply(1:13, FUN = function(i) mean(z[group == i]) )
}


#Load list of genotypes with geographic coordinates
samples <-  read.table("./list-glab-georef.txt",header=TRUE,na.strings = "")

# Load datasets
genotype <- read.table("~/Riz_data/thin-OgOb-all-MSU7-CHR1.vcf.012",row.names = 1)

for(i in c(2:12)){
  genotype <- cbind(genotype,read.table(paste("~/Riz_data/thin-OgOb-all-MSU7-CHR",
                                              i,".vcf.012",sep=""),row.names = 1))
}

# import the list of genotypes in the order of the vcf files
list_vcf <- read.table("./list_indiv_vcf.txt",colClasses = "character",na.strings = "")

# adding this info to the genotype object
row.names(genotype) <- list_vcf

# Filter glaberrima genotypes
glab.genotypes <- genotype[which(row.names(genotype) %in% samples$code_vcf),]
# Cleaning up
rm(genotype)

# Drop MR genotype
glab.genotypes <- glab.genotypes[which(row.names(glab.genotypes)!="MR"),]

# Load group definition from previous analysis
group = (read.table("./og_groups_final.txt", na.strings="-9"))
group = group[which(group$V1 != "MR"),]
group = as.numeric(as.vector(group[,2]))

# Data transformation
## Tranform df into matrix for faster replacements
glab.genotypes.mat <- as.matrix(glab.genotypes)

## When heterozygote, randomly chose one or the other allele to code as homozygote
glab.genotypes.mat[glab.genotypes.mat == 1] = sample(c(0,2),
                                                     sum(glab.genotypes.mat == 1),
                                                     replace=T)
glab.genotypes.mat[glab.genotypes.mat == -1] = rep(9, sum(glab.genotypes.mat == -1))
glab.genotypes.mat[glab.genotypes.mat == 2] = rep(1, sum(glab.genotypes.mat == 2))



#### Compute SFS on cultivated african ####
#cult.african <- read.table("/Users/cubry/Documents/Mil_analysis/real_dataset_analysis/snps_polym_cultafrican.txt",header=T)
spect = apply(glab.genotypes.mat, MARGIN = 2, FUN = maf3   )
plot(table(spect),xlab="",ylab="",main = "folded SFS using MAF")
# Filtering monomorphic loci
spect.poly <- spect[which(spect!=0)]
plot(table(spect.poly),xlab="",ylab="",main = "folded SFS using MAF",las=2)

write.table(spect,"./rice_folded_sfs_cult.txt")



#### Compute summary statistics
stat.obs <- NULL

lst1 = which(spect == 1)
geno1 = glab.genotypes.mat[,lst1]

ind1 = apply(geno1, MARGIN = 2, FUN = function(x) {
  y = x[x != 9]
  i1 = which(x != 9)[1]
  i = which( (x != y[1]) & (x != 9) )
  if (length(i) > 1) i1 else i 
})

table(ind1) 
z = sapply(1:nrow(glab.genotypes.mat),
           FUN = function(x) sum(ind1 == x) )

stat.obs = as.data.frame(t(c(histo.bin(spect.poly)/sum(histo.bin(spect.poly)), z.groupm(z)/sum(z.groupm(z)))))
colnames(stat.obs) <- c(paste("SFS",seq(1,7,1), sep=""),paste("RareVariants",seq(1,13,1), sep=""))

write.table(stat.obs,"./cult_rice_stat.obs_moinsMR.txt")
