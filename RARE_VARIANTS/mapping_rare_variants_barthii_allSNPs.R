###################################################################################################################################
#
# Copyright 2017 IRD
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
# Written by Philippe Cubry, Olivier François
#
###################################################################################################################################

# Following what we've made for glaberrima, we can also plot the rare variants repartition for barthii.

# Initialisation
source("quickKrig.R")

##### Function to compute allele frequency, missing-data wisely #####
# Simple MAF function
maf3 = function(x){ 
  x = x[x != 9]
  f = sum(x) ; min(f , length(x) - f ) }

# Simple DAF function
daf3 = function(x){ 
    x = x[x != 9]
    sum(x)
}

#Load list of genotypes with geographic coordinates
samples <-  read.table("list-barth-georef.txt",header=TRUE,na.strings = "")

# Initalize variables
sfs.unfolded <- NULL
sfs.folded <- NULL
barth.genotypes <- NULL
geno1 <- NULL
geno1.unfolded <- NULL

# Read HAPLOIDIZED VCF
for(c in 1:12){
  input <- paste("chr",c,".barthii.recode.HAPLO.NEW.vcf",sep="")
  input <- paste("chr",c,".barthii.recode.HAPLO.NEW.vcf",sep="")

  #### Open the connection and determine the length (in row) of the header ####
  con <- file(input, open='r')
  i <- 1
  test <- readLines(con,n=1)
  while(grepl(pattern="#", test)==TRUE){
    test <- readLines(con,n=1)
    i <- i+1
  }
  close(con)

  #### Open a new connection to read only the last line of the header ####
  con <- file(input, open='r')
  i <- i-1
  header <- readLines(con,n=i)
  close(con)
  head <- unlist(strsplit(header[length(header)], split="\t"))

#### Read the data and create the header ####
  vcf_haplo <- read.table(file = input, header = FALSE,comment.char = "#",colClasses =  "character")
  colnames(vcf_haplo) <- head
  vcf_haplo <- vcf_haplo[10:length(head)]
  vcf_haplo <- vcf_haplo[,which(colnames(vcf_haplo) %in% samples$code.vcf)]
  vcf_haplo[vcf_haplo == "."] <- 9
  vcf_haplo <- apply(vcf_haplo,MARGIN = 2,as.numeric)

  barth.genotypes <- rbind(barth.genotypes,vcf_haplo)

# Compute the SFS
  spect.unfolded <- apply(vcf_haplo, MARGIN = 1, FUN = daf3)
  spect.folded <- apply(vcf_haplo, MARGIN = 1, FUN = maf3)
  sfs.unfolded <- c(sfs.unfolded,spect.unfolded)
  sfs.folded <- c(sfs.folded,spect.folded)


# Retrieve singletons in the SFS
  lst1 = which(spect.folded == 1)
  geno1 <- rbind(geno1, vcf_haplo[lst1,])

  lst1.unfolded = which(spect.unfolded == 1)
  geno1.unfolded <- rbind(geno1.unfolded,vcf_haplo[lst1.unfolded,])

  }

write.table(sfs.folded,"barthii.sfs.folded")
write.table(sfs.unfolded,"barthii.sfs.unfolded")
write.table(geno1,"barthii.singletons")
write.table(geno1.unfolded,"barthii.singletons.unfolded")


# We use these matrix to compute the number of rare variants for each genotype in the dataset.
# We first list the index of which genotype exhibit the rare variants
#for each SNPs in the singletons class of the SFS
ind1 = apply(geno1, MARGIN = 1, FUN = function(x) {
  y = x[x != 9]
  i1 = which(x != 9)[1]
  i = which( (x != y[1]) & (x != 9) )
  if (length(i) > 1) i1 else i 
})

ind1.unfolded = apply(geno1.unfolded, MARGIN = 1, FUN = function(x) {
  y = x[x != 9]
  i1 = which(x != 9)[1]
  i = which( (x != y[1]) & (x != 9) )
  if (length(i) > 1) i1 else i
})

# Then we add all the rare variants for a given genotype
barth.genotypes.mat <- as.matrix(t(barth.genotypes))

z = sapply(1:nrow(barth.genotypes.mat), FUN = function(x) sum(ind1 == x) )
z.unfolded = sapply(1:nrow(barth.genotypes.mat), FUN = function(x) sum(ind1.unfolded == x) )

# We finally store these as vectors
rare_variants <- as.data.frame(z, row.names = row.names(barth.genotypes.mat))
rare_variants.unfolded <- as.data.frame(z.unfolded, row.names = row.names(barth.genotypes.mat))

write.table(rare_variants,"rare_variants_barth.txt",col.names = FALSE,quote=FALSE)
write.table(rare_variants.unfolded,"rare_variants_unfolded_barth.txt",col.names = FALSE,quote=FALSE)

# Geographic representation
rare_variants <- read.table("rare_variants_barth.txt",na.strings = "")
rare_variants.unfolded <- read.table("rare_variants_unfolded_barth.txt",na.strings = "")
colnames(rare_variants) <- c("ID","z")
colnames(rare_variants.unfolded) <- c("ID","z.unfolded")
combined <- merge(samples,rare_variants,by.x="code_vcf",by.y="ID")
combined <- merge(combined,rare_variants.unfolded,by.x="code_vcf",by.y="ID")

par(mfrow = c(1,2))
quickKrig2(coordinates = cbind(combined$Long_,combined$lat_),
          data_to_Krig = combined,datacolumn = "z",
          asciifile="/Users/cubry/Documents/GIS/africa_india.lowres.asc",
          cell_value_min = 0,
          main="Rare variants for\n wild genotypes,\n folded SFS",
          xlab="Longitude",
          ylab="Latitude",
          pts.size=.5,
          pts.shape="+",
          xlim=c(-20,60),ylim=c(-40,40))
quickKrig2(coordinates = cbind(combined$Long_,combined$lat_),
          data_to_Krig = combined,datacolumn = "z.unfolded",
          asciifile="/Users/cubry/Documents/GIS/africa_india.lowres.asc",
          cell_value_min = 0,
          main="Rare variants for\n wild genotypes,\n unfolded SFS",
          xlab="Longitude",
          ylab="Latitude",
          pts.size=.5,
          pts.shape="+",
          xlim=c(-20,60),ylim=c(-40,40))


# We can now quickly see how rare variants are spanned over all genotypes.
par(mfrow = c(1,2))
plot(combined$z,x = combined$code_vcf,las=2,
     main = "Rare variants\n per genotype,\n folded SFS")
plot(combined$z.unfolded,x = combined$code_vcf,las=2,
     main="Rare variants\n per genotype,\n unfolded SFS")

save_plot(filename = "rare_var_per_genotypes_barth_folded.pdf",plot = 
qplot(combined$code_vcf,combined$z,ylab = "Count of rare variants",xlab="Individual labels")+theme(axis.text.x = element_text(angle = 90, hjust = 1)),base_aspect_ratio = 3.5)

library(ggplot2)
library(cowplot)
rare_var_per_gen <- ggplot(combined,aes(code_vcf,z))+geom_point()+ guides(colour=guide_legend(title=""))+xlab("")+ylab("Number of singletons")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("plot_singletons_per_genotype_barth.pdf",rare_var_per_gen,scale=2)
