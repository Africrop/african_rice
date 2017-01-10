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
# Written by Philippe Cubry, Olivier François
#
###################################################################################################################################


# Initialisation
source("/Users/cubry/Documents/scripts/quickKrig.R")

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
samples <-  read.table("./list-glab-georef.txt",header=TRUE,na.strings = "")

sfs.unfolded <- NULL
sfs.folded <- NULL
glab.genotypes <- NULL
geno1 <- NULL
geno1.unfolded <- NULL
lst.snps <- NULL
lst.snps.singletons.folded <- NULL


# Read HAPLOIDIZED VCFs
for(c in 1:12){
  setwd("C:/Users/cubry/Documents/Riz_analysis/msmc/VCF chr recoded/")
  input <- paste("chr",c,".recode.HAPLO.NEW.vcf",sep="")
  
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
  lst.snps.tmp <- paste(vcf_haplo[,1],"_",vcf_haplo[,2],sep="")
  vcf_haplo <- vcf_haplo[10:length(head)]
  vcf_haplo <- vcf_haplo[,which(colnames(vcf_haplo) %in% samples$code_vcf)]
  vcf_haplo[vcf_haplo == "."] <- 9
  vcf_haplo <- apply(vcf_haplo,MARGIN = 2,as.numeric)

  glab.genotypes <- rbind(glab.genotypes,vcf_haplo)
  lst.snps <- c(lst.snps,lst.snps.tmp)

  # Compute the SFS  
  spect.unfolded <- apply(vcf_haplo, MARGIN = 1, FUN = daf3)
  plot(table(spect.unfolded),xlab="",ylab="",main = paste("unfolded SFS using DAF, chr ",c,sep=""))
  spect.folded <- apply(vcf_haplo, MARGIN = 1, FUN = maf3)
  plot(table(spect.folded),xlab="",ylab="",main = paste("unfolded SFS using DAF, chr ",c,sep=""))
  sfs.unfolded <- c(sfs.unfolded,spect.unfolded)
  sfs.folded <- c(sfs.folded,spect.folded)
  

# Retrieve singletons in the SFS
  lst1 = which(spect.folded == 1)
  geno1 <- rbind(geno1, vcf_haplo[lst1,])

  lst1.unfolded = which(spect.unfolded == 1)
  geno1.unfolded <- rbind(geno1.unfolded,vcf_haplo[lst1.unfolded,])

  lst.snps.singletons.folded <- c(lst.snps.singletons.folded,(lst.snps.tmp[lst1]))
  }

write.table(sfs.folded,"sfs.folded")
write.table(sfs.unfolded,"sfs.unfolded")
write.table(geno1,"singletons")
write.table(geno1.unfolded,"singletons.unfolded")
write.table(lst.snps,"snp.list.glab")
write.table(lst.snps.singletons.folded,"singletonssnp.list.glab")


# We use these matrix to compute the number of rare variants for each genotype in the dataset.
geno1 <- read.table("singletons")
lst.snps.singletons.folded <- read.table("singletonssnp.list.glab")
row.names(geno1) <- lst.snps.singletons.folded$x

# We first list the index of which genotype exhibit the rare variants
#for each SNPs in the singletons class of the SFS
ind1 = apply(geno1, MARGIN = 1, FUN = function(x) {
  y = x[x != 9]
  i1 = which(x != 9)[1]
  i = which( (x != y[1]) & (x != 9) )
  if (length(i) > 1) i1 else i 
})

ind1.lab = sapply(ind1, FUN = function(x) {
  samples[x,"code_vcf"]
})
ind1.lab <- as.data.frame(ind1.lab)
ind1.lab$chr_pos <- row.names(ind1.lab)
ind1.lab$chr <- apply(ind1.lab,MARGIN=1,FUN=function(x){unlist(strsplit(x["chr_pos"],"_"))[1]})
ind1.lab$pos <- apply(ind1.lab,MARGIN=1,FUN=function(x){unlist(strsplit(x["chr_pos"],"_"))[2]})

ind1.unfolded = apply(geno1.unfolded, MARGIN = 1, FUN = function(x) {
  y = x[x != 9]
  i1 = which(x != 9)[1]
  i = which( (x != y[1]) & (x != 9) )
  if (length(i) > 1) i1 else i
})

# Then we add all the rare variants for a given genotype
glab.genotypes.mat <- as.matrix(t(glab.genotypes))
z = sapply(1:nrow(glab.genotypes.mat),
           FUN = function(x) sum(ind1 == x) )
z.unfolded = sapply(1:nrow(glab.genotypes.mat),
                    FUN = function(x) sum(ind1.unfolded == x) )

# We finally store these as vectors
rare_variants <- as.data.frame(z, row.names = colnames(glab.genotypes))
rare_variants.unfolded <- as.data.frame(z.unfolded,
                                        row.names = colnames(glab.genotypes))

write.table(rare_variants,"rare_variants_glab.txt",col.names = FALSE,quote=FALSE)
write.table(rare_variants.unfolded,"rare_variants_unfolded_glab.txt",col.names = FALSE,quote=FALSE)

# Geographic representation
#We now want to represent the rare variant statistic on geographic map. For that with first build an object with all the statistic we'll want to use, including rare variants calculated here and rare variants exhibited by the chloroplast genomes for the sake of comparisons.

rare_variants <- read.table("rare_variants_glab.txt",na.strings = "")
rare_variants.unfolded <- read.table("rare_variants_unfolded_glab.txt",na.strings = "")
colnames(rare_variants) <- c("ID","z")
colnames(rare_variants.unfolded) <- c("ID","z.unfolded")
combined <- merge(samples,rare_variants,by.x="code_vcf",by.y="ID")
combined <- merge(combined,rare_variants.unfolded,by.x="code_vcf",by.y="ID")

# import rare variants on chlorotypes
rare_variants_chloro <- read.delim(
  "~/scripts/Riz/Riz/Riz_scripts_R/rare_variants_chloro_nora.txt",na.strings = "")

combined <- merge(combined,rare_variants_chloro,by.x="code_vcf",by.y="Code.avec.duplicat")

par(mfrow = c(1,2))
quickKrig2(coordinates = cbind(combined$Long_,combined$lat_),
          data_to_Krig = combined,datacolumn = "z",
          asciifile="/Users/cubry/Documents/GIS/africa_india.lowres.asc",
          cell_value_min = 0,
          main="Rare variants for\n cultivated genotypes,\n folded SFS",
          xlab="Longitude",
          ylab="Latitude",
          pts.size=.5,
          pts.shape="+",
          xlim=c(-20,60),ylim=c(-40,40))
quickKrig2(coordinates = cbind(combined$Long_,combined$lat_),
          data_to_Krig = combined,datacolumn = "z.unfolded",
          asciifile="/Users/cubry/Documents/GIS/africa_india.lowres.asc",
          cell_value_min = 0,
          main="Rare variants for\n cultivated genotypes,\n unfolded SFS",
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

pdf("rare_var_per_genotypes_glab_folded_unfolded.pdf")
par(mfrow = c(1,2))
plot(combined$z,x = combined$code_vcf,las=2,
     main = "Rare variants\n per genotype,\n folded SFS")
plot(combined$z.unfolded,x = combined$code_vcf,las=2,
     main="Rare variants\n per genotype,\n unfolded SFS")
dev.off()

library(ggplot2)
library(cowplot)

save_plot(filename = "rare_var_per_genotypes_glab_folded.pdf",plot = 
qplot(combined$code_vcf,combined$z,ylab = "Count of rare variants",xlab="Individual labels")+theme(axis.text.x = element_text(angle = 90, hjust = 1)),base_aspect_ratio = 3.5)

# There is one individual which appears to be odd in the repartition, that's the one which gave the dark red point on the previously drawn maps. We try to redraw the maps without it.

par(mfrow = c(1,2))
quickKrig2(coordinates = cbind(combined[-98,]$Long_,combined[-98,]$lat_),
          data_to_Krig = combined[-98,],datacolumn = "z",
          asciifile="/Users/cubry/Documents/GIS/africa_india.lowres.asc",
          cell_value_min = 0,
          main="Rare variants for\n cultivated genotypes",
          xlab="Longitude",
          ylab="Latitude",
          pts.size=.5,
          pts.shape="+",
          xlim=c(-20,60),ylim=c(-40,40))
quickKrig2(coordinates = cbind(combined[-98,]$Long_,combined[-98,]$lat_),
          data_to_Krig = combined[-98,],datacolumn = "z.unfolded",
          asciifile="/Users/cubry/Documents/GIS/africa_india.lowres.asc",
          cell_value_min = 0,
          main="Rare variants for\n cultivated genotypes,\n unfolded SFS",
          xlab="Longitude",
          ylab="Latitude",
          pts.size=.5,
          pts.shape="+",
          xlim=c(-20,60),ylim=c(-40,40))