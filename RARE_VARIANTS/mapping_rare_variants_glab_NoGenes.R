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
#source("/Users/cubry/Documents/scripts/quickKrig.R")
source("quickKrig.R")
library(data.table)

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
samples <-  read.table("list-glab-georef.txt",header=TRUE,na.strings = "")

sfs.unfolded.noGenes <- NULL
sfs.folded.noGenes <- NULL
glab.genotypes.noGenes <- NULL
geno1.noGenes <- NULL
geno1.unfolded.noGenes <- NULL
lst.snps.noGenes <- NULL
lst.snps.singletons.folded.noGenes <- NULL


# Read HAPLOIDIZED VCFs
  setwd("/data2/projects/africrop_models/Riz_vcf/")
  input <- "glab_allChr.recode.HAPLO.NEW.No_GENES_HEADER.vcf"

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
  vcf_haplo.noGenes <- fread(file = input, header = FALSE,colClasses =  "character",skip=i)
  colnames(vcf_haplo.noGenes) <- head
  lst.snps.tmp.noGenes <- paste(as.matrix(vcf_haplo.noGenes[,1]),"_",as.matrix(vcf_haplo.noGenes[,2]),sep="")
  vcf_haplo.noGenes <- (vcf_haplo.noGenes[,10:length(head)])
  vcf_haplo.noGenes <- vcf_haplo.noGenes[,as.character(samples$code_vcf),with=FALSE]
  vcf_haplo.noGenes <- as.data.frame(vcf_haplo.noGenes)
  vcf_haplo.noGenes[vcf_haplo.noGenes == "."] <- 9
  vcf_haplo.noGenes <- apply(vcf_haplo.noGenes,MARGIN = 2,as.numeric)

  glab.genotypes.noGenes <- rbind(glab.genotypes.noGenes,vcf_haplo.noGenes)
  lst.snps.noGenes <- c(lst.snps.noGenes,lst.snps.tmp.noGenes)

  # Compute the SFS  
  spect.unfolded.noGenes <- apply(vcf_haplo.noGenes, MARGIN = 1, FUN = daf3)
  plot(table(spect.unfolded.noGenes),xlab="",ylab="",main = paste("unfolded SFS using DAF, chr ",c,sep=""))
  spect.folded.noGenes <- apply(vcf_haplo.noGenes, MARGIN = 1, FUN = maf3)
  plot(table(spect.folded.noGenes),xlab="",ylab="",main = paste("unfolded SFS using DAF, chr ",c,sep=""))
  sfs.unfolded.noGenes <- c(sfs.unfolded.noGenes,spect.unfolded.noGenes)
  sfs.folded.noGenes <- c(sfs.folded.noGenes,spect.folded.noGenes)
  

# Retrieve singletons in the SFS
  lst1.noGenes = which(spect.folded.noGenes == 1)
  geno1.noGenes <- rbind(geno1.noGenes, vcf_haplo.noGenes[lst1.noGenes,])

  lst1.unfolded.noGenes = which(spect.unfolded.noGenes == 1)
  geno1.unfolded.noGenes <- rbind(geno1.unfolded.noGenes,vcf_haplo.noGenes[lst1.unfolded.noGenes,])

  lst.snps.singletons.folded.noGenes <- c(lst.snps.singletons.folded.noGenes,(lst.snps.tmp.noGenes[lst1.noGenes]))

write.table(sfs.folded.noGenes,"sfs.folded.noGenes")
write.table(sfs.unfolded.noGenes,"sfs.unfolded.noGenes")
write.table(geno1.noGenes,"singletons.noGenes")
write.table(geno1.unfolded.noGenes,"singletons.unfolded.noGenes")
write.table(lst.snps.noGenes,"snp.list.glab.noGenes")
write.table(lst.snps.singletons.folded.noGenes,"singletonssnp.list.glab.noGenes")


# We use these matrix to compute the number of rare variants for each genotype in the dataset.
geno1.noGenes <- read.table("singletons.noGenes")
lst.snps.singletons.folded.noGenes <- read.table("singletonssnp.list.glab.noGenes")
row.names(geno1.noGenes) <- lst.snps.singletons.folded.noGenes$x

# We first list the index of which genotype exhibit the rare variants
#for each SNPs in the singletons class of the SFS
ind1.noGenes = apply(geno1.noGenes, MARGIN = 1, FUN = function(x) {
  y = x[x != 9]
  i1 = which(x != 9)[1]
  i = which( (x != y[1]) & (x != 9) )
  if (length(i) > 1) i1 else i 
})

ind1.lab.noGenes = sapply(ind1.noGenes, FUN = function(x) {
  samples[x,"code_vcf"]
})
ind1.lab.noGenes <- as.data.frame(ind1.lab.noGenes)
ind1.lab.noGenes$chr_pos <- row.names(ind1.lab.noGenes)
ind1.lab.noGenes$chr <- apply(ind1.lab.noGenes,MARGIN=1,FUN=function(x){unlist(strsplit(x["chr_pos"],"_"))[1]})
ind1.lab.noGenes$pos <- apply(ind1.lab.noGenes,MARGIN=1,FUN=function(x){unlist(strsplit(x["chr_pos"],"_"))[2]})

ind1.unfolded.noGenes = apply(geno1.unfolded.noGenes, MARGIN = 1, FUN = function(x) {
  y = x[x != 9]
  i1 = which(x != 9)[1]
  i = which( (x != y[1]) & (x != 9) )
  if (length(i) > 1) i1 else i
})

# Then we add all the rare variants for a given genotype
glab.genotypes.mat.noGenes <- as.matrix(t(glab.genotypes.noGenes))
z.noGenes = sapply(1:nrow(glab.genotypes.mat.noGenes),
           FUN = function(x) sum(ind1.noGenes == x) )
z.unfolded.noGenes = sapply(1:nrow(glab.genotypes.mat.noGenes),
                    FUN = function(x) sum(ind1.unfolded.noGenes == x) )

# We finally store these as vectors
rare_variants.noGenes <- as.data.frame(z.noGenes, row.names = colnames(glab.genotypes.noGenes))
rare_variants.unfolded.noGenes <- as.data.frame(z.unfolded.noGenes,
                                        row.names = colnames(glab.genotypes.noGenes))

write.table(rare_variants.noGenes,"rare_variants_glab.noGenes.txt",col.names = FALSE,quote=FALSE)
write.table(rare_variants.unfolded.noGenes,"rare_variants_unfolded_glab.noGenes.txt",col.names = FALSE,quote=FALSE)

# Geographic representation
#We now want to represent the rare variant statistic on geographic map. For that with first build an object with all the statistic we'll want to use, including rare variants calculated here and rare variants exhibited by the chloroplast genomes for the sake of comparisons.

rare_variants.noGenes <- read.table("rare_variants_glab.noGenes.txt",na.strings = "")
rare_variants.unfolded.noGenes <- read.table("rare_variants_unfolded_glab.noGenes.txt",na.strings = "")
colnames(rare_variants.noGenes) <- c("ID","z.noGenes")
colnames(rare_variants.unfolded.noGenes) <- c("ID","z.unfolded.noGenes")
combined <- merge(samples,rare_variants.noGenes,by.x="code_vcf",by.y="ID")
combined <- merge(combined,rare_variants.unfolded.noGenes,by.x="code_vcf",by.y="ID")

#Graphical representation
pdf(file = "rare_variants.noGenes.map.withMR.pdf",paper = "a4r")
par(mfrow = c(1,2))
quickKrig2(coordinates = cbind(combined$Long_,combined$lat_),
          data_to_Krig = combined,datacolumn = "z.noGenes",
          asciifile="africa_india.lowres.asc",
          cell_value_min = 0,
          main="Rare variants for\n cultivated genotypes,\n folded SFS, \n .noGenes",
          xlab="Longitude",
          ylab="Latitude",
          pts.size=.5,
          pts.shape="+",
          xlim=c(-20,60),ylim=c(-40,40))
quickKrig2(coordinates = cbind(combined$Long_,combined$lat_),
          data_to_Krig = combined,datacolumn = "z.unfolded.noGenes",
          asciifile="africa_india.lowres.asc",
          cell_value_min = 0,
          main="Rare variants for\n cultivated genotypes,\n unfolded SFS, \n .noGenes",
          xlab="Longitude",
          ylab="Latitude",
          pts.size=.5,
          pts.shape="+",
          xlim=c(-20,60),ylim=c(-40,40))
dev.off()

# We can now quickly see how rare variants are spanned over all genotypes.
# par(mfrow = c(1,2))
# plot(combined$z,x = combined$code_vcf,las=2,
#      main = "Rare variants\n per genotype,\n folded SFS")
# plot(combined$z.unfolded,x = combined$code_vcf,las=2,
#      main="Rare variants\n per genotype,\n unfolded SFS")
# 
# pdf("rare_var_per_genotypes_glab_folded_unfolded.pdf")
# par(mfrow = c(1,2))
# plot(combined$z,x = combined$code_vcf,las=2,
#      main = "Rare variants\n per genotype,\n folded SFS")
# plot(combined$z.unfolded,x = combined$code_vcf,las=2,
#      main="Rare variants\n per genotype,\n unfolded SFS")
# dev.off()

library(ggplot2)
library(cowplot)

save_plot(filename = "rare_var_per_genotypes_glab_folded.noGenes.pdf",plot = 
qplot(combined$code_vcf,combined$z.noGenes,ylab = "Count of rare variants.noGenes",xlab="Individual labels")+theme(axis.text.x = element_text(angle = 90, hjust = 1)),base_aspect_ratio = 3.5)

# There is one individual which appears to be odd in the repartition, that's the one which gave the dark red point on the previously drawn maps. We try to redraw the maps without it.
pdf(file = "rare_variants.noGenes.map.withoutMR.pdf",paper = "a4r")
par(mfrow = c(1,2))
quickKrig2(coordinates = cbind(combined[-105,]$Long_,combined[-105,]$lat_),
          data_to_Krig = combined[-105,],datacolumn = "z.noGenes",
          asciifile="africa_india.lowres.asc",
          cell_value_min = 0,
          main="Rare variants for\n cultivated genotypes \n .noGenes",
          xlab="Longitude",
          ylab="Latitude",
          pts.size=.5,
          pts.shape="+",
          xlim=c(-20,60),ylim=c(-40,40))
quickKrig2(coordinates = cbind(combined[-105,]$Long_,combined[-105,]$lat_),
          data_to_Krig = combined[-105,],datacolumn = "z.unfolded.noGenes",
          asciifile="africa_india.lowres.asc",
          cell_value_min = 0,
          main="Rare variants for\n cultivated genotypes,\n unfolded SFS \n .noGenes",
          xlab="Longitude",
          ylab="Latitude",
          pts.size=.5,
          pts.shape="+",
          xlim=c(-20,60),ylim=c(-40,40))
dev.off()

pdf(file = "rare_variants.noGenes.map.withoutMR.foldedOnly.pdf",paper = "a4r")
quickKrig2(coordinates = cbind(combined[-105,]$Long_,combined[-105,]$lat_),
           data_to_Krig = combined[-105,],datacolumn = "z.noGenes",
           asciifile="africa_india.lowres.asc",
           cell_value_min = 0,
           main="Singletons for\n cultivated genotypes \n excluding genes",
           xlab="Longitude",
           ylab="Latitude",
           pts.size=.5,
           pts.shape="+",
           xlim=c(-20,60),ylim=c(-40,40))
dev.off()