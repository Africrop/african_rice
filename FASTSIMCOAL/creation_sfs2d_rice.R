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
# Intellectual property belongs to IRD 
#
# Written by Philippe Cubry, Yves Vigouroux
#
###################################################################################################################################

###############################################
# Script to generate 2D-SFS from lfmmb data   #
# lfmmb = format identical to lfmm but with   #
# headers and row names for individual and    #
# SNP name information. This script use a     #
# sub-sampling method to estimate 2D-SFS      #
# from datasets with missing data (threshold) #
# to be adjusted.                             #
# Auth: Philippe Cubry with ideas from Yves   #
# Vigouroux and Olivier Francois.             #
# April 2016                                  #
# Requirement: a script to source with        #
# functions used - a script that define       #
# groups to be considered - a file with SNP   #
# info (lfmmb format but can easily be adapted#
# to other formats)                           #
# CAUTION : Version to estimate unfolded SFS  #
# based on polarised SNPs                     #
###############################################

library(ggplot2)
library(reshape2)
library(gridExtra)
library(gtable)
library(cowplot)
library(LEA)
library(grid)

# Load SNPs data from polarized biallelic SNPs file
snps.geno <- NULL
print("read dataset 1")
# temp <- read.table("/Users/cubry/Documents/Mil_data/genotypes_215/polarized_biallelic/chr1_215samples_polarized_SNPs_data_biallelic.txt",header=T)
temp <- read.table("/data2/projects/africrop_models/polarized_snps_mil/polarized_biallelic/chr1_215samples_polarized_SNPs_data_biallelic_CORRECTED.txt",header=T)
row.names(temp) <- paste(temp[,1],"_",temp[,2],sep="")
temp <- temp[,-c(1,2)]
snps.geno <- temp ; rm(temp)

print("read other datasets and combine")
for(i in 2:7){
  # temp <- read.table(paste("/Users/cubry/Documents/Mil_data/genotypes_215/polarized_biallelic/chr",i,"_215samples_polarized_SNPs_data_biallelic.txt",sep=""),header=T)
  temp <- read.table(paste("/data2/projects/africrop_models/polarized_snps_mil/polarized_biallelic/chr",i,"_215samples_polarized_SNPs_data_biallelic_CORRECTED.txt",sep=""),header=T)
  row.names(temp) <- paste(temp[,1],"_",temp[,2],sep="")
  temp <- temp[,-c(1,2)]
  snps.geno <- rbind(snps.geno,temp) ; rm(temp)
}
# Source used R scripts
# source("/Users/cubry/Documents/scripts/fastsimcoal/R/groups_definition.R") # Load groups definition
# source("/Users/cubry/Documents/scripts/fastsimcoal/R/functions.R") # Load used functions
source("/home/cubry/scripts/prepare_data/groups_definition.R") # Load groups definition
source("/home/cubry/scripts/prepare_data/functions.R") # Load used functions

#Sampling individuals within cultivated groups
cult_c <- sample(cult_c,size = 25,replace = FALSE)
cult_i <- sample(cult_i,size = 25,replace = FALSE)
cult_s <- sample(cult_s,size = 25,replace = FALSE)
write.table(cult_c,"/data2/projects/africrop_models/polarized_snps_mil/results_b/cult_c_sampled.txt")
write.table(cult_i,"/data2/projects/africrop_models/polarized_snps_mil/results_b/cult_i_sampled.txt")
write.table(cult_s,"/data2/projects/africrop_models/polarized_snps_mil/results_b/cult_s_sampled.txt")

# Saving data
print("save lfmm and geno objects")
write.lfmm(t(snps.geno),"/data2/projects/africrop_models/polarized_snps_mil/results_b/snps_polarized_biallelic_215_CORRECTED_sampling.lfmm")
lfmm2geno("/data2/projects/africrop_models/polarized_snps_mil/results_b/snps_polarized_biallelic_215_CORRECTED_sampling.lfmm")
# write.lfmm(t(snps.geno),"snps_polarized_biallelic_215.lfmm")
# lfmm2geno("snps_polarized_biallelic_215.lfmm")

# reformat data
print("reformat data")
snps.lfmm <- t(snps.geno) ; rm(snps.geno)
snps.lfmm <- merge(snps.lfmm,all_ind,by.x="row.names",by.y="Id",sort=F)
row.names(snps.lfmm) <- snps.lfmm$Row.names ; snps.lfmm$Row.names<-NULL

# Retaining only polymorphic loci over the set of 215 genotypes
print("filtering monomorphic loci over the whole dataset")
lst <- which(lapply(apply(snps.lfmm[,-length(snps.lfmm)],2,function(x){x=x[x!=9] ; unique(x)}),length)!=1) #List polymorphic loci
sfs <- apply(snps.lfmm[,lst],2,daf4)
snps.poly <- names(snps.lfmm[,lst])
print("save list of polymorphic loci")
write.table(snps.poly,"/data2/projects/africrop_models/polarized_snps_mil/results_b/snps.poly_CORRECTED_sampling.txt")

# Calculating SFS on cultivated genotypes
print("calculating cultivated genotypes SFS")
lst.cult <- which(row.names(snps.lfmm) %in% as.character(all_ind[-grep(all_ind[,2],pattern = "wild"),]$Id))
lst.cult.poly <- which(lapply(apply(snps.lfmm[lst.cult,-length(snps.lfmm)],2,function(x){x=x[x!=9] ; unique(x)}),length)!=1) #List polymorphic loci
sfs.cult.poly <- sfs.daf4(snps.lfmm[lst.cult,lst.cult.poly]) # Calculating SFS with a sub_sampling approach
print("save SFS for cultivated genotypes")
write.table(as.data.frame(table(sfs.cult.poly)),"/data2/projects/africrop_models/polarized_snps_mil/results_b/sfs.cult_CORRECTED_sampling.poly",quote = FALSE)

# Calculating group-specific SFS
print("calculating pop specifc SFSs")
sfs.bypop <- by(snps.lfmm[,snps.poly],snps.lfmm$Pop,FUN = sfs.daf4)
id <- names(sfs.bypop)
pop.size <- by(snps.lfmm$Pop,snps.lfmm$Pop,function(x){as.integer(0.75*length(x))})
write(sfs.bypop,"/data2/projects/africrop_models/polarized_snps_mil/results_b/sfs.bypop_Rformat_CORRECTED_sampling")
write.table(pop.size,"/data2/projects/africrop_models/polarized_snps_mil/results_b/pop.size_CORRECTED_sampling.txt")

# Transforming in a matrix format
print("transform SFS by pop in a matrix form more convenient to use")
sfs.bypop.b <- NULL
for(i in 1:length(sfs.bypop)){ sfs.bypop.b <- cbind(sfs.bypop.b,sfs.bypop[[i]])}
colnames(sfs.bypop.b) <- id
print("saving non filtered pop specific SFSs")
write(sfs.bypop.b,"/data2/projects/africrop_models/polarized_snps_mil/results_b/sfs.bypop_CORRECTED_sampling")

# Identifying monomorphic SNPs over the whole dataset and filtering them
lst <- which(apply(sfs.bypop.b,1,sum) == 0 |apply(sfs.bypop.b,1,sum) == 215 )
sfs.bypop.b <- sfs.bypop.b[-lst,]
print("saving monomorphic filtered pop specific SFSs")
write.table(sfs.bypop.b,"/data2/projects/africrop_models/polarized_snps_mil/results_b/sfs.bypop.poly_only_CORRECTED_sampling.txt")


# Generate pairwise 2D-SFS
print("generating SFS-2D")
for(i in 1:(ncol(sfs.bypop.b)-1)){
  for(j in (i+1):(ncol(sfs.bypop.b))){
toto <- table(factor(sfs.bypop.b[,i],levels = c(0:pop.size[[i]])),factor(sfs.bypop.b[,j],levels=c(0:pop.size[[j]])))
colnames(toto) <- paste("d",j-1,"_",colnames(toto),sep="")
row.names(toto) <- paste("d",i-1,"_",row.names(toto),sep="")

# Write 2D SFS to a file
print("saving SFS-2D")
conn <- file(paste("/data2/projects/africrop_models/polarized_snps_mil/results_b/Mil_jointDAFpop_CORRECTED_sampling",(j-1),"_",(i-1),".obs",sep=""),"w")
cat("1 observations\n",file=conn)
cat("\t",file=conn)
write.table(toto,conn,quote = F,sep = "\t")
close(conn)

# Plot the 2D-SFS
print("plotting SFS-2D")
{sfs2d.b <- melt(toto); names(sfs2d.b) <- c(names(pop.size[i]),names(pop.size[j]),"value")

sfs2d.plot <- ggplot(sfs2d.b, aes(sfs2d.b[,1],sfs2d.b[,2])) +
  geom_tile(aes(fill = log(value)), colour = "white") +
  scale_fill_continuous(low = "white", high = "black",na.value = "white") +
  theme(legend.position = "none",
        panel.border=element_blank(),
        axis.line=element_blank(),
        axis.text.x  = element_text(angle=270, vjust=0),
        plot.margin=unit(c(0,0,0,0),"mm"),
        panel.margin = unit(0,"null")) +
  labs(x = names(pop.size[i]),y = names(pop.size[j]))

hist1 <- ggplot() +
  geom_bar(aes(na.exclude(factor(sfs.bypop.b[,i],levels = c(0:pop.size[[i]])))),col="black",fill="gray75") +
  theme_bw() + theme(panel.grid = element_blank(),panel.border=element_blank(),
                     axis.title.x = element_blank(),axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     plot.margin=unit(c(0,0,0,0),"mm"),
                     panel.margin = unit(0,"null"))

hist2 <- ggplot() +
  geom_bar(aes(na.exclude(factor(sfs.bypop.b[,j],levels = c(0:pop.size[[j]])))),col="black",fill="gray75") +
  theme_bw() + theme(panel.grid = element_blank(),panel.border=element_blank(),
                     axis.title.y = element_blank(),axis.text.y = element_blank(),
                     axis.text.x  = element_text(angle=270, vjust=0),
                     axis.ticks.y = element_blank(),
                     plot.margin=unit(c(0,0,0,0),"mm"),
                     panel.margin = unit(0,"null")) +

  coord_flip()

p4 <- ggplot() +
  geom_blank() +
  theme(plot.margin=unit(c(0,0,0,0),"mm"))
g.sfs2d.plot <- ggplotGrob(sfs2d.plot)
g.hist1 <- ggplotGrob(hist1)
g.hist2 <- ggplotGrob(hist2)
g.p4 <- ggplotGrob(p4)

maxWidth = unit.pmax(g.sfs2d.plot$widths[2:3], g.hist2$widths[2:3],
                     g.hist1$widths[2:3])
g.sfs2d.plot$widths[2:3] <- maxWidth
g.hist2$widths[2:3] <- maxWidth
g.hist1$widths[2:3] <- maxWidth

pdf(file = paste("/data2/projects/africrop_models/polarized_snps_mil/results_b/Mil_jointDAFpop_",names(pop.size[j]),"_",names(pop.size[i]),"-with0_CORRECTED_sampling.pdf",sep=""))
print(plot_grid(g.hist1,g.p4,g.sfs2d.plot,g.hist2,ncol=2,align="hv",rel_widths = c(3,1),rel_heights = c(1,3)))
dev.off()

# The same without the 0 class
sfs2d.b <- melt(toto[-1,-1]); names(sfs2d.b) <- c(names(pop.size[i]),names(pop.size[j]),"value")

sfs2d.plot <- ggplot(sfs2d.b, aes(sfs2d.b[,1],sfs2d.b[,2])) +
  geom_tile(aes(fill = log(value)), colour = "white") +
  scale_fill_continuous(low = "white", high = "black",na.value = "white") +
  theme(legend.position = "none",
        panel.border=element_blank(),
        axis.line=element_blank(),
        axis.text.x  = element_text(angle=270, vjust=0),
        plot.margin=unit(c(0,0,0,0),"mm"),
        panel.margin = unit(0,"null")) +
  labs(x = names(pop.size[i]),y = names(pop.size[j]))

hist1 <- ggplot() +
  geom_bar(aes(na.exclude(factor(sfs.bypop.b[,i],levels = c(1:pop.size[[i]])))),col="black",fill="gray75") +
  theme_bw() + theme(panel.grid = element_blank(),panel.border=element_blank(),
                     axis.title.x = element_blank(),axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     plot.margin=unit(c(0,0,0,0),"mm"),
                     panel.margin = unit(0,"null"))

hist2 <- ggplot() +
  geom_bar(aes(na.exclude(factor(sfs.bypop.b[,j],levels = c(1:pop.size[[j]])))),col="black",fill="gray75") +
  theme_bw() + theme(panel.grid = element_blank(),panel.border=element_blank(),
                     axis.title.y = element_blank(),axis.text.y = element_blank(),
                     axis.text.x  = element_text(angle=270, vjust=0),
                     axis.ticks.y = element_blank(),
                     plot.margin=unit(c(0,0,0,0),"mm"),
                     panel.margin = unit(0,"null")) +

  coord_flip()

g.sfs2d.plot <- ggplotGrob(sfs2d.plot)
g.hist1 <- ggplotGrob(hist1)
g.hist2 <- ggplotGrob(hist2)
g.p4 <- ggplotGrob(p4)

maxWidth = unit.pmax(g.sfs2d.plot$widths[2:3], g.hist2$widths[2:3],
                     g.hist1$widths[2:3])
g.sfs2d.plot$widths[2:3] <- maxWidth
g.hist2$widths[2:3] <- maxWidth
g.hist1$widths[2:3] <- maxWidth
p <- plot_grid(g.hist1,g.p4,g.sfs2d.plot,g.hist2,ncol=2,align="hv",rel_widths = c(3,1),rel_heights = c(1,3))

pdf(file = paste("/data2/projects/africrop_models/polarized_snps_mil/results_b/Mil_jointDAFpop_",names(pop.size[i]),"_",names(pop.size[j]),"-without0_CORRECTED_sampling.pdf",sep=""),width =15,height=10)
print(p)
dev.off()}

  }
}

# Test for fastsimcoal : calculating SFS2D for cultivated and wild center
test <- (table(factor(sfs.bypop.n$cult_c,levels = c(0:25)),factor(sfs.bypop.n$wild_c,levels=c(0:9))))
colnames(test) <- paste("d0_",colnames(test),sep="")
row.names(test) <- paste("d1_",row.names(test),sep="")
test
conn <- file("/data2/projects/africrop_models/polarized_snps_mil/results_b/Mil_test_Wildc_Cultc_jointDAFpop1_0_CORRECTED_sampling.obs","w")
cat("1 observations\n",file=conn)
cat("\t",file=conn)
write.table(test,conn,quote = F,sep = "\t")
close(conn)
