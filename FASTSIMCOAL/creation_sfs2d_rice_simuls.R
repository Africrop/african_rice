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

#######################################################
# Script to generate 2D-SFS from ms simulation data of#
# 100 inidividuals (50 wild + 50 cult)                #
# Requirement: a script to source with functions used #
# - a file with SNP info                              #
# CAUTION : Version to estimate unfolded SFS based on #
# polarised SNPs                                      #
#######################################################

library(ggplot2)
library(reshape2)
library(gridExtra)
library(gtable)
library(cowplot)
library(LEA)
library(grid)

# Load SNPs data from polarized biallelic SNPs file
load('./FASTSIMCOAL/out_ms_lfmm')

# Source used R scripts
source("./FASTSIMCOAL/functions_daf.R") # Load used functions

# Calculating SFS on wild genotypes
print("calculating wild genotypes SFS")
sfs.wild <- sfs.daf.count(datlf[1:50,])
plot(table(sfs.wild))

# Calculating SFS on cultivated genotypes
print("calculating cultivated genotypes SFS")
sfs.cult <- sfs.daf.count(datlf[51:100,])
plot(table(sfs.cult))

# Calculating group-specific SFS
datlf <- as.data.frame(datlf)
datlf$Pop <- as.vector(c(rep("wild",50),rep("cult",50)))
print("calculating pop specifc SFSs")
sfs.bypop <- by(datlf[,-ncol(datlf)],as.factor(datlf$Pop),FUN = sfs.daf.count)
id <- names(sfs.bypop)
pop.size <- c(50,50)

# Transforming in a matrix format
print("transform SFS by pop in a matrix form more convenient to use")
sfs.bypop.b <- NULL
for(i in 1:length(sfs.bypop)){ sfs.bypop.b <- cbind(sfs.bypop.b,sfs.bypop[[i]])}
colnames(sfs.bypop.b) <- id

# Generate pairwise 2D-SFS
print("generating SFS-2D")
for(i in 1:(ncol(sfs.bypop.b)-1)){
  for(j in (i+1):(ncol(sfs.bypop.b))){
toto <- table(factor(sfs.bypop.b[,i],levels = c(0:pop.size[[i]])),factor(sfs.bypop.b[,j],levels=c(0:pop.size[[j]])))
colnames(toto) <- paste("d",j-1,"_",colnames(toto),sep="")
row.names(toto) <- paste("d",i-1,"_",row.names(toto),sep="")

  }}

# Plot the 2D-SFS
print("plotting SFS-2D")
sfs2d.b <- melt(toto)

sfs2d.plot <- ggplot(sfs2d.b, aes(sfs2d.b[,1],sfs2d.b[,2])) +
  geom_tile(aes(fill=log(value)),col="white" ) +
  scale_fill_continuous(low = "white", high = "black",na.value = "white") +
  theme(legend.position = "none",
        panel.border=element_blank(),
        axis.line=element_blank(),
        axis.text.x  = element_text(angle=270, vjust=0),
        plot.margin=unit(c(0,0,0,0),"mm"),
        panel.margin = unit(0,"null")) +
  labs(x = names(pop.size[i]),y = names(pop.size[j]))
sfs2d.plot

# The same without the 0 class
sfs2d.b <- melt(toto[-1,-1])

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
sfs2d.plot