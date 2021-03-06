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

library(ggplot2)
library(reshape2)
library(gridExtra)
library(gtable)
library(cowplot)
library(LEA)
library(grid)
library(data.table)

#setwd("/data2/projects/africrop_models/Riz_fastsimcoal/")

# Load SNPs data from polarized biallelic SNPs file
snps.geno <- NULL
print("read dataset")
temp <- read.table("allChr.recode.HAPLO.NEW.vcf",header=T,skip=38,comment.char="",colClasses = "character")
names(temp)[1]<- "CHROM"
row.names(temp) <- paste(temp[,1],"_",temp[,2],sep="")
temp <- temp[,-c(1:9)]
temp[temp=="."]<- "9"
names(temp)[colnames(temp)=="ID1"] <- "ID"
names(temp)[colnames(temp)=="NA."] <- "NA"
snps.geno <- temp ; rm(temp)

# Load list of individuals
indiv <- read.table("100_ind_Mcswan.txt",na.strings="NaN"); names(indiv) <- c("Id","Species")

# Source used R scripts
source("functions_daf.R") # Load used functions

# reformat data
print("reformat data")
snps.lfmm <- t(snps.geno) ; rm(snps.geno)
snps.lfmm <- merge(snps.lfmm,indiv,by.x="row.names",by.y="Id",sort=F)
row.names(snps.lfmm) <- snps.lfmm$Row.names ; snps.lfmm$Row.names<-NULL
snps.lfmm.mat <- as.matrix(snps.lfmm[,-length(snps.lfmm)])

# Retaining only polymorphic loci over the set of 215 genotypes
print("filtering monomorphic loci over the whole dataset")
lst <- which(lapply(apply(snps.lfmm.mat[,],2,function(x){x=x[x!=9] ; unique(x)}),length)!=1) #List polymorphic loci
sfs <- apply(snps.lfmm.mat[,lst],2,daf.count)
print("save list of polymorphic loci")
write.table(names(lst),"listSNPs_recode.HAPLO.NEW.Mcswan.poly.txt")

# Calculating SFS on cultivated genotypes
print("calculating cultivated genotypes SFS")
lst.cult <- which(row.names(snps.lfmm.mat) %in% as.character(indiv[-grep(indiv[,2],pattern = "barth"),]$Id))
lst.cult.poly <- which(lapply(apply(snps.lfmm.mat[lst.cult,],2,function(x){x=x[x!=9] ; unique(x)}),length)!=1) #List polymorphic loci
sfs.cult.poly <- sfs.daf.count(snps.lfmm.mat[lst.cult,lst.cult.poly]) # Calculating SFS with a sub_sampling approach
print("save SFS for cultivated genotypes")
write.table(as.data.frame(table(sfs.cult.poly)),"allChr.recode.HAPLO.NEW.sfs.glabMcswan.poly",quote = FALSE)

# Calculating group-specific SFS
print("calculating pop specifc SFSs")
pop.size <- by(snps.lfmm$Species,snps.lfmm$Species,function(x){length(x)})
write.table(rbind(pop.size),"allChr.recode.HAPLO.NEW.popMcswan.size_CORRECTED.txt")
sfs.bypop <- by(snps.lfmm.mat[,lst],snps.lfmm$Species,FUN = sfs.daf.count)
id <- names(sfs.bypop)
save(sfs.bypop,file="allChr.recode.HAPLO.NEW.sfsMcswan.bypop_Rformat_CORRECTED")

# Transforming in a matrix format
print("transform SFS by pop in a matrix form more convenient to use")
sfs.bypop.b <- NULL
for(i in 1:length(sfs.bypop)){ sfs.bypop.b <- cbind(sfs.bypop.b,sfs.bypop[[i]])}
colnames(sfs.bypop.b) <- id
print("saving non filtered pop specific SFSs")
write.table(sfs.bypop.b,"allChr.recode.HAPLO.NEW.sfsMcswan.bypop_CORRECTED")

# Generate pairwise 2D-SFS
print("generating SFS-2D")
for(i in 1:(ncol(sfs.bypop.b)-1)){
  for(j in (i+1):(ncol(sfs.bypop.b))){
    toto <- table(factor(sfs.bypop.b[,i],levels = c(0:pop.size[[i]])),factor(sfs.bypop.b[,j],levels=c(0:pop.size[[j]])))
    colnames(toto) <- paste("d",j-1,"_",colnames(toto),sep="")
    row.names(toto) <- paste("d",i-1,"_",row.names(toto),sep="")
    
    # Write 2D SFS to a file
    print("saving SFS-2D")
    conn <- file(paste("allChr.recode.HAPLO.NEW._jointDAFpopMcswan_CORRECTED",(j-1),"_",(i-1),".obs",sep=""),"w")
    cat("1 observations\n",file=conn)
    cat("\t",file=conn)
    write.table(toto,conn,quote = F,sep = "\t")
    close(conn)
    
    # Plot the 2D-SFS
    print("plotting SFS-2D")
    
    sfs2d.b <- melt(toto)
    names(sfs2d.b) <- c(names(pop.size[i]),names(pop.size[j]),"value")
    write.table(sfs2d.b,"sfs2d.Mcswan.txt")
    
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
    
    pdf(file = paste("allChr.recode.HAPLO.NEW._jointDAFpopMcswan_",names(pop.size[j]),"_",names(pop.size[i]),"-with0_CORRECTED.pdf",sep=""))
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
    
    pdf(file = paste("allChr.recode.HAPLO.NEW._jointDAFpopMcswan_",names(pop.size[i]),"_",names(pop.size[j]),"-without0_CORRECTED.pdf",sep=""),width =15,height=10)
    print(p)
    dev.off()
    
  }
}
