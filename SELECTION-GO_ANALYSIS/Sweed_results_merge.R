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
# Written by Concetta Burgarella
#
###################################################################################################################################

##########################################################################################################################
### Concetta Burgarella
###
### January - Mars 2016   
###
### Work on SweeD results (detection of completed selctive sweeps - Pavlidis et al. 2012)
### There is a file per chromosome, per cultivated group and per grid size (i.e. nb of genomic position associated to a CLR)
###   
### Steps of the script:
###
### 1) Read the results per chromosome and merge them for each group analysed.  
###    Retain only the 0.1% highest CLR values per group, ie. the same nb of retained positions for all groups,            
###    but different CLR threshold per group. Put them in a csv file. 
###    
### 2) Merge selected positions if pairwise distances < 100 Kb, to identify the genomic regions/fragments with signatures of 
###    selection across groups (here selected positions are not separated per group). Put them in a csv file.
###
### 3) Venn diagram of "identified" fragments (plot also fragment size distribution)   
###    
###########################################################################################################################




### 1) Read the results per chromosome and merge them for each cultivated group.  
###    Retain only the 0.1% highest CLR values per group, ie. the same nb of retained positions for all groups,    
###    but different CLR threshold per group.

setwd("/.../.../Pathway/.../Folder/")

# list folders (one folder per analysed group)
dirs <- list.dirs(path = '.' , full.names = F, recursive = TRUE)

# loop on folders (ie. apply CLR threshold per group after merging individual chrs)
for(dir in dirs[-1]){
  print(dir)
  
  # select result files per group folder
  files1 <- list.files(path = dir, pattern = "Report")         # files with "Report" in its name
  files2 <- Filter(function(x) grepl("10", x), files1)        # files with "10" in its name (ie. 10e4 grid)
  
  print(files2)
  
  
  # loop on files within group folder
  
   for(file in files2){
     print(file)
     # read files one by one
     content <- read.table(paste(dir,"/",file,sep=""), skip=3)
 
     names(content) <- c("position","CLR","alpha")
     # add column info on chr and group
     chr <- substr(file, 18,22)
     temp <- cbind(chr,content,dir)
     # merge files within group
     if(file==files2[1]){
       merged_table <- temp
     }else{
       merged_table <- rbind(merged_table,temp)
     }
     
   }   
     
  # select 1% highest CLR values
  CLR_thres <- min(tail(sort(merged_table$CLR), nrow(merged_table)/100))
  useful_table <- merged_table[which(merged_table$CLR >= CLR_thres),]
  
  # merge files for all groups
  if(dir==dirs[2]){
    useful_table_all <- useful_table
  }else{
    useful_table_all <- rbind(useful_table_all,useful_table)
  }
  
  write.csv(useful_table_all,"SweeD_results_CLR_thres_0.01.csv")
  
}






### 2) Merge selected positions if pairwise distances < 100 Kb, to identify the genomic regions/fragments with signatures of 
###    selection across groups (here selected positions are not separated per group)

# if starting from this point:
useful_table_all <- read.csv("/.../Pathway/.../SweeD_results_CLR_thres_0.01.csv")



frag_table <- data.frame()  # creates df to store results
table_venn <- data.frame()  # creates df to store results

# plot the distribution of fragment size per cromosome and over all chrs in a separate window 
X11()
op <- par(mfrow=c(2,4))

# work per chr
for(chr in c("CHR1_","CHR2_","CHR3_","CHR4_","CHR5_","CHR6_","CHR7_","CHR8_","CHR9_","CHR10","CHR11","CHR12")){
 
  temp <- useful_table_all[which(useful_table_all$chr == chr),]

  if(nrow(temp)>0){
  
  # sort row by position within chr
  temp2 <- temp[order(temp$position),]

  # set fragment delimitation (ie within 100000 max distance)
  fragment <- array("end", nrow(temp2))  # creates the array to stock the fragment ID
  
  # for each position checks if the distance to the following one is < 100kb, 
  # in this case assign the two positions to the same fragment
  for(i in 1:(nrow(temp2)-1)){
    diff <- temp2$position[i+1] - temp2$position[i]
    if(diff < 100000){
      fragment[i] <- 1
    }
  }

  temp3 <- cbind(temp2,fragment)

  # Record the starting and ending position of fragment interval
  end   <- temp3$position[which(temp3$fragment =="end")]
  start1 <- temp3$position[1]
  start <- head(n=length(end),c(start1, temp3$position[which(temp3$fragment =="end")+1]))
  
  fragment_lenght <- end-start
  hist(fragment_lenght, breaks=20,main = chr)

  # Sort out what cultivated groups are interested by each fragment (i.e. had CLR positions within the defined fragment)
  frag_table_temp <- data.frame()
  for(i in 1:length(start)){
    temp4 <- temp3[intersect(which(temp3$position>=start[i]), which(temp3$position<=end[i])),]
    groups <- as.character(unique(temp4$dir)); print(c(chr,i,start[i],groups))
    frag_lines <- cbind(chr, i, start[i], end[i],groups)  # one line per fragment and group
    frag_table_temp <- rbind(frag_table_temp,frag_lines)
  }
  names(frag_table_temp) <- c("chr","frag_id_x_chr","start_pos","end_pos","group")
  
  # change to numeric the two columns of factors
  frag_table_temp[,3:4] <- apply(frag_table_temp[,3:4], MARGIN=2, FUN= function(x) as.numeric(as.character(x)))
  
  # increase the fragment length up to 100 Kb if it is < 100Kb
  for(i in 1:nrow(frag_table_temp)){
    
    # calculate fragment length
    frag_length <- frag_table_temp$end_pos[i] - frag_table_temp$start_pos[i]
    if(frag_length < 100000){
      print(c(i, frag_length))
      # calculate the bp increment to get 100kb
      increment <- 100000 - frag_length
      print(c("increment",increment))
      
      # set the new start and end position of SD fragment
      frag_table_temp$start_pos[i] <- frag_table_temp$start_pos[i] - increment/2
      frag_table_temp$end_pos[i] <- frag_table_temp$end_pos[i] + increment/2
    }
  }
  
  # fill tables
  table_venn_chr <- as.matrix(table(frag_table_temp$start_pos,frag_table_temp$group))
  table_venn <- rbind(table_venn, table_venn_chr) 

  frag_table <- rbind(frag_table,frag_table_temp)
  } # end of "if"
}  # end of the loop on chrs

names(frag_table) <- c("chr","frag_id_x_chr","start_pos","end_pos","group")

# plot fragment size across chromosomes
png(filename = "frgmt-size-distribution.png", width = 20, height = 15, units = "cm", res = 300)
size <- as.numeric(as.character(frag_table$end_pos)) - as.numeric(as.character(frag_table$start_pos))
hist(log10(size), main = "all chromosomes")
summary(size)
dev.off()


write.csv(frag_table,"SweeD_fragments_CLR_thres_0.01.csv")








### 3) Venn diagram 

#source("https://bioconductor.org/biocLite.R")

#biocLite("limma")

library(limma)

# if starting from this point
# frag_table <- read.csv("C:/.../Pathway/.../SweeD_fragments_cult_CLR_thres_0.001.csv")




table_venn <- table(frag_table$start_pos,frag_table$group)


a <- vennCounts(table_venn, include="both")
X11()
vennDiagram(a, include = "both", names=c("barthii", "glaberrima"), circle.col=c("green","red"),cex=1,lwd=5,show.include=F)

png(filename = "venn-diagram.png", width = 20, height = 15, units = "cm", res = 300)
vennDiagram(a, include = "both", names=c("barthii", "glaberrima"), circle.col=c("green","red"),cex=1,lwd=5,show.include=F)
dev.off()









