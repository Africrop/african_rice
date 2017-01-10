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

#######################################
# Script to tranform a vcf containing #
# only GT format and biallelic loci   #
# with several diploids to haploids   #
#######################################

#### define the file to be processed ####
print("initialising analysis")
#setwd("C:/Users/cubry/Documents/Riz_analysis/msmc/VCF chr recoded/")
setwd("/scratch/cubry/vcf_recoded/")
for(i in 1:12){
input <- paste("chr",i,".recode.vcf",sep="")
output <- paste("chr",i,".recode.HAPLO.NEW.vcf",sep="")
print(paste("process VCF",i,"of 12"))

#### Open the connection and determine the length (in row) of the header ####
print("open the connection")
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
print("reading data")
vcf_init <- read.table(file = input, header = FALSE,comment.char = "#",colClasses =  "character")
colnames(vcf_init) <- head
indiv <- head[10:length(head)]
print("extracting only GT format")
vcf_init <- apply(vcf_init,MARGIN=2,FUN=function(y){ y <- cbind(sapply(y, FUN=function(x){x <- unlist(strsplit(x,split=":"))[1]}))})

#### Making haploids from the diploid VCF using the following rules for replacement ####
print("Performing haploidization")
vcf_init[,10:ncol(vcf_init)][vcf_init[,10:ncol(vcf_init)] == "0/0"] <- 0
vcf_init[,10:ncol(vcf_init)][vcf_init[,10:ncol(vcf_init)] =="1/1"] <- 1
vcf_init[,10:ncol(vcf_init)][vcf_init[,10:ncol(vcf_init)] =="0/1"] <- sample(c(1,0),length(vcf_init[,10:ncol(vcf_init)][vcf_init[,10:ncol(vcf_init)] =="0/1"]),replace=TRUE)
vcf_init[,10:ncol(vcf_init)][vcf_init[,10:ncol(vcf_init)] =="./."] <- "."

#### Saving the resulting dataframe in a new vcf-like file ####
print("saving")
con <- file(output, open="w")
writeLines(text = header,con = con)
close(con)

write.table(vcf_init, file = output, append=TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE,sep="\t")
}


#### The same for barthii ####
#FOR GEOREFERENCED GENOTYPES ONLY
#### define the file to be processed ####
print("initialising analysis")
setwd("C:/Users/cubry/Documents/Riz_analysis/msmc/VCF chr recoded/")
#setwd("/scratch/cubry/vcf_recoded/")
for(i in 1:12){
  input <- paste("chr",i,".barthii.recode.vcf",sep="")
  output <- paste("chr",i,".barthii.recode.HAPLO.NEW.vcf",sep="")
  print(paste("process VCF",i,"of 12"))
  
  #### Open the connection and determine the length (in row) of the header ####
  print("open the connection")
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
  print("reading data")
  vcf_init <- read.table(file = input, header = FALSE,comment.char = "#",colClasses =  "character")
  colnames(vcf_init) <- head
  indiv <- head[10:length(head)]
  print("extracting only GT format")
  vcf_init <- apply(vcf_init,MARGIN=2,FUN=function(y){ y <- cbind(sapply(y, FUN=function(x){x <- unlist(strsplit(x,split=":"))[1]}))})
  
  #### Making haploids from the diploid VCF using the following rules for replacement ####
  print("Performing haploidization")
  vcf_init[,10:ncol(vcf_init)][vcf_init[,10:ncol(vcf_init)] == "0/0"] <- 0
  vcf_init[,10:ncol(vcf_init)][vcf_init[,10:ncol(vcf_init)] =="1/1"] <- 1
  vcf_init[,10:ncol(vcf_init)][vcf_init[,10:ncol(vcf_init)] =="0/1"] <- sample(c(1,0),length(vcf_init[,10:ncol(vcf_init)][vcf_init[,10:ncol(vcf_init)] =="0/1"]),replace=TRUE)
  vcf_init[,10:ncol(vcf_init)][vcf_init[,10:ncol(vcf_init)] =="./."] <- "."
  
  #### Saving the resulting dataframe in a new vcf-like file ####
  print("saving")
  con <- file(output, open="w")
  writeLines(text = header,con = con)
  close(con)
  
  write.table(vcf_init, file = output, append=TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE,sep="\t")
}


#FOR BARTHII GENOTYPES WITHOUT ADVENTICES
#### define the file to be processed ####
print("initialising analysis")
setwd("C:/Users/cubry/Documents/Riz_analysis/msmc/VCF chr recoded/")
#setwd("/scratch/cubry/vcf_recoded/")
for(i in 1:12){
  input <- paste("chr",i,".barthii.allwithtadv.recode.vcf",sep="")
  output <- paste("chr",i,".barthii.allwithtadv.recode.HAPLO.NEW.vcf",sep="")
  print(paste("process VCF",i,"of 12"))
  
  #### Open the connection and determine the length (in row) of the header ####
  print("open the connection")
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
  print("reading data")
  vcf_init <- read.table(file = input, header = FALSE,comment.char = "#",colClasses =  "character")
  colnames(vcf_init) <- head
  indiv <- head[10:length(head)]
  print("extracting only GT format")
  vcf_init <- apply(vcf_init,MARGIN=2,FUN=function(y){ y <- cbind(sapply(y, FUN=function(x){x <- unlist(strsplit(x,split=":"))[1]}))})
  
  #### Making haploids from the diploid VCF using the following rules for replacement ####
  print("Performing haploidization")
  vcf_init[,10:ncol(vcf_init)][vcf_init[,10:ncol(vcf_init)] == "0/0"] <- 0
  vcf_init[,10:ncol(vcf_init)][vcf_init[,10:ncol(vcf_init)] =="1/1"] <- 1
  vcf_init[,10:ncol(vcf_init)][vcf_init[,10:ncol(vcf_init)] =="0/1"] <- sample(c(1,0),length(vcf_init[,10:ncol(vcf_init)][vcf_init[,10:ncol(vcf_init)] =="0/1"]),replace=TRUE)
  vcf_init[,10:ncol(vcf_init)][vcf_init[,10:ncol(vcf_init)] =="./."] <- "."
  
  #### Saving the resulting dataframe in a new vcf-like file ####
  print("saving")
  con <- file(output, open="w")
  writeLines(text = header,con = con)
  close(con)
  
  write.table(vcf_init, file = output, append=TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE,sep="\t")
}




#FOR ALL BARTHII GENOTYPES
#### define the file to be processed ####
print("initialising analysis")
setwd("C:/Users/cubry/Documents/Riz_analysis/msmc/VCF chr recoded/")
#setwd("/scratch/cubry/vcf_recoded/")
for(i in 1:12){
  input <- paste("chr",i,".barthii.all.recode.vcf",sep="")
  output <- paste("chr",i,".barthii.all.recode.HAPLO.NEW.vcf",sep="")
  print(paste("process VCF",i,"of 12"))
  
  #### Open the connection and determine the length (in row) of the header ####
  print("open the connection")
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
  print("reading data")
  vcf_init <- read.table(file = input, header = FALSE,comment.char = "#",colClasses =  "character")
  colnames(vcf_init) <- head
  indiv <- head[10:length(head)]
  print("extracting only GT format")
  vcf_init <- apply(vcf_init,MARGIN=2,FUN=function(y){ y <- cbind(sapply(y, FUN=function(x){x <- unlist(strsplit(x,split=":"))[1]}))})
  
  #### Making haploids from the diploid VCF using the following rules for replacement ####
  print("Performing haploidization")
  vcf_init[,10:ncol(vcf_init)][vcf_init[,10:ncol(vcf_init)] == "0/0"] <- 0
  vcf_init[,10:ncol(vcf_init)][vcf_init[,10:ncol(vcf_init)] =="1/1"] <- 1
  vcf_init[,10:ncol(vcf_init)][vcf_init[,10:ncol(vcf_init)] =="0/1"] <- sample(c(1,0),length(vcf_init[,10:ncol(vcf_init)][vcf_init[,10:ncol(vcf_init)] =="0/1"]),replace=TRUE)
  vcf_init[,10:ncol(vcf_init)][vcf_init[,10:ncol(vcf_init)] =="./."] <- "."
  
  #### Saving the resulting dataframe in a new vcf-like file ####
  print("saving")
  con <- file(output, open="w")
  writeLines(text = header,con = con)
  close(con)
  
  write.table(vcf_init, file = output, append=TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE,sep="\t")
}