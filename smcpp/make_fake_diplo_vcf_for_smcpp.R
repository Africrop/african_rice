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
# Written by Philippe Cubry
#
###################################################################################################################################

#######################################
# Script to tranform a vcf-like file  #
# containing haploid individuals in   #
# one fake diploid in a new VCF file  #
#######################################

#### Sample genotypes for smc++ analysis, excluding eastern individuals and MR
glab <- read.table("Rice_smcpp/Liste-individus-riz-vcf-08-04-2016_glaberima.txt",
                   header=FALSE,na.strings = "-",
                   colClasses = "character")[,1]
gen.list <- sample(glab[which(glab!="MR"&
                                        glab!="KC"&
                                        glab!="IK"&
                                        glab!="EK"&
                                        glab!="ND")])

#### Make a list of fake diplos to be created
diplos.list <- paste(gen.list[seq(1,length(gen.list),2)],gen.list[seq(2,length(gen.list),2)],sep="_")
write.table(diplos.list,file = "Rice_smcpp/list_of_pseudo_diploids_Og.txt")

#### define the file to be processed ####
for(c in 1:12){
input <- paste("Riz_vcf/chr",c,".recode.HAPLO.NEW.vcf.gz",sep="")
output <- paste("Rice_smcpp/chr",c,"_",sep="")

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
indiv <- head[10:length(head)]

#### Making the diploidization with replacement rules ####
for (i in 1:length(diplos.list)){
# Gather ind data
    a <- unlist(strsplit(diplos.list[i],split="_"))
# Pasting in a phased fake diploid
    temp <- paste(vcf_haplo[,a[1]],vcf_haplo[,a[2]],sep="|")
# Replacing missing data
    temp[temp == "0|."| temp == ".|0"| temp == "1|."|
       temp == ".|1"| temp == ".|."] <- ".|."

#### Generating a vcf
    temp <- cbind(vcf_haplo[1:9], temp,stringsAsFactors=FALSE) ; colnames(temp)[10] <- diplos.list[i]
assign(paste("vcf_",colnames(temp)[10],sep=""),temp)

# Saving the vcf file
con = file(paste(output,colnames(temp)[10],".vcf",sep=""),open="w")

# Writing header
writeLines(header[-length(header)],con = con)
writeLines(paste(names(temp),collapse="\t"),con = con)

# Writing only SNPs that fulfill the following conditions
for(l in 1:nrow(temp)){
 writeLines(paste(temp[l,],collapse="\t"),con = con)

} 
close(con)

}

}

rm(list = ls())

#### Making the same for barthii

#### Sample genotypes for smc++ analysis
barth <- read.table("Rice_smcpp/Liste-individus-riz-vcf-08-04-2016_barthii.txt",
                    header=FALSE,na.strings = "-",
                    colClasses = "character")[,1]
gen.list <- sample(barth)

#### Make a list of fake diplos to be created
diplos.list <- paste(gen.list[seq(1,length(gen.list),2)],gen.list[seq(2,length(gen.list),2)],sep="_")
write.table(diplos.list,file = "Rice_smcpp/list_of_pseudo_diploids_Ob.txt")


#### define the file to be processed ####
for(c in 1:12){
  input <- paste("Riz_vcf/chr",c,".barthii.all.recode.HAPLO.NEW.vcf.gz",sep="")
  output <- paste("Rice_smcpp/chr",c,"_barthii_",sep="")
  
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
  indiv <- head[10:length(head)]
  
  #### Making the diploidization with replacement rules ####
  for (i in 1:length(diplos.list)){
    # Gather ind data
    a <- unlist(strsplit(diplos.list[i],split="_"))
    # Pasting in a phased fake diploid
    temp <- paste(vcf_haplo[,a[1]],vcf_haplo[,a[2]],sep="|")
    # Replacing missing data
    temp[temp == "0|."| temp == ".|0"| temp == "1|."|
           temp == ".|1"| temp == ".|."] <- ".|."
    
    #### Generating a vcf
    temp <- cbind(vcf_haplo[1:9], temp,stringsAsFactors=FALSE) ; colnames(temp)[10] <- diplos.list[i]
    assign(paste("vcf_",colnames(temp)[10],sep=""),temp)
    
    # Saving the vcf file
    con = file(paste(output,colnames(temp)[10],".vcf",sep=""),open="w")

    # Writing header
    writeLines(header[-length(header)],con = con)
    writeLines(paste(names(temp),collapse="\t"),con = con)

    # Writing only SNPs that fulfill the following conditions
    for(l in 1:nrow(temp)){
       writeLines(paste(temp[l,],collapse="\t"),con = con)
      } 
    close(con)
    
  }
  
}
