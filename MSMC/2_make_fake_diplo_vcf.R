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

#### Sample 8 genotypes for msmc analysis, excluding eastern individuals and MR
glab <- read.table("list-glab-georef.txt",header=TRUE,na.strings = "-")
gen.list <- sample(glab$code_vcf[which(glab$code_vcf!="MR"&
                                        glab$code_vcf!="KC"&
                                        glab$code_vcf!="IK"&
                                        glab$code_vcf!="EK"&
                                        glab$code_vcf!="ND")],8)

#### Make a list of fake diplos to be created
diplos.list <- paste(gen.list[seq(1,length(gen.list),2)],gen.list[seq(2,length(gen.list),2)],sep="_")


#### define the file to be processed ####
for(c in 1:12){
setwd("msmc/VCF chr recoded/")
input <- paste("chr",c,".recode.HAPLO.NEW.vcf",sep="")
output <- paste("chr",c,"_",sep="")

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

# #### Defining the couples of genotypes to be merge in a fake diploid ####
# # from group 5
 g5_EA_FP <- c("EA","FP") ; g5_FT_GC <- c("FT","GC") ; g5_GP_HE <- c("GP","HE")
 g5_MI_NG <- c("MI","NG")
# 
# # from group 3
 g3_EI_FI <- c("EI","FI") ; g3_FL_GK <- c("FL","GK") ; g3_GR_HF <- c("GR","HF")
 g3_HM_KF <- c("HM","KF") ; g3_LA_MN <- c("LA","MN")
# 
# # from group 7
 g7_EN_FC<- c("EN","FC") ; g7_FS_GB <- c("FS","GB") ; g7_HB_HC <- c("HB","HC")
 g7_IM_IQ <- c("IM","IQ") ; g7_KM_LB <- c("KM","LB") ; g7_LR_MG <- c("LR","MG")
 g7_MN_MP <- c("MN","MP")
# 
# #### Defining the list of vcf to be created ####
 list <- c("g5_EA_FP", "g5_FT_GC", "g5_GP_HE", "g5_MI_NG", "g3_EI_FI", "g3_FL_GK", "g3_GR_HF", "g3_HM_KF", "g3_LA_MN",
 "g7_EN_FC", "g7_FS_GB", "g7_HB_HC", "g7_IM_IQ", "g7_KM_LB", "g7_LR_MG", "g7_MN_MP")

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
con.b = file(paste(output,colnames(temp)[10],"missing.vcf",sep=""),open="w")

# Writing header
writeLines(header[-length(header)],con = con)
writeLines(paste(names(temp),collapse="\t"),con = con)
writeLines(header[-length(header)],con = con.b)
writeLines(paste(names(temp),collapse="\t"),con = con.b)

# Writing only SNPs that fulfill the following conditions
for(l in 1:nrow(temp)){
if(temp[l,10]=="0|1"|temp[l,10]=="1|0"|temp[l,10]=="1|1"){ writeLines(paste(temp[l,],collapse="\t"),con = con)
  } else if(temp[l,10]==".|."){writeLines(paste(temp[l,],collapse="\t"),con = con.b) }  
} 
close(con)
close(con.b)
}

}

rm(list = ls())

#### Making the same for barthii

#### Sample 16 genotypes for msmc analysis
barth <- read.table("C:/Users/cubry/Documents/scripts/Riz/Riz/Riz_scripts_R/list-barth-georef.txt",header=TRUE,na.strings = "-")
gen.list <- sample(barth$code_vcf,16)

#### Make a list of fake diplos to be created
diplos.list <- paste(gen.list[seq(1,length(gen.list),2)],gen.list[seq(2,length(gen.list),2)],sep="_")


#### define the file to be processed ####
for(c in 1:12){
  setwd("C:/Users/cubry/Documents/Riz_analysis/msmc/VCF chr recoded/")
  input <- paste("chr",c,".barthii.recode.HAPLO.NEW.vcf",sep="")
  output <- paste("chr",c,"_barthii_",sep="")
  
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
    con.b = file(paste(output,colnames(temp)[10],"missing.vcf",sep=""),open="w")
    
    # Writing header
    writeLines(header[-length(header)],con = con)
    writeLines(paste(names(temp),collapse="\t"),con = con)
    writeLines(header[-length(header)],con = con.b)
    writeLines(paste(names(temp),collapse="\t"),con = con.b)
    
    # Writing only SNPs that fulfill the following conditions
    for(l in 1:nrow(temp)){
      if(temp[l,10]=="0|1"|temp[l,10]=="1|0"|temp[l,10]=="1|1"){ writeLines(paste(temp[l,],collapse="\t"),con = con)
      } else if(temp[l,10]==".|."){writeLines(paste(temp[l,],collapse="\t"),con = con.b) }  
    } 
    close(con)
    close(con.b)
  }
  
}


#### Computing fake diploids for barthii groups - ordered genotypes within groups ####
rm(list = ls())

#### Sample genotypes for msmc analysis

for(g in 1:5){
barth <- read.table(paste("C:/Users/cubry/Documents/Riz_analysis/msmc/Barthii-bestK6_for_PSMC_group",g,".txt",sep=""),header=FALSE,na.strings = "-",colClasses = "character")
gen.list <- NULL
for(i in seq(1,(nrow(barth)%/%2*2),2)){
  gen.list <- rbind(gen.list,paste(barth[i,],"_",barth[i+1,],sep=""))
}


#### Make a list of fake diplos to be created
diplos.list <- gen.list


#### define the file to be processed ####
for(c in 1:12){
  setwd("C:/Users/cubry/Documents/Riz_analysis/msmc/VCF chr recoded/")
  input <- paste("chr",c,".barthii.all.recode.HAPLO.NEW.vcf",sep="")
  output <- paste("chr",c,"_barthii_",sep="")
  
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
    con = file(paste(output,colnames(temp)[10],"_G",g,".vcf",sep=""),open="w")
    con.b = file(paste(output,colnames(temp)[10],"_G",g,"missing.vcf",sep=""),open="w")
    
    # Writing header
    writeLines(header[-length(header)],con = con)
    writeLines(paste(names(temp),collapse="\t"),con = con)
    writeLines(header[-length(header)],con = con.b)
    writeLines(paste(names(temp),collapse="\t"),con = con.b)
    
    # Writing only SNPs that fulfill the following conditions
    for(l in 1:nrow(temp)){
      if(temp[l,10]=="0|1"|temp[l,10]=="1|0"|temp[l,10]=="1|1"){ writeLines(paste(temp[l,],collapse="\t"),con = con)
      } else if(temp[l,10]==".|."){writeLines(paste(temp[l,],collapse="\t"),con = con.b) }  
    } 
    close(con)
    close(con.b)
  }
  
}
}



#### Computing fake diploids for barthii groups - using random genotypes within groups ####
rm(list = ls())

#### Sample genotypes for msmc analysis

for(g in 1:5){
  barth <- read.table(paste("C:/Users/cubry/Documents/Riz_analysis/msmc/Barthii-bestK6_for_PSMC_group",g,".txt",sep=""),header=FALSE,na.strings = "-",colClasses = "character")
  barth <- as.data.frame(sample(barth$V1,nrow(barth),FALSE))
  gen.list <- NULL
  for(i in seq(1,(nrow(barth)%/%2*2),2)){
    gen.list <- rbind(gen.list,paste(barth[i,],"_",barth[i+1,],sep=""))
  }
  
  #### Make a list of fake diplos to be created
  diplos.list <- gen.list
  
  
  #### define the file to be processed ####
  for(c in 1:12){
    setwd("C:/Users/cubry/Documents/Riz_analysis/msmc/VCF chr recoded/")
    input <- paste("chr",c,".barthii.all.recode.HAPLO.NEW.vcf",sep="")
    output <- paste("chr",c,"_barthii_",sep="")
    
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
      con = file(paste(output,colnames(temp)[10],"_G",g,".vcf",sep=""),open="w")
      con.b = file(paste(output,colnames(temp)[10],"_G",g,"missing.vcf",sep=""),open="w")
      
      # Writing header
      writeLines(header[-length(header)],con = con)
      writeLines(paste(names(temp),collapse="\t"),con = con)
      writeLines(header[-length(header)],con = con.b)
      writeLines(paste(names(temp),collapse="\t"),con = con.b)
      
      # Writing only SNPs that fulfill the following conditions
      for(l in 1:nrow(temp)){
        if(temp[l,10]=="0|1"|temp[l,10]=="1|0"|temp[l,10]=="1|1"){ writeLines(paste(temp[l,],collapse="\t"),con = con)
        } else if(temp[l,10]==".|."){writeLines(paste(temp[l,],collapse="\t"),con = con.b) }  
      } 
      close(con)
      close(con.b)
    }
    
  }
}