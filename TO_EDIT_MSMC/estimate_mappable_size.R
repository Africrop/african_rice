### Script to estimate mappable genome size based on masks generated for msmc

length <- 0

for(c in 1:12){
  # Get the mask file for each chromosome for the MG genotype
  input = paste("C:/Users/cubry/Documents/Riz_analysis/msmc/msmc_mappable_masks/chr",
              c,"/",list.files(path=paste("C:/Users/cubry/Documents/Riz_analysis/msmc/msmc_mappable_masks/chr",
                        c,"/",sep=""),pattern="MG"),sep="")
  con <- gzfile(input,open = "r")
  
  #Read the file
  while(length(oneline <- readLines(con,n=1, warn = FALSE))>0){
    oneline <- unlist(strsplit(oneline,split="\t"))
    tmp.length <- as.numeric(oneline[3])-as.numeric(oneline[2])+1
    length <- length + tmp.length
  }
  
close(con)  
  
  }


### Do the same for other genotypes for the sake of comparisons
##For "LR"
length.LR <- 0

for(c in 1:12){
  # Get the mask file for each chromosome for the LR genotype
  input = paste("C:/Users/cubry/Documents/Riz_analysis/msmc/msmc_mappable_masks/chr",
                c,"/",list.files(path=paste("C:/Users/cubry/Documents/Riz_analysis/msmc/msmc_mappable_masks/chr",
                                            c,"/",sep=""),pattern="LR"),sep="")
  con <- gzfile(input,open = "r")
  
  #Read the file
  while(length(oneline <- readLines(con,n=1, warn = FALSE))>0){
    oneline <- unlist(strsplit(oneline,split="\t"))
    tmp.length <- as.numeric(oneline[3])-as.numeric(oneline[2])+1
    length.LR <- length.LR + tmp.length
  }
  
  close(con)  
  
}
##For "GR"
length.GR <- 0

for(c in 1:12){
  # Get the mask file for each chromosome for the GR genotype
  input = paste("C:/Users/cubry/Documents/Riz_analysis/msmc/msmc_mappable_masks/chr",
                c,"/",list.files(path=paste("C:/Users/cubry/Documents/Riz_analysis/msmc/msmc_mappable_masks/chr",
                                            c,"/",sep=""),pattern="GR"),sep="")
  con <- gzfile(input,open = "r")
  
  #Read the file
  while(length(oneline <- readLines(con,n=1, warn = FALSE))>0){
    oneline <- unlist(strsplit(oneline,split="\t"))
    tmp.length <- as.numeric(oneline[3])-as.numeric(oneline[2])+1
    length.GR <- length.GR + tmp.length
  }
  
  close(con)  
  
}

##For "GB"
length.GB <- 0

for(c in 1:12){
  # Get the mask file for each chromosome for the GR genotype
  input = paste("C:/Users/cubry/Documents/Riz_analysis/msmc/msmc_mappable_masks/chr",
                c,"/",list.files(path=paste("C:/Users/cubry/Documents/Riz_analysis/msmc/msmc_mappable_masks/chr",
                                            c,"/",sep=""),pattern="GB"),sep="")
  con <- gzfile(input,open = "r")
  
  #Read the file
  while(length(oneline <- readLines(con,n=1, warn = FALSE))>0){
    oneline <- unlist(strsplit(oneline,split="\t"))
    tmp.length <- as.numeric(oneline[3])-as.numeric(oneline[2])+1
    length.GB <- length.GB + tmp.length
  }
  
  close(con)  
  
}



##### BARTHII #####

##For "BQ"
length.BQ <- 0

for(c in 1:12){
  # Get the mask file for each chromosome for the GR genotype
  input = paste("C:/Users/cubry/Documents/Riz_analysis/msmc/msmc_mappable_masks/chr",
                c,"/",list.files(path=paste("C:/Users/cubry/Documents/Riz_analysis/msmc/msmc_mappable_masks/chr",
                                            c,"/",sep=""),pattern="BQ"),sep="")
  con <- gzfile(input,open = "r")
  
  #Read the file
  while(length(oneline <- readLines(con,n=1, warn = FALSE))>0){
    oneline <- unlist(strsplit(oneline,split="\t"))
    tmp.length <- as.numeric(oneline[3])-as.numeric(oneline[2])+1
    length.BQ <- length.BQ + tmp.length
  }
  
  close(con)  
  
}

##For "CE"
length.CE <- 0

for(c in 1:12){
  # Get the mask file for each chromosome for the GR genotype
  input = paste("C:/Users/cubry/Documents/Riz_analysis/msmc/msmc_mappable_masks/chr",
                c,"/",list.files(path=paste("C:/Users/cubry/Documents/Riz_analysis/msmc/msmc_mappable_masks/chr",
                                            c,"/",sep=""),pattern="CE"),sep="")
  con <- gzfile(input,open = "r")
  
  #Read the file
  while(length(oneline <- readLines(con,n=1, warn = FALSE))>0){
    oneline <- unlist(strsplit(oneline,split="\t"))
    tmp.length <- as.numeric(oneline[3])-as.numeric(oneline[2])+1
    length.CE <- length.CE + tmp.length
  }
  
  close(con)  
  
}

##For "CV"
length.CV <- 0

for(c in 1:12){
  # Get the mask file for each chromosome for the GR genotype
  input = paste("C:/Users/cubry/Documents/Riz_analysis/msmc/msmc_mappable_masks/chr",
                c,"/",list.files(path=paste("C:/Users/cubry/Documents/Riz_analysis/msmc/msmc_mappable_masks/chr",
                                            c,"/",sep=""),pattern="CV"),sep="")
  con <- gzfile(input,open = "r")
  
  #Read the file
  while(length(oneline <- readLines(con,n=1, warn = FALSE))>0){
    oneline <- unlist(strsplit(oneline,split="\t"))
    tmp.length <- as.numeric(oneline[3])-as.numeric(oneline[2])+1
    length.CV <- length.CV + tmp.length
  }
  
  close(con)  
  
}

##For "DK"
length.DK <- 0

for(c in 1:12){
  # Get the mask file for each chromosome for the GR genotype
  input = paste("C:/Users/cubry/Documents/Riz_analysis/msmc/msmc_mappable_masks/chr",
                c,"/",list.files(path=paste("C:/Users/cubry/Documents/Riz_analysis/msmc/msmc_mappable_masks/chr",
                                            c,"/",sep=""),pattern="DK"),sep="")
  con <- gzfile(input,open = "r")
  
  #Read the file
  while(length(oneline <- readLines(con,n=1, warn = FALSE))>0){
    oneline <- unlist(strsplit(oneline,split="\t"))
    tmp.length <- as.numeric(oneline[3])-as.numeric(oneline[2])+1
    length.DK <- length.DK + tmp.length
  }
  
  close(con)  
  
}

##For "DN"
length.DN <- 0

for(c in 1:12){
  # Get the mask file for each chromosome for the GR genotype
  input = paste("C:/Users/cubry/Documents/Riz_analysis/msmc/msmc_mappable_masks/chr",
                c,"/",list.files(path=paste("C:/Users/cubry/Documents/Riz_analysis/msmc/msmc_mappable_masks/chr",
                                            c,"/",sep=""),pattern="DN"),sep="")
  con <- gzfile(input,open = "r")
  
  #Read the file
  while(length(oneline <- readLines(con,n=1, warn = FALSE))>0){
    oneline <- unlist(strsplit(oneline,split="\t"))
    tmp.length <- as.numeric(oneline[3])-as.numeric(oneline[2])+1
    length.DN <- length.DN + tmp.length
  }
  
  close(con)  
  
}