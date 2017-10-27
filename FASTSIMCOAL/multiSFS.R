############
# Script to produce multi-dimensional SFS
# from a dataframe of named count data per population
############

library(ggplot2)
library(reshape2)
library(gridExtra)
library(gtable)
library(cowplot)
library(LEA)
library(grid)

# Read count data by pop for polarized mil SNPs
sfs.bypop <- read.table("/Users/cubry/Documents/mil_polarized_sfs/sfs.bypop.poly_only_CORRECTED.txt",header=TRUE)

# Getting population sizes
pop.size <- apply((na.exclude(sfs.bypop)),MARGIN = 2,max) ; names(pop.size)<- colnames(sfs.bypop)

# Producing multi-dimensional SFS, here 4D SFS with always three wild groups + one cultivated
# For fastsimcoal, the pop "0" will be the first called in factors list, pop "3" the last one

## cult center
multisfs.cult_c <- (as.data.frame(na.omit(table((sfs.bypop$wild_c),
                                                (sfs.bypop$wild_e),
                                                (sfs.bypop$wild_w),
                                                (sfs.bypop$cult_c),
                                                dnn = c("wild_c","wild_e","wild_w","cult_c")
))))

write.table(multisfs.cult_c,"/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs.cult_c")
write.table("1 observations. No. of demes and sample sizes are on next line","/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs_cult_c_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE)
write.table(t(c("4",pop.size["cult_c"],pop.size["wild_w"],pop.size["wild_e"],pop.size["wild_c"])),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs_cult_c_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)
write.table(t(multisfs.cult_c$Freq),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs_cult_c_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)


## cult east
multisfs.cult_e <- (as.data.frame(na.omit(table((sfs.bypop$wild_c),
                                                (sfs.bypop$wild_e),
                                                (sfs.bypop$wild_w),
                                                (sfs.bypop$cult_e),
                                                dnn = c("wild_c","wild_e","wild_w","cult_e")
))))
write.table(multisfs.cult_e,"/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs.cult_e")
write.table("1 observations. No. of demes and sample sizes are on next line","/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs_cult_e_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE)
write.table(t(c("4",pop.size["cult_e"],pop.size["wild_w"],pop.size["wild_e"],pop.size["wild_c"])),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs_cult_e_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)
write.table(t(multisfs.cult_e$Freq),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs_cult_e_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)


## cult west
multisfs.cult_w <- (as.data.frame(na.omit(table((sfs.bypop$wild_c),
                                                (sfs.bypop$wild_e),
                                                (sfs.bypop$wild_w),
                                                (sfs.bypop$cult_w),
                                                dnn = c("wild_c","wild_e","wild_w","cult_w")
))))
write.table(multisfs.cult_w,"/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs.cult_w")
write.table("1 observations. No. of demes and sample sizes are on next line","/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs_cult_w_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE)
write.table(t(c("4",pop.size["cult_w"],pop.size["wild_w"],pop.size["wild_e"],pop.size["wild_c"])),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs_cult_w_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)
write.table(t(multisfs.cult_w$Freq),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs_cult_w_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)


## cult india
multisfs.cult_i <- (as.data.frame(na.omit(table((sfs.bypop$wild_c),
                                                (sfs.bypop$wild_e),
                                                (sfs.bypop$wild_w),
                                                (sfs.bypop$cult_i),
                                                dnn = c("wild_c","wild_e","wild_w","cult_i")
))))
write.table(multisfs.cult_i,"/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs.cult_i")
write.table("1 observations. No. of demes and sample sizes are on next line","/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs_cult_i_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE)
write.table(t(c("4",pop.size["cult_i"],pop.size["wild_w"],pop.size["wild_e"],pop.size["wild_c"])),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs_cult_i_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)
write.table(t(multisfs.cult_i$Freq),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs_cult_i_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)

## cult south
multisfs.cult_s <- (as.data.frame(na.omit(table((sfs.bypop$wild_c),
                                                (sfs.bypop$wild_e),
                                                (sfs.bypop$wild_w),
                                                (sfs.bypop$cult_s),
                                                dnn = c("wild_c","wild_e","wild_w","cult_s")
))))
write.table(multisfs.cult_s,"/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs.cult_s")
write.table("1 observations. No. of demes and sample sizes are on next line","/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs_cult_s_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE)
write.table(t(c("4",pop.size["cult_s"],pop.size["wild_w"],pop.size["wild_e"],pop.size["wild_c"])),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs_cult_s_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)
write.table(t(multisfs.cult_s$Freq),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs_cult_s_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)

# Producing a 3D-SFS for wild groups only
multisfs.wild <- (as.data.frame(na.omit(table((sfs.bypop$wild_c),
                                                (sfs.bypop$wild_e),
                                                (sfs.bypop$wild_w),
                                                dnn = c("wild_c","wild_e","wild_w")
))))
write.table(multisfs.wild,"/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs.wild")
write.table("1 observations. No. of demes and sample sizes are on next line","/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs_wild_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE)
write.table(t(c("3",pop.size["wild_w"],pop.size["wild_e"],pop.size["wild_c"])),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs_wild_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)
write.table(t(multisfs.wild$Freq),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs_wild_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)


# Producing a 3D-SFS to investigate indian origin
multisfs.cult <- (as.data.frame(na.omit(table((sfs.bypop$cult_e),
                                              (sfs.bypop$cult_s),
                                              (sfs.bypop$cult_i),
                                              dnn = c("cult_e","cult_s","cult_i")
))))
write.table(multisfs.cult,"/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs.cult_esi")
write.table("1 observations. No. of demes and sample sizes are on next line","/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs_cult_esi_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE)
write.table(t(c("3",pop.size["cult_e"],pop.size["cult_s"],pop.size["cult_i"])),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs_cult_esi_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)
write.table(t(multisfs.cult$Freq),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS/multisfs_cult_esi_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)
