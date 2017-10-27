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
sfs.bypop_sample <- read.table("/Users/cubry/Documents/mil_polarized_sfs/sfs.bypop.poly_only_CORRECTED_sampling.txt",header=TRUE)

# Getting population sizes
pop.size_sample <- apply((na.exclude(sfs.bypop_sample)),MARGIN = 2,max) ; names(pop.size_sample)<- colnames(sfs.bypop_sample)

# Producing multi-dimensional SFS, here 4D SFS with always three wild groups + one cultivated
# For fastsimcoal, the pop "0" will be the first called in factors list, pop "3" the last one

## cult center
multisfs.cult_c_sample <- (as.data.frame(na.omit(table((sfs.bypop_sample$wild_c),
                                                (sfs.bypop_sample$wild_e),
                                                (sfs.bypop_sample$wild_w),
                                                (sfs.bypop_sample$cult_c),
                                                dnn = c("wild_c","wild_e","wild_w","cult_c"))
)))

write.table(multisfs.cult_c_sample,"/Users/cubry/Documents/mil_polarized_sfs/multiSFS_sample/multisfs.cult_c_sample")
write.table("1 observations. No. of demes and sample sizes are on next line","/Users/cubry/Documents/mil_polarized_sfs/multiSFS_sample/multisfs_cult_c_sample_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE)
write.table(t(c("4",pop.size_sample["cult_c"],pop.size_sample["wild_w"],pop.size_sample["wild_e"],pop.size_sample["wild_c"])),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS_sample/multisfs_cult_c_sample_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)
write.table(t(multisfs.cult_c_sample$Freq),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS_sample/multisfs_cult_c_sample_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)

## cult east
multisfs.cult_e_sample <- (as.data.frame(na.omit(table((sfs.bypop_sample$wild_c),
                                                (sfs.bypop_sample$wild_e),
                                                (sfs.bypop_sample$wild_w),
                                                (sfs.bypop_sample$cult_e),
                                                dnn = c("wild_c","wild_e","wild_w","cult_e")
))))
write.table(multisfs.cult_e_sample,"/Users/cubry/Documents/mil_polarized_sfs/multiSFS_sample/multisfs.cult_e_sample")
write.table("1 observations. No. of demes and sample sizes are on next line","/Users/cubry/Documents/mil_polarized_sfs/multiSFS_sample/multisfs_cult_e_sample_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE)
write.table(t(c("4",pop.size_sample["cult_e"],pop.size_sample["wild_w"],pop.size_sample["wild_e"],pop.size_sample["wild_c"])),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS_sample/multisfs_cult_e_sample_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)
write.table(t(multisfs.cult_e_sample$Freq),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS_sample/multisfs_cult_e_sample_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)

## cult west
multisfs.cult_w_sample <- (as.data.frame(na.omit(table((sfs.bypop_sample$wild_c),
                                                (sfs.bypop_sample$wild_e),
                                                (sfs.bypop_sample$wild_w),
                                                (sfs.bypop_sample$cult_w),
                                                dnn = c("wild_c","wild_e","wild_w","cult_w")
))))
write.table(multisfs.cult_w_sample,"/Users/cubry/Documents/mil_polarized_sfs/multiSFS_sample/multisfs.cult_w_sample")
write.table("1 observations. No. of demes and sample sizes are on next line","/Users/cubry/Documents/mil_polarized_sfs/multiSFS_sample/multisfs_cult_w_sample_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE)
write.table(t(c("4",pop.size_sample["cult_w"],pop.size_sample["wild_w"],pop.size_sample["wild_e"],pop.size_sample["wild_c"])),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS_sample/multisfs_cult_w_sample_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)
write.table(t(multisfs.cult_w_sample$Freq),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS_sample/multisfs_cult_w_sample_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)

## cult india
multisfs.cult_i_sample <- (as.data.frame(na.omit(table((sfs.bypop_sample$wild_c),
                                                (sfs.bypop_sample$wild_e),
                                                (sfs.bypop_sample$wild_w),
                                                (sfs.bypop_sample$cult_i),
                                                dnn = c("wild_c","wild_e","wild_w","cult_i")
))))
write.table(multisfs.cult_i_sample,"/Users/cubry/Documents/mil_polarized_sfs/multiSFS_sample/multisfs.cult_i_sample")
write.table("1 observations. No. of demes and sample sizes are on next line","/Users/cubry/Documents/mil_polarized_sfs/multiSFS_sample/multisfs_cult_i_sample_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE)
write.table(t(c("4",pop.size_sample["cult_i"],pop.size_sample["wild_w"],pop.size_sample["wild_e"],pop.size_sample["wild_c"])),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS_sample/multisfs_cult_i_sample_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)
write.table(t(multisfs.cult_i_sample$Freq),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS_sample/multisfs_cult_i_sample_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)

## cult south
multisfs.cult_s_sample <- (as.data.frame(na.omit(table((sfs.bypop_sample$wild_c),
                                                (sfs.bypop_sample$wild_e),
                                                (sfs.bypop_sample$wild_w),
                                                (sfs.bypop_sample$cult_s),
                                                dnn = c("wild_c","wild_e","wild_w","cult_s")
))))
write.table(multisfs.cult_s_sample,"/Users/cubry/Documents/mil_polarized_sfs/multiSFS_sample/multisfs.cult_s_sample")
write.table("1 observations. No. of demes and sample sizes are on next line","/Users/cubry/Documents/mil_polarized_sfs/multiSFS_sample/multisfs_cult_s_sample_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE)
write.table(t(c("4",pop.size_sample["cult_s"],pop.size_sample["wild_w"],pop.size_sample["wild_e"],pop.size_sample["wild_c"])),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS_sample/multisfs_cult_s_sample_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)
write.table(t(multisfs.cult_s_sample$Freq),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS_sample/multisfs_cult_s_sample_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)

# Producing a 3D-SFS for wild groups only
multisfs.wild_sample <- (as.data.frame(na.omit(table((sfs.bypop_sample$wild_c),
                                                (sfs.bypop_sample$wild_e),
                                                (sfs.bypop_sample$wild_w),
                                                dnn = c("wild_c","wild_e","wild_w")
))))
write.table(multisfs.wild_sample,"/Users/cubry/Documents/mil_polarized_sfs/multiSFS_sample/multisfs_sample.wild")
write.table("1 observations. No. of demes and sample sizes are on next line","/Users/cubry/Documents/mil_polarized_sfs/multiSFS_sample/multisfs_wild_sample_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE)
write.table(t(c("3",pop.size_sample["wild_w"],pop.size_sample["wild_e"],pop.size_sample["wild_c"])),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS_sample/multisfs_wild_sample_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)
write.table(t(multisfs.wild_sample$Freq),"/Users/cubry/Documents/mil_polarized_sfs/multiSFS_sample/multisfs_wild_sample_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)
