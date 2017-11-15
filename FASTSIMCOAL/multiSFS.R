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
sfs.bypop <- read.table("allChr.recode.HAPLO.NEW.sfs.bypop_CORRECTED.withmonomorphs",header=TRUE)

# Getting population sizes
pop.size <- apply((na.exclude(sfs.bypop)),MARGIN = 2,max) ; names(pop.size)<- colnames(sfs.bypop)

# Producing multi-dimensional SFS, here 2D
# For fastsimcoal, the pop "0" will be the first called in factors list
multisfs.rice <- (as.data.frame(na.omit(table((sfs.bypop$barthii),
                                              (sfs.bypop$glab),
                                              dnn = c("barthii","glab")
))))

write.table(multisfs.rice,"/data2/projects/africrop_models/Riz_fastsimcoal/multisfs_rice_DSFS.obs")
write.table("1 observations. No. of demes and sample sizes are on next line","/data2/projects/africrop_models/Riz_fastsimcoal/multisfs_rice_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE)
write.table(t(c("2",pop.size["glab"],pop.size["barthii"])),"/data2/projects/africrop_models/Riz_fastsimcoal/multisfs_rice_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)
write.table(t(multisfs.rice),"/data2/projects/africrop_models/Riz_fastsimcoal/multisfs_rice_DSFS.obs", col.names = FALSE, row.names = FALSE,quote=FALSE,append=TRUE)