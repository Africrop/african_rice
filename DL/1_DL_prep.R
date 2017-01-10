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


###################################################################################################################################
#
# This script calculate mean r2 values for each position from vcftools generated files
# Files produced are used to feed DL_plotting.R
#
###################################################################################################################################


#### Perform analysis for glaberrima ####
# Initiate variable
glab.DL <- NULL

# Read the vcftools produced table of r2 values
tmp <- read.table("glab_all.recode.hap.ld.hap.ld",header=TRUE)

# Calculate distance between considered polymorphisms
tmp$dist <- tmp$POS2 - tmp$POS1
glab.DL <- cbind(tmp$R.2,tmp$dist)
colnames(glab.DL) <- c("R.2","dist")

# Save the resulting file
write.table(glab.DL, "glab.R.2.txt")

# Calculate mean r2 for each considered distance
glab.DL <- as.data.frame(glab.DL)
glab.mean.R.2 <- NULL
glab.mean.R2 <- cbind(as.data.frame(tapply(glab.DL$R.2, glab.DL$dist,mean)), as.data.frame(tapply(tmp$R.2, glab.DL$dist,var)) ,as.data.frame(table(glab.DL$dist)))
# Save the resulting file
write.table(file="glab.mean.R.2.txt",x= glab.mean.R2)


#### Perform analysis for glaberrima ####
# Initiate variable
barth.DL <- NULL

# Read the vcftools produced table of r2 values
tmp <- read.table("barth_all.recode.hap.ld.hap.ld",header=TRUE)

# Calculate distance between considered polymorphisms
tmp$dist <- tmp$POS2 - tmp$POS1
barth.DL <- cbind(tmp$R.2,tmp$dist)
colnames(barth.DL) <- c("R.2","dist")
barth.DL <- as.data.frame(barth.DL)

# Save the resulting file
con <- file("barth.R.2.NEW.txt",open="wt")
for(i in nrow(barth.DL)){
  writeLines(con, barth.DL[i,])
}
close(con)

# Calculate mean r2 for each considered distance
barth.mean.R2 <- NULL
barth.mean.R2 <- cbind(as.data.frame(tapply(barth.DL $R.2, barth.DL $dist,mean)), as.data.frame(tapply(barth.DL $R.2, barth.DL $dist,var)),as.data.frame(table(barth.DL $dist)))

# Save the resulting file
con <- file("barthii.mean.R.2.new.txt",open="wt")
for(i in nrow(barth.mean.R2)){
  writeLines(con, barth.mean.R2[i,])
}
close(con)
