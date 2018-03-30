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
# Written by Philippe Cubry, Yves Vigouroux and Olivier François
#
###################################################################################################################################

# Load simulations from ABC retained ones after empirical data test
load("./Riz_analysis/splatche/nn10pop")

# Load group definition
groups <- read.table("./scripts/african_rice/ABC/og_groups_final.txt",na.strings = "NA.")
groups <- groups[groups$V1!="MR",]

# Randomly attribute 2000 singletons to 111 genotypes
groups$random <- as.data.frame(table(sample(111,2000,replace = T)))$Freq


#Retrieve some simulations
simuls <- cbind(
  nn10pop$unadj.values[row.names(nn10pop$unadj.values)%in%c("14105","16799","19834","19779","32696","53785","57960","64518","94649","100600"),],
  nn10pop$ss[row.names(nn10pop$ss)%in%c("14105","16799","19834","19779","32696","53785","57960","64518","94649","100600"),])

#Alterating SFS classes
simuls[,11:17] <- simuls[,11:17]*10000
simuls[,11] <- simuls[,11]+2000
simuls[,11:17] <- simuls[,11:17]/12000

#Modifying Rare variants classes
simuls[,18:30] *
sapply(1:13, FUN = function(i) length(which(groups$V2 == i)))


#
z.group = function(z){
  sapply(1:13, FUN = function(i) sum(z[groups$V2 == i]) )
}

#
z.groupm = function(z){
  sapply(1:13, FUN = function(i) mean(z[groups$V2 == i]) )
}

z <- groups$random
z.groupm(z)/sum(z.groupm(z))
