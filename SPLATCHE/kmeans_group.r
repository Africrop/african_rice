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
# Intellectual property belongs to  IRD and Grenoble-Alpes University
#
# Written by Olivier François and Philippe Cubry
#
###################################################################################################################################

####################
# Script to group genotypes based on their coordinates
# using a k-means procedure with a rejection algorithm
####################

# 


#' Takes a matrix of coordinates and output a group composition
#'
#' @param coord a matrix of coordinates
#' @param c the number of groups to build
#' @param m the minimal number of genotype in one group
#'
#' @return an object of class kmeans after completion of the rejection algorithm
#' @export
#'
#' @examples
kmeans_grouping = function(coord = coord, c = 10, m = 8) { grps = kmeans(coord, c) # Initial grouping, consider 10 groups (adjustable)
ta = table(grps$cl); tam = min(ta) # Compute the minimal number of ind per group
 while( tam < m){ #Loop to redo the kmeans until at least 8 individuals are presents in each group
   grps = kmeans(coord, c); ta = table(grps$cl); tam = min(ta)
 }
return(grps)
}