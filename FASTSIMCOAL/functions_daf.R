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

### Set of functions to use to estimate SFS
### Warning : only usefull with O/1 coded lfmm-style matrix
### Adaptation might be needed to use with other kind of data
### Especially regarding the way MAF/DAF is estimated

#' Function to compute the number of missing data in a sample
#'
#' @param x # the sample to investigate
#'
#' @return number of missing data coded as "9" in the sample
#'
missing_data_computation <- function(x){
  length(which(x==9))
}


#' Function to compute MAF wih regards to missing data
#'
#' @param x # the matrix with data to use
#'
#' @return the MAF (estimated frequency) for the considered data
#'
maf.freq = function(x){
  x = x[x != 9]
  f = sum(x) ; (min(f , length(x) - f ))/length(x)
  }

#' Function to compute MAF wih regards to missing data.
#' CAUTION this statistics is intended to compute MAF within
#' the sample under investigation without regards to reference
#' allele.
#'
#' @param x # the matrix with data to use
#'
#' @return the MAF as allele counts for the considered data
#'
maf.count = function(x){
  x = x[x != 9]
  f = sum(x) ; min(f , length(x) - f ) }



#' Function to compute DAF wih regards to missing data
#'
#' @param x # the matrix with data to use
#'
#' @return the DAF (observed frequency) for the considered data
#'
daf.freq = function(x){
  x = x[x != 9]
  (length(x[x==1]))/length(x)
  }

#' Function to compute DAF wih regards to missing data.
#'
#' @param x # the matrix with data to use
#'
#' @return the DAF as allele counts for the considered data
#'
daf.count = function(x){
  x = x[x != 9]
  f = length(x[x==1])
  }
  
#' Function to compute SFS based on daf.count results
#'
#' @param x # the matrix with data to use
#'
#' @return the SFS with DAF as allele counts for the considered data
#'
sfs.daf.count = function(x){ apply(x,2,daf.count)}
