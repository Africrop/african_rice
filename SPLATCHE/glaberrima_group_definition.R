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

# Defining groups of genotypes for the computation of summary statistics to be used for spatial modelisation

## Loading required packages and datasets

library(ggmap)
library(cowplot)
library(ggrepel)
library(rgdal)
source("kmeans_group.r")

# We then launch the list of georeferenced genotypes.

list_og <- read.table("list-glab-georef.txt",header=TRUE,na.strings = "")


# We import shapefiles that will be used as the map background and convert them in a format usable with ggplot2.
lakes <- readOGR("naturalearth",layer = "ne_10m_lakes")
playas <- readOGR("naturalearth",layer = "ne_50m_playas")
ocean <- readOGR("naturalearth",layer = "ne_50m_ocean")

shapefile <- readOGR("naturalearth",layer = "ne_50m_admin_0_countries")
data <- fortify(shapefile)
data.lakes <- fortify(lakes)
data.playas <- fortify(playas)
data.ocean <- fortify(ocean)

## Geographical repartition of the samples, We draw a map of the genotypes.
map <- ggplot()+  
  geom_polygon(aes(x = long, y = lat, group = group),
               fill="lightblue", data = data.ocean, size = .3) +
  geom_polygon(aes(x = long, y = lat, group = group), data = data,
               colour ='antiquewhite4', fill = "antiquewhite", size = .3) +
  geom_polygon(aes(x = long, y = lat, group = group), data = data.lakes,
               colour ='lightblue',fill="lightblue", size = .3) +
  geom_polygon(aes(x = long, y = lat, group = group), data = data.playas,
               colour ='lightblue',fill="lightblue", size = .3) +
  geom_point(data = list_og, aes(x = Long_, y = lat_), size = 2) +
  geom_label_repel(data = list_og, aes(x = Long_, y = lat_),
                   label = list_og$code_vcf,
                   fontface = 'bold',
                   color = 'black',
                   size = 2,
                   box.padding = unit(0.005, "lines"),
                   point.padding = unit(0.5, "lines") ) +
  ggtitle('Geographic position of samples') +
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        legend.key = element_blank(),
        panel.border = element_rect(colour = "black",
                                    size = .5,
                                    linetype = "solid"
        )
  ) +
  coord_quickmap(xlim = c(-26,45),ylim =c(-20,25))
map


# To improve the kmeans clustering speed, we discard these 4 genotypes and apply a kmeans procedure with a rejection algorithm to the remaining genotypes. We used a customised function to perform this task, indicating the number of groups and the minimum number of genotypes we want.
list_og_without_outliers <- list_og[-which(list_og$code_vcf %in% c("IK","KC","EK","ND")),]

# We then perform the kmeans clustering. After some trial runs, we found that a good compromise between number of groups and minimum number of individuals per group was respectively 12 and 6. The code indicated below is for the sake of the report, it was ran only once to produce a final clustering that will be retained for analysis.
kmeans_grouping(coord = cbind(list_og_without_outliers$lat_,list_og_without_outliers$Long_),c = 12, m=6)
