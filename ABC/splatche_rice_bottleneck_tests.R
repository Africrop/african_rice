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

# load libraries
library(abc)
library(ggplot2)
library(raster)
library(rgdal)
library(cowplot)


# Load map elements - TO BE DOWNLOADED FROM NATURALEARTH website : http://www.naturalearthdata.com/
nat.earth <- stack('naturalearth/HYP_50M_SR_W.tif')
lakes <- readOGR("naturalearth",layer = "ne_10m_lakes")
rivers <- readOGR("naturalearth",layer = "ne_50m_rivers_lake_centerlines")
playas <- readOGR("/naturalearth",layer = "ne_50m_playas")
ocean <- readOGR("naturalearth",layer = "ne_10m_ocean")
shapefile <- readOGR("naturalearth",layer = "ne_50m_admin_0_countries")
data <- fortify(shapefile)
data.lakes <- fortify(lakes)
data.playas <- fortify(playas)
data.ocean <- fortify(ocean)
data.rivers <- fortify(rivers)
# Define extent and load raster for the map
extent <- c(-25,57,-5,40)
nat.crop <- crop(nat.earth, y=extent(extent))
rast.table <- data.frame(xyFromCell(nat.crop, 1:ncell(nat.crop)),
                         getValues(nat.crop/255))
rast.table$rgb <- with(rast.table, rgb(HYP_50M_SR_W.1,
                                       HYP_50M_SR_W.2,
                                       HYP_50M_SR_W.3,
                                      1))

rm(nat.earth,nat.crop)

# functions definitions
densities_plot <- function(prior = prior, post = post, title = "") {ggplot()  +
            #Plot the density graphs
            geom_density(aes(prior,fill="priors"),alpha=0.5) +
            geom_density(aes(post,fill="posteriors"),alpha=0.5) +
            #guide
            scale_fill_manual(values=c('#E69F00','#999999')) +
            #labels
            labs(title= title,x="") +
            #Define the theme  
            #theme_bw() +
            theme(panel.grid=element_blank(),legend.title=element_blank())
}

posteriorposition_map <- function(posteriorposition = posteriorposition,title="",n=50,h=25) {
  ggplot()+  
    geom_polygon(aes(x = long, y = lat, group = group),
                 fill="lightblue", data = data.ocean, size = .3) +
    geom_polygon(aes(x = long, y = lat, group = group), data = data,
                 colour ='antiquewhite4', fill = "white", size = .3) +
    geom_polygon(aes(x = long, y = lat, group = group), data = data.lakes,
                 colour ='lightblue',fill="lightblue", size = .3) +
    stat_density2d(data = posteriorposition,
                   aes(x=Long1,y=Lat1,  fill = ..level.., alpha = ..level..),
                   size = 0.1, geom = 'polygon',n = n,h=h) +
     geom_path(aes(x = long, y = lat, group = group), data = data.rivers,
                 colour ='lightblue') +
    scale_fill_gradient(low = "yellow", high = "red", guide = FALSE) +
    scale_alpha(range = c(0, 0.25), guide = FALSE) +
    ggtitle(title) +
    theme_bw() +
    theme(axis.line = element_blank(),
          axis.title = element_blank(),
          legend.key = element_blank(),
          panel.border = element_rect(colour = "black",
                                      size = .5,
                                      linetype = "solid"
          )
    ) + coord_quickmap(xlim = c(-20,50),ylim =c(-0,35))
}


posteriorposition_map_raster <- function(posteriorposition = posteriorposition,title="",n=50,h=25) {
  ggplot(data=rast.table,aes(x=x,y=y))+
    geom_tile(fill= rast.table$rgb)+
    geom_polygon(aes(x = long, y = lat, group = group), data = data,
                 colour ='antiquewhite4', fill = "NA", size = .3) +
    geom_polygon(aes(x = long, y = lat, group = group), data = data.lakes,
                 colour ='lightblue',fill="lightblue", size = .3) +
    stat_density2d(data = posteriorposition,
                   aes(x=Long1,y=Lat1,  fill = ..level.., alpha = ..level..),
                   size = 0.1, geom = 'polygon',n = n,h=h) +
    scale_fill_gradient(low = "green", high = "red",guide=FALSE) +
    scale_alpha(range = c(0.00, 0.25), guide = FALSE) +
    geom_point(aes(x=-4.540278,y=13.890556),size=1)+
    geom_path(aes(x = long, y = lat, group = group), data = data.rivers,
                 colour ='lightblue') +
    ggtitle(title) +
    theme_bw() +
    theme(axis.line = element_blank(),
          axis.title = element_blank(),
          legend.key = element_blank(),
          panel.border = element_rect(colour = "black",
                                      size = .5,
                                      linetype = "solid"
          )
    ) + coord_quickmap(xlim = c(-20,50),ylim =c(-0,35))
}

estimate_mode <- function(x) { # Function to estimate the mode of a distribution
  d <- density(x)
  d$x[which.max(d$y)]
}


# read simulation sumarry statistics results for model with 13 groups
param.moinsMR = rbind(
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_1.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_2.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_3.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_4.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_5.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_6.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_7.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_8.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_9.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_10.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_11.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_12.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_13.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_14.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_15.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_16.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_17.moinsMR_sum.stat.txt", header = T),
    read.table("~/Riz_analysis/splatche/riz_newpriorsB_18.moinsMR_sum.stat.txt", header = T),
    read.table("~/Riz_analysis/splatche/riz_newpriorsB_19.moinsMR_sum.stat.txt", header = T),
    read.table("~/Riz_analysis/splatche/riz_newpriorsB_20.moinsMR_sum.stat.txt", header = T),
    read.table("~/Riz_analysis/splatche/riz_newpriorsB_21.moinsMR_sum.stat.txt", header = T),
    read.table("~/Riz_analysis/splatche/riz_newpriorsB_22.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_23.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_24.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_25.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_26.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_29.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_30.moinsMR_sum.stat.txt", header = T),  
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_31.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_32.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_33.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_34.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_35.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_36.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_37.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_38.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_39.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_40.moinsMR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_41.moinsMR_sum.stat.txt", header = T)
  )

# Read corresponding observed summary statistics
stat.obs_moinsMR = read.table("~/Riz_analysis/splatche/cult_rice_stat.obs_moinsMR.txt",header = T)

# boxplot of simulated summary statistics with comparison with empirical summary statistics
boxplot(param.moinsMR[,(length(param.moinsMR)-19):(length(param.moinsMR))],
        las=2, main = "Simulated vs Observed sum stat",xaxt="n", ylim=c(0,1))
axis(1, at= 1:20,c(paste("SFS",seq(1:7)),paste("RV",seq(1:13))), las=2)
points(1:20,stat.obs_moinsMR,pch=3,col="red")
legend("topright","Observed",pch=3,col="red")

# ABC analysis
nn10pop <- abc(target = stat.obs_moinsMR[c(1:20)],sumstat = param.moinsMR[,c(13:32)],param = param.moinsMR[,c(1:9,11)],method = "neuralnet",tol = 0.001,numnet = 500)
save(nn10pop,file ="nn10pop")

# Plotting posteriors of the parameters
long <- densities_plot(param.moinsMR$Long1,nn10pop$adj.values[,"Long1"],title = "Longitude")
lat <- densities_plot(param.moinsMR$Lat1,nn10pop$adj.values[,"Lat1"],title = "Latitude")
gen <- densities_plot(param.moinsMR$generations,nn10pop$adj.values[,"generations"],title = "Exp_length")
acc <- densities_plot(param.moinsMR$accroissement,nn10pop$adj.values[,"accroissement"],title = "Growth rate")
mig <- densities_plot(param.moinsMR$migration,nn10pop$adj.values[,"migration"],title = "Migration rate")
sizeBeforeExp <- densities_plot(param.moinsMR$SizeBeforeExpansion,nn10pop$adj.values[,"SizeBeforeExpansion"],title = "Bott_size")
TimeBott <- densities_plot(param.moinsMR$TimeOfBottleneck,nn10pop$adj.values[,"TimeOfBottleneck"],title = "Bott_length")
TimeRecBott <- densities_plot(param.moinsMR$TimeForRecentBott,nn10pop$adj.values[,"TimeForRecentBott"],title = "Rec_bott_time")
AncSize <- densities_plot(param.moinsMR$AncestralSize,nn10pop$adj.values[,"AncestralSize"],title = "Anc_Size")
K <- densities_plot(param.moinsMR$MainCarryingCapacity,nn10pop$adj.values[,"MainCarryingCapacity"],title = "K")

prow<-plot_grid(lat + theme(legend.position="none"),
                long + theme(legend.position="none"),
                gen + theme(legend.position="none"),
                mig + theme(legend.position="none"),
                acc + theme(legend.position="none"),
                sizeBeforeExp + theme(legend.position="none"),
                TimeBott + theme(legend.position="none"),
                TimeRecBott + theme(legend.position="none"),
                AncSize + theme(legend.position="none"),
                K + theme(legend.position="none")
                )
grobs <- ggplotGrob(lat)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
plot_grid( prow, legend, rel_widths = c(3, .5))
save_plot(filename = "/Users/cubry/Dropbox/Tracing the origin and diversity of african rice using 263 wild and cultivated depthly resequenced accessions/modelisation/graphes geographic inference/Posteriors_13groups_0001.pdf",plot_grid( prow, legend, rel_widths = c(3, .5)),base_height = 8,base_aspect_ratio = 2)


```{r nn10pop_map4}
# Repr?sentation des densit?s sur la carte


map <- as.data.frame(nn10pop$adj.values[,1:2])
maps <- plot_grid(posteriorposition_map(map,n=100,h=30))
maps
#maps_raster <- plot_grid(posteriorposition_map_raster(map,n=100,h=30))
#ggsave("/Users/cubry/Desktop/map_b.png",maps_raster)

long <- qplot(nn10pop$adj.values[,1],geom="histogram",binwidth=4,xlab = "Posterior of longitude")
lat <- qplot(nn10pop$adj.values[,2],geom="histogram",binwidth=4,xlab = "Posterior of latitude")
save_plot(filename = "/Users/cubry/Dropbox/Tracing the origin and diversity of african rice using 263 wild and cultivated depthly resequenced accessions/modelisation/graphes geographic inference/Histos_coords_13groups_0001.pdf",plot_grid( long,lat),base_height = 4,base_aspect_ratio = 2)

```



```{r read_newpriorsB_withoutEast, echo=FALSE}
param.moinsEstMR = rbind(
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_1.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_2.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_3.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_4.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_5.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_6.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_7.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_7.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_8.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_9.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_10.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_11.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_12.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_13.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_14.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_15.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_16.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_17.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_18.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_19.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_20.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_21.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_22.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_23.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_24.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/riz_newpriorsB_25.moinsEst_MR_sum.stat.txt", header = T),
 read.table("~/Riz_analysis/splatche/riz_newpriorsB_26.moinsEst_MR_sum.stat.txt", header = T),
 read.table("~/Riz_analysis/splatche/riz_newpriorsB_29.moinsEst_MR_sum.stat.txt", header = T),
 read.table("~/Riz_analysis/splatche/riz_newpriorsB_30.moinsEst_MR_sum.stat.txt", header = T),
 read.table("~/Riz_analysis/splatche/riz_newpriorsB_31.moinsEst_MR_sum.stat.txt", header = T),
 read.table("~/Riz_analysis/splatche/riz_newpriorsB_32.moinsEst_MR_sum.stat.txt", header = T),
 read.table("~/Riz_analysis/splatche/riz_newpriorsB_33.moinsEst_MR_sum.stat.txt", header = T),
 read.table("~/Riz_analysis/splatche/riz_newpriorsB_34.moinsEst_MR_sum.stat.txt", header = T),
 read.table("~/Riz_analysis/splatche/riz_newpriorsB_35.moinsEst_MR_sum.stat.txt", header = T),
 read.table("~/Riz_analysis/splatche/riz_newpriorsB_36.moinsEst_MR_sum.stat.txt", header = T),
 read.table("~/Riz_analysis/splatche/riz_newpriorsB_37.moinsEst_MR_sum.stat.txt", header = T),
 read.table("~/Riz_analysis/splatche/riz_newpriorsB_38.moinsEst_MR_sum.stat.txt", header = T),
 read.table("~/Riz_analysis/splatche/riz_newpriorsB_39.moinsEst_MR_sum.stat.txt", header = T),
 read.table("~/Riz_analysis/splatche/riz_newpriorsB_40.moinsEst_MR_sum.stat.txt", header = T),
 read.table("~/Riz_analysis/splatche/riz_newpriorsB_41.moinsEst_MR_sum.stat.txt", header = T)
  
 )

stat.obs_moinsEst_MR = read.table("~/Riz_analysis/splatche/cult_rice_stat.obs_moinsEast_MR.txt",header = T)

boxplot(param.moinsEstMR[,(length(param.moinsEstMR)-23):(length(param.moinsEstMR)-5)],
        las=2, main = "Simulated vs Observed sum stat",xaxt="n", ylim=c(0,1))
axis(1, at= 1:19,c(paste("SFS",seq(1:7)),paste("RV",seq(1:12))), las=2)
points(1:19,stat.obs_moinsEst_MR,pch=3,col="red")
legend("topright","Obs",pch=3,col="red")
```
```{r nn11pop, include=FALSE}

#nn11pop <- abc(target = stat.obs_moinsEst_MR[c(1:19)],sumstat = param.moinsEstMR[,c(13:31)],param = param.moinsEstMR[,c(1:9,11)],method = "neuralnet",tol = 0.001,numnet = 500)
#save(nn11pop,file ="/Users/cubry/Documents/Riz_analysis/splatche/nn11pop")
load('/Users/cubry/Documents/Riz_analysis/splatche/nn11pop')
```
```{r nn11pop_posts, echo=FALSE}
long <- densities_plot(param.moinsEstMR$Long1,nn11pop$adj.values[,"Long1"],title = "Longitude")
lat <- densities_plot(param.moinsEstMR$Lat1,nn11pop$adj.values[,"Lat1"],title = "Latitude")
gen <- densities_plot(param.moinsEstMR$generations,nn11pop$adj.values[,"generations"],title = "Exp_length")
acc <- densities_plot(param.moinsEstMR$accroissement,nn11pop$adj.values[,"accroissement"],title = "Growth rate")
mig <- densities_plot(param.moinsEstMR$migration,nn11pop$adj.values[,"migration"],title = "Migration rate")
sizeBeforeExp <- densities_plot(param.moinsEstMR$SizeBeforeExpansion,nn11pop$adj.values[,"SizeBeforeExpansion"],title = "Bott_size")
TimeBott <- densities_plot(param.moinsEstMR$TimeOfBottleneck,nn11pop$adj.values[,"TimeOfBottleneck"],title = "Bott_length")
TimeRecBott <- densities_plot(param.moinsEstMR$TimeForRecentBott,nn11pop$adj.values[,"TimeForRecentBott"],title = "Rec_bott_time")
AncSize <- densities_plot(param.moinsEstMR$AncestralSize,nn11pop$adj.values[,"AncestralSize"],title = "Anc_size")
K <- densities_plot(param.moinsEstMR$MainCarryingCapacity,nn11pop$adj.values[,"MainCarryingCapacity"],title = "K")

prow<-plot_grid(lat + theme(legend.position="none"),
                long + theme(legend.position="none"),
                gen + theme(legend.position="none"),
                mig + theme(legend.position="none"),
                acc + theme(legend.position="none"),
                sizeBeforeExp + theme(legend.position="none"),
                TimeBott + theme(legend.position="none"),
                TimeRecBott + theme(legend.position="none"),
                AncSize + theme(legend.position="none"),
                K + theme(legend.position="none")
                )
grobs <- ggplotGrob(lat)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
plot_grid( prow, legend, rel_widths = c(3, .5))

save_plot(filename = "/Users/cubry/Dropbox/Tracing the origin and diversity of african rice using 263 wild and cultivated depthly resequenced accessions/modelisation/graphes geographic inference/Posteriors_12groups_0001.pdf",plot_grid( prow, legend, rel_widths = c(3, .5)),base_height = 8,base_aspect_ratio = 2)


```

```{r nn11pop_map5}
# Repr?sentation des densit?s sur la carte


map <- as.data.frame(nn11pop$adj.values[,1:2])



maps <- plot_grid(posteriorposition_map(map,n=100,h=30))
maps
#maps_raster <- plot_grid(posteriorposition_map_raster(map,n=100,h=30))
#ggsave("/Users/cubry/Desktop/map_withtEast.png",maps_raster)

long <- qplot(nn11pop$adj.values[,1],geom="histogram",binwidth=4,xlab = "Posterior of longitude")
lat <- qplot(nn11pop$adj.values[,2],geom="histogram",binwidth=4,xlab = "Posterior of latitude")
save_plot(filename = "/Users/cubry/Dropbox/Tracing the origin and diversity of african rice using 263 wild and cultivated depthly resequenced accessions/modelisation/graphes geographic inference/Histos_coords_12groups_0001.pdf",plot_grid( long,lat),base_height = 4,base_aspect_ratio = 2)

```



```{r validation 12 groups, include=FALSE}
# Choosing a simulation from each of the three possible regions for african rice domestication, only ran once !!!

# senegal <- which(param.moinsEstMR$Long1< -14 & param.moinsEstMR$Long1> -17 & param.moinsEstMR$Lat1 >11 & param.moinsEstMR$Lat1 <16)
#  senegal_position <- sample(senegal,size = 1)
#  while(length(which(is.na(param.moinsEstMR[senegal_position,])))!=0){
#    senegal_position <- sample(senegal,size = 1)
#  }
# 
# stat.obs_simul_west <- param.moinsEstMR[senegal_position,]
#write.table(stat.obs_simul_west,"/Users/cubry/Documents/Riz_analysis/splatche/simul_west_retained_12grps")
stat.obs_simul_west <-read.table("/Users/cubry/Documents/Riz_analysis/splatche/simul_west_retained_12grps")


# chad <- which(param.moinsEstMR$Long1>13 & param.moinsEstMR$Long1<15 & param.moinsEstMR$Lat1 >13 & param.moinsEstMR$Lat1 <15)
# chad_position <- sample(chad,size = 1)
# while(length(which(is.na(param.moinsEstMR[chad_position,])))!=0){
#   chad_position <- sample(chad,size = 1)
# }
#stat.obs_simul_east <- param.moinsEstMR[chad_position,]
#write.table(stat.obs_simul_east,"/Users/cubry/Documents/Riz_analysis/splatche/simul_east_retained_12grps")
stat.obs_simul_east <-read.table("/Users/cubry/Documents/Riz_analysis/splatche/simul_east_retained_12grps")

# mali <- which(param.moinsEstMR$Long1> -3 & param.moinsEstMR$Long1< -1 & param.moinsEstMR$Lat1 >16 & param.moinsEstMR$Lat1 <18)
# mali_position <- sample(mali,size = 1)
# while(length(which(is.na(param.moinsEstMR[mali_position,])))!=0){
#   mali_position <- sample(mali,size = 1)
# }
#stat.obs_simul_central <- param.moinsEstMR[mali_position,]
#write.table(stat.obs_simul_central,"/Users/cubry/Documents/Riz_analysis/splatche/simul_central_retained_12grps")
stat.obs_simul_central <- read.table("/Users/cubry/Documents/Riz_analysis/splatche/simul_central_retained_12grps")


nnpop.simul_west <- abc(target = stat.obs_simul_west[13:31],sumstat = param.moinsEstMR[,c(13:31)],param = param.moinsEstMR[,c(1:9,11)],method = "neuralnet",tol = 0.001,numnet = 50)
map <- as.data.frame(nnpop.simul_west$adj.values[,1:2])
maps <- plot_grid(posteriorposition_map(map,n=100,h=30)+geom_point(aes(x= stat.obs_simul_west$Long1,y=stat.obs_simul_west$Lat1)))
maps
ggsave("/Users/cubry/Dropbox/Tracing the origin and diversity of african rice using 263 wild and cultivated depthly resequenced accessions/modelisation/graphes geographic inference/map_simul_west.png",maps)

nnpop.simul_central <- abc(target = stat.obs_simul_central[13:31],sumstat = param.moinsEstMR[,c(13:31)],param = param.moinsEstMR[,c(1:9,11)],method = "neuralnet",tol = 0.001,numnet = 50)
map <- as.data.frame(nnpop.simul_central$adj.values[,1:2])
maps <- plot_grid(posteriorposition_map(map,n=100,h=30)+geom_point(aes(x= stat.obs_simul_central$Long1,y=stat.obs_simul_central$Lat1)))
maps
ggsave("/Users/cubry/Dropbox/Tracing the origin and diversity of african rice using 263 wild and cultivated depthly resequenced accessions/modelisation/graphes geographic inference/map_simul_central.png",maps)

nnpop.simul_east <- abc(target = stat.obs_simul_east[13:31],sumstat = param.moinsEstMR[,c(13:31)],param = param.moinsEstMR[,c(1:9,11)],method = "neuralnet",tol = 0.001,numnet = 50)
map <- as.data.frame(nnpop.simul_east$adj.values[,1:2])
maps <- plot_grid(posteriorposition_map(map,n=100,h=30)+geom_point(aes(x= stat.obs_simul_east$Long1,y=stat.obs_simul_east$Lat1)))
maps
ggsave("/Users/cubry/Dropbox/Tracing the origin and diversity of african rice using 263 wild and cultivated depthly resequenced accessions/modelisation/graphes geographic inference/map_simul_east.png",maps)
```

```{r read_ppc_12groups, echo=FALSE}
param.PPC_12groups = rbind(
  read.table("~/Riz_analysis/splatche/ppc_riz_12_post.moinsEst_MR_sum.stat.txt", header = T),
  read.table("~/Riz_analysis/splatche/ppc_riz_12_post_2.moinsEst_MR_sum.stat.txt", header = T)
  )

stat.obs_moinsEst_MR = read.table("~/Riz_analysis/splatche/cult_rice_stat.obs_moinsEast_MR.txt",header = T)

param.PPC_12groups.reshape <- reshape2::melt(param.PPC_12groups[,13:31])


ppc_12_groups <- qplot(data=param.PPC_12groups.reshape[-which(param.PPC_12groups.reshape[,"value"]==0|param.PPC_12groups.reshape[,"value"]>0.9),],
      x=variable,y=value,
      geom="violin") + geom_point(data = reshape2::melt(stat.obs_moinsEst_MR),aes(x=variable, y=value))+ geom_point(data = reshape2::melt(stat.obs_moinsEst_MR),aes(x=variable, y=value),col="red")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ylab("Distributions of summary statistics")

# boxplot(param.PPC_12groups[,13:31],
#         las=2, main = "Simulated vs Observed sum stat",xaxt="n", ylim=c(0,1))
# axis(1, at= 1:19,c(paste("SFS",seq(1:7)),paste("RV",seq(1:12))), las=2)
# points(1:19,stat.obs_moinsEst_MR,pch=3,col="red")
# legend("topright","Obs",pch=3,col="red")

ppc_12_ci <- apply(X = param.PPC_12groups[,13:31],MARGIN = 2,FUN=function(x){
  y <- x[x>0]
  y <- y[y<0.9]
  c(quantile(y,probs = c(0.001315789,1-0.001315789),na.rm=TRUE))
})
ppc_12_ci <- rbind(stat.obs_moinsEst_MR,ppc_12_ci)

row.names(ppc_12_ci) <- c("obs. value","lower 95% c.i. boundary, Bonferroni corrected","upper 95% c.i. boundary, Bonferroni corrected")
write.table(t(ppc_12_ci),"/Users/cubry/Dropbox/Tracing the origin and diversity of african rice using 263 wild and cultivated depthly resequenced accessions/modelisation/graphes geographic inference/PPC_12_CI.txt",quote=FALSE)
ggsave("/Users/cubry/Dropbox/Tracing the origin and diversity of african rice using 263 wild and cultivated depthly resequenced accessions/modelisation/graphes geographic inference/PPC_12_groups.png",ppc_12_groups)
```

```{r read_ppc_13groups, echo=FALSE}
param.PPC_13groups = rbind(
  read.table("~/Riz_analysis/splatche/ppc_riz_13_post_.moinsMR_sum.stat.txt", header = T)
  )

stat.obs_moinsMR = read.table("~/Riz_analysis/splatche/cult_rice_stat.obs_moinsMR.txt",header = T)

param.PPC_13groups.reshape <- reshape2::melt(param.PPC_13groups[,13:32])


ppc_13_groups <- qplot(data=param.PPC_13groups.reshape[-which(param.PPC_13groups.reshape[,"value"]==0|param.PPC_13groups.reshape[,"value"]>0.9),],
      x=variable,y=value,
      geom="violin") + geom_point(data = reshape2::melt(stat.obs_moinsMR),aes(x=variable, y=value))+ geom_point(data = reshape2::melt(stat.obs_moinsMR),aes(x=variable, y=value),col="red")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ylab("Distributions of summary statistics")

# boxplot(param.PPC_13groups[,13:32],
#         las=2, main = "Simulated vs Observed sum stat",xaxt="n", ylim=c(0,1))
# axis(1, at= 1:20,c(paste("SFS",seq(1:7)),paste("RV",seq(1:13))), las=2)
# points(1:20,stat.obs_moinsMR,pch=3,col="red")
# legend("topright","Obs",pch=3,col="red")

ppc_13_ci <- apply(X = param.PPC_13groups[,13:32],MARGIN = 2,FUN=function(x){
  y <- x[x>0]
  y <- y[y<0.9]
  c(quantile(y,probs = c(0.00125,1-0.00125),na.rm=TRUE))
})
ppc_13_ci <- rbind(stat.obs_moinsMR,ppc_13_ci)

row.names(ppc_13_ci) <- c("obs. value","lower 95% c.i. boundary, Bonferroni corrected","upper 95% c.i. boundary, Bonferroni corrected")
write.table(t(ppc_13_ci),"/Users/cubry/Dropbox/Tracing the origin and diversity of african rice using 263 wild and cultivated depthly resequenced accessions/modelisation/graphes geographic inference/PPC_13_CI.txt",quote=FALSE)
ggsave("/Users/cubry/Dropbox/Tracing the origin and diversity of african rice using 263 wild and cultivated depthly resequenced accessions/modelisation/graphes geographic inference/PPC_13_groups.png",ppc_13_groups)
```