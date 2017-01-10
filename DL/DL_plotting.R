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
# Written by Philippe Cubry, Benedicte Rhone
#
###################################################################################################################################


###################################################################################################################################
#
# This script performs the graphical representation of linkage disequilibrium
# from files containing mean r2 values for each position
#
###################################################################################################################################

glab.dl <- read.table("/Users/cubry/Documents/Riz_analysis/DL/glab.mean.R.2.txt",header=TRUE)
barth.dl <- read.table("/Users/cubry/Documents/Riz_analysis/DL/barthii.mean.R.2.txt",header=TRUE)

library(ggplot2)
library(cowplot)

##### Computation of Hill and Weir expectations - script adapted from Benedicte Rhone #####
n.barth<-83
n.glab<-163
 distance.barth<-barth.dl$Var1
 LD.data.barth<-barth.dl$tapply.barth.DL.R.2..barth.DL.dist..mean.

 distance.glab<-glab.dl$Var1
 LD.data.glab<-glab.dl$tapply.glab.DL.R.2..glab.DL.dist..mean.



HW.st<-c(C=0.2)                   
 HW.nonlinear.barth<-nls(LD.data.barth~(((10+C*distance.barth)/
                                          ((2+C*distance.barth)*(11+C*distance.barth)))*
                                          (1+((3+C*distance.barth)*(12+12*C*distance.barth+
                                          (C*distance.barth)^2))/(n.barth*(2+C*distance.barth)*
                                          (11+C*distance.barth)))), 
                                          start=HW.st, control=nls.control(maxiter=100),
                                          weights = barth.dl$Freq)

 HW.nonlinear.barth.prediction <- predict(HW.nonlinear.barth,list(distance = barth.dl$Var1))
 
summary(HW.nonlinear.barth)
tt<-summary(HW.nonlinear.barth)
new.rho<-tt$parameters[1]
new.rho

HW.nonlinear.glab<-nls(LD.data.glab~((10+C*distance.glab)/
                                       (((2+C*distance.glab)*(11+C*distance.glab)))*
                                       (1+((3+C*distance.glab)*(12+12*C*distance.glab+
                                       (C*distance.glab)^2))/(n.glab*(2+C*distance.glab)*
                                       (11+C*distance.glab)))),start=HW.st,
                                       control=nls.control(maxiter=100),
                                       weights = glab.dl$Freq)
HW.nonlinear.glab.prediction <- predict(HW.nonlinear.glab,list(distance = glab.dl$Var1))


##### Représenation graphique en fonction de la distance entre paires de SNP (en pb)
smoothness <- 30

p <- ggplot() +
  stat_smooth(aes(barth.dl$Var1,barth.dl$tapply.barth.DL.R.2..barth.DL.dist..mean.,
                  weight=barth.dl$Freq,color="blue",linetype = "dotted"),method = "gam",formula = y ~ s(x,bs="cs",k=smoothness),se = FALSE)+
  geom_line(aes(barth.dl$Var1,HW.nonlinear.barth.prediction,color="red",linetype = "dotted"))+
  stat_smooth(aes(x=glab.dl$Var1,y=glab.dl$tapply.glab.DL.R.2..glab.DL.dist..mean.,
                  weight=glab.dl$Freq,color="blue",linetype = "solid"),method = "gam",formula = y ~ s(x,bs="cs",k=smoothness),se = FALSE)+
  geom_line(aes(glab.dl$Var1,HW.nonlinear.glab.prediction),color="red",linetype = "solid")+
  scale_color_manual(name="Statistic",values=c("blue","red"),labels=c("Smoothed mean r2","Hill and Weir (1988) expectation"))+
  scale_linetype_manual(name="Species",values=c("dotted","solid"),labels=c("O. barthii","O. glaberrima"))+
  labs(x="Distance (in bp)",y= "r2")

ggsave(filename = "Genome wide LD Decay.png",p)
