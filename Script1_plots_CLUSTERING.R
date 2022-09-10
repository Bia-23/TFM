### LIBRARIES ######
library(shapes)
library(factoextra)
library(pca3d)
library(rgl)
library(ggplot2)
library(corrplot)
library(geomorph)
library(gridExtra)
library(rpart)
library(rpart.plot)
library(pracma)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(MASS)
library(candisc)
library(psych)
library(caret)
library(rstatix)
library(FactoMineR)
library(cluster)
####################

# Import data
data = read.csv("0_957_2022_left_arm.csv", sep=";", header=T, dec = ",")

# Correction of exchanged coordinates:
aux=data[933:957,c(22:25,30:33)]
data[933:957,22]=aux[5] # 22 -> 30
data[933:957,23]=aux[6] # 23 -> 31
data[933:957,24]=aux[7] # 24 -> 32
data[933:957,25]=aux[8] # 25 -> 33
data[933:957,30]=aux[1] # 30 -> 22
data[933:957,31]=aux[2] # 31 -> 23
data[933:957,32]=aux[3] # 32 -> 24
data[933:957,33]=aux[4] # 33 -> 25

######################
### DATA CLEANING #### 
######################

# - Age > 60 months
# - z-scores >5, <-5
# - outliers
data = data[(data$agemons<=60),]
data = data[(data$zlen<=5 & data$zlen>=-5),]
data = data[(data$zwei<=5 & data$zwei>=-5),]
data = data[(data$zwfl<=5 & data$zwfl>=-5),]
data = data[(data$zbmi<=5 & data$zbmi>=-5),]
data = data[(data$zac<=5 & data$zac>=-5),]

boxplot(muac)$out
boxplot(weight)$out
boxplot(height)$out

# Categorical variables
data$sex = as.factor(data$sex)
data$measure = as.factor(data$measure)
data$class_MUAC = as.factor(data$class_MUAC)
data$class_WFL = as.factor(data$class_WFL)
data$class_global = as.factor(data$class_global)

### Generalized Procrustes Analysis ####
n=nrow(data)
# 20 landmarks, 2 dim
X_raw = array(0,dim=c(20,2,n))
X_Proc = array(0,dim=c(20,2,n))
for (i in 1:n){
  X_raw[,,i] = cbind(t(data[i,seq(2,40,2)]),-t(data[i,seq(3,41,2)]))
  X_Proc[,,i] = cbind(t(data[i,seq(44,82,2)]),-t(data[i,seq(45,83,2)]))
}

X_Proc_R = procGPA(X_raw)

par(mfrow=c(1,1))
plotAllSpecimens(X_Proc_R$rotated[,,1:2],mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                                  4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),label=TRUE,
                 plot.param = list(pt.cex=.0,mean.cex=.6,link.lwd=1,txt.cex=2,txt.col="red",txt.pos="1"))

shapepca(X_Proc_R, pcno = c(1,2,3), type = "r", mag = 1, joinline=c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),
         project=c(1,2),scores3d=FALSE,color=2,axes3=FALSE,rglopen=TRUE,zslice=0)

Proc_coords_R=as.data.frame(cbind(t(X_Proc_R$rotated[1,,]),t(X_Proc_R$rotated[2,,]),t(X_Proc_R$rotated[3,,]),
                                  t(X_Proc_R$rotated[4,,]),t(X_Proc_R$rotated[5,,]),t(X_Proc_R$rotated[6,,]),
                                  t(X_Proc_R$rotated[7,,]),t(X_Proc_R$rotated[8,,]),t(X_Proc_R$rotated[9,,]),
                                  t(X_Proc_R$rotated[10,,]),t(X_Proc_R$rotated[11,,]),t(X_Proc_R$rotated[12,,]),
                                  t(X_Proc_R$rotated[13,,]),t(X_Proc_R$rotated[14,,]),t(X_Proc_R$rotated[15,,]),
                                  t(X_Proc_R$rotated[16,,]),t(X_Proc_R$rotated[17,,]),t(X_Proc_R$rotated[18,,]),
                                  t(X_Proc_R$rotated[19,,]),t(X_Proc_R$rotated[20,,])))
data[,44:83]=Proc_coords_R

attach(data)


# Create new categorical variables
data$age=ifelse(agemons<=12,1,
                ifelse(agemons>12&agemons<=24,2,
                       ifelse(agemons>24&agemons<=36,3,
                              ifelse(agemons>36&agemons<=48,4,5))))
data$age=as.factor(data$age)
data$zlen_cat=ifelse(zlen<=(-3),'SAM',
                     ifelse(zlen>(-3)&zlen<=(-2),'MAM',
                            ifelse(zlen>(-2)&zlen<=(-1),'RIS','NOR')))
data$zlen_cat=as.factor(data$zlen_cat)
data$zwei_cat=ifelse(zwei<=(-3),'SAM',
                     ifelse(zwei>(-3)&zwei<=(-2),'MAM',
                            ifelse(zwei>(-2)&zwei<=(-1),'RIS','NOR')))
data$zwei_cat=as.factor(data$zwei_cat)
data$zwfl_cat=ifelse(zwfl<=(-3),'SAM',
                     ifelse(zwfl>(-3)&zwfl<=(-2),'MAM',
                            ifelse(zwfl>(-2)&zwfl<=(-1),'RIS','NOR')))
data$zwfl_cat=as.factor(data$zwfl_cat)
data$zbmi_cat=ifelse(zbmi<=(-3),'SAM',
                     ifelse(zbmi>(-3)&zbmi<=(-2),'MAM',
                            ifelse(zbmi>(-2)&zbmi<=(-1),'RIS','NOR')))
data$zbmi_cat=as.factor(data$zbmi_cat)
data$zac_cat=ifelse(zac<=(-3),'SAM',
                    ifelse(zac>(-3)&zac<=(-2),'MAM',
                           ifelse(zac>(-2)&zac<=(-1),'RIS','NOR')))
data$zac_cat=as.factor(data$zac_cat)

# 3 nutritional statuses
data$class_MUAC_3=ifelse(class_MUAC=="SAM",'SAM+MAM',
                         ifelse(class_MUAC=="MAM",'SAM+MAM',
                                ifelse(class_MUAC=="RIS",'RIS','NOR')))
data$class_MUAC_3=as.factor(data$class_MUAC_3)
data$class_WFL_3=ifelse(class_WFL=="SAM",'SAM+MAM',
                        ifelse(class_WFL=="MAM",'SAM+MAM',
                               ifelse(class_WFL=="RIS",'RIS','NOR')))
data$class_WFL_3=as.factor(data$class_WFL_3)
data$class_global_3=ifelse(class_global=="SAM",'SAM+MAM',
                           ifelse(class_global=="MAM",'SAM+MAM',
                                  ifelse(class_global=="RIS",'RIS','NOR')))
data$class_global_3=as.factor(data$class_global_3)

# 2 nutritional statuses
data$class_MUAC_2=ifelse(class_MUAC=="SAM",'SAM+MAM',
                         ifelse(class_MUAC=="MAM",'SAM+MAM',
                                ifelse(class_MUAC=="RIS",'RIS+NOR','RIS+NOR')))
data$class_MUAC_2=as.factor(data$class_MUAC_2)
data$class_WFL_2=ifelse(class_WFL=="SAM",'SAM+MAM',
                        ifelse(class_WFL=="MAM",'SAM+MAM',
                               ifelse(class_WFL=="RIS",'RIS+NOR','RIS+NOR')))
data$class_WFL_2=as.factor(data$class_WFL_2)
data$class_global_2=ifelse(class_global=="SAM",'SAM+MAM',
                           ifelse(class_global=="MAM",'SAM+MAM',
                                  ifelse(class_global=="RIS",'RIS+NOR','RIS+NOR')))
data$class_global_2=as.factor(data$class_global_2)

# Reordering of levels
data$class_MUAC = factor(data$class_MUAC, levels=c("NOR","RIS","MAM","SAM"))
data$class_WFL = factor(data$class_WFL, levels=c("NOR","RIS","MAM","SAM"))
data$class_global = factor(data$class_global, levels=c("NOR","RIS","MAM","SAM"))
data$zlen_cat = factor(data$zlen_cat, levels=c("NOR","RIS","MAM","SAM"))
data$zwei_cat = factor(data$zwei_cat, levels=c("NOR","RIS","MAM","SAM"))
data$zwfl_cat = factor(data$zwfl_cat, levels=c("NOR","RIS","MAM","SAM"))
data$zbmi_cat = factor(data$zbmi_cat, levels=c("NOR","RIS","MAM","SAM"))
data$zac_cat = factor(data$zac_cat, levels=c("NOR","RIS","MAM","SAM"))

attach(data)

##################################
### GRAPHICAL REPRESENTATION #####
##################################

v = seq(1,40,by=2)
vv = seq(2,40,by=2)

# All coordinates
coord_x = as.matrix(Proc_coords_R[,v]) # x
coord_y = as.matrix(Proc_coords_R[,vv]) # y
plot(coord_x,coord_y,main="All coordinates for the 957 observations")


# Coordinates by groups
with(Proc_coords_R, plot(coord_x,coord_y, col=data$sex))
with(Proc_coords_R, plot(coord_x,coord_y, col=data$class_WFL)) # possible groups
with(Proc_coords_R, plot(coord_x,coord_y, col=data$class_WFL_3)) 
with(Proc_coords_R, plot(coord_x,coord_y, col=data$class_WFL_2)) 

with(Proc_coords_R, plot(coord_x,coord_y, col=data$class_MUAC))
with(Proc_coords_R, plot(coord_x,coord_y, col=data$class_MUAC_3))
with(Proc_coords_R, plot(coord_x,coord_y, col=data$class_MUAC_2)) 

with(Proc_coords_R, plot(coord_x,coord_y, col=data$class_global)) # possible groups
with(Proc_coords_R, plot(coord_x,coord_y, col=data$class_global_3)) 
with(Proc_coords_R, plot(coord_x,coord_y, col=data$class_global_2)) 

with(Proc_coords_R, plot(coord_x,coord_y, col=data$age)) 
with(Proc_coords_R, plot(coord_x,coord_y, col=data$measure)) # clear differences


# Histograms of Procrustes coordinates
par(mfrow=c(4,5))
sapply(names(data[,c(43:84)]),function(cname){hist(data[,c(43:84)][[cname]],main=cname,col="lightblue")})



# Histograms of numerical variables
par(mfrow=c(1,3))
sapply(names(data[,c(84,86:87)]),function(cname){hist(data[,c(84,86:87)][[cname]],main=cname,col="lightblue")})

# Barplots
barplot(table(data$class_MUAC_2),col=c("palegreen","indianred1"),main="class_MUAC_2")
barplot(table(data$class_WFL_2),col=c("palegreen","indianred1"),main="class_WFL_2")
barplot(table(data$class_global_2),col=c("palegreen","indianred1"),main="class_global_2")

barplot(table(data$class_MUAC_3),col=c("palegreen","indianred1","lightsteelblue"),main="class_MUAC_3")
barplot(table(data$class_WFL_3),col=c("palegreen","indianred1","lightsteelblue"),main="class_WFL_3")
barplot(table(data$class_global_3),col=c("palegreen","indianred1","lightsteelblue"),main="class_global_3")

# datasets male vs female
male=data[sex==1,]
fem=data[sex==2,]

mean(fem$weight)
mean(male$weight)

# histograms: weight by sex
df = data[,c('sex','weight')]
ggplot(df, aes(x=weight, color=sex)) +
  geom_histogram(fill="white", alpha=0.5, position="identity",bins=25)

# histograms: height by sex
df = data[,c('sex','height')]
ggplot(df, aes(x=height, color=sex)) +
  geom_histogram(fill="white", alpha=0.5, position="identity",bins=25)

range(fem$height)
range(male$height)

mean(fem$muac)
mean(male$muac)

# weight vs age
par(mfrow=c(1,1))
plot(agemons,weight)
abline(lm(weight~agemons),col="red")

# histogram muac by sex
df = data[,c('sex','muac')]
ggplot(df, aes(x=muac, color=sex)) +
  geom_histogram(fill="white", alpha=0.5, position="identity",bins=25)

par(mfrow=c(1,1))
plot(muac,weight)
abline(lm(weight~muac), col="red")
with(data,plot(muac,weight,col=age))

# OUTLIERS
plot(muac,zwfl,xlab = "MUAC", ylab="Weight for length (ZWFL)",cex.lab=1.3)
abline(lm(zwfl~muac), col="red")
abline(v=12.5)
abline(h=-2)

# age groups
with(data,plot(muac,zwfl,col=age))
with(data,plot(muac,zwfl,col=class_global))
with(data,plot(muac,zwfl,col=class_MUAC))

# muac vs height
plot(muac,height)
abline(lm(height~muac), col="red")


#####################
### CORRELATIONS ####
#####################

corrplot(cor(data[,c(44:84,86:88)]),order='hclust',type="lower",tl.pos = 'l',tl.cex=0.66)

###########################
### K-MEANS CLUSTERING ####
###########################

scaled_coords=scale(Proc_coords_R,center = TRUE,scale = TRUE)

k2 = kmeans(scaled_coords, centers = 2, nstart = 25)
k3 = kmeans(scaled_coords, centers = 3, nstart = 25)
k4 = kmeans(scaled_coords, centers = 4, nstart = 25)
k5 = kmeans(scaled_coords, centers = 5, nstart = 25)
# plots
p1 = fviz_cluster(k2, geom = "point", data = scaled_coords) + ggtitle("k = 2")
p2 = fviz_cluster(k3, geom = "point",  data = scaled_coords) + ggtitle("k = 3")
p3 = fviz_cluster(k4, geom = "point",  data = scaled_coords) + ggtitle("k = 4")
p4 = fviz_cluster(k5, geom = "point",  data = scaled_coords) + ggtitle("k = 5")

grid.arrange(p1, p2, p3, p4, nrow = 2) 

# k optimal
fviz_nbclust(scaled_coords, kmeans, method = "silhouette") # 2 clusters

data[,'cluster']=k2$cluster # numero de cluster correspondiente a cada observación
data$cluster=as.factor(data$cluster)

# datos del custer 1 vs cluster 2 por separado:
clust1=data[k2$cluster==1,]
clust2=data[k2$cluster==2,]

# Mean shape of each cluster:
n1=nrow(clust1) # 594
n2=nrow(clust2) # 341
X_Proc_c1 = array(0,dim=c(20,2,n1))
X_Proc_c2 = array(0,dim=c(20,2,n2))
for (i in 1:n1){
  X_Proc_c1[,,i] = cbind(t(clust1[i,seq(44,82,2)]),t(clust1[i,seq(45,83,2)]))
}

for (i in 1:n2){
  X_Proc_c2[,,i] = cbind(t(clust2[i,seq(44,82,2)]),t(clust2[i,seq(45,83,2)]))
}


plotAllSpecimens(X_Proc_c1,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                    4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0.5,mean.cex=0.5,link.lwd=0.5))
plotAllSpecimens(X_Proc_c2,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                    4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0.5,mean.cex=0.5,link.lwd=0.5))


# media de todas las observaciones y encima media del cluster 1:
C1=mshape(X_Proc_R$rotated[,,k2$cluster==1])
C2=mshape(X_Proc_R$rotated[,,k2$cluster==2])
joinline=c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1)
par(mfrow=c(1,2))
#,link.lty="dotted"
plotAllSpecimens(X_Proc_R$rotated,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                           4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0,mean.cex=0.5,link.lwd=0.5))
lines(C1[joinline,],col="lightsteelblue3",lwd=2) # mean shape cluster 1 
title(main="Mean shape cluster 1",cex.main=1.7)


plotAllSpecimens(X_Proc_R$rotated,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                           4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0,mean.cex=0.5,link.lwd=0.5))
lines(C2[joinline,],col="lightsteelblue3",lwd=2) # mean shape cluster 2 
title(main="Mean shape cluster 2",cex.main=1.7)


# GRAPHICAL ANALYSIS
# Factor variables
par(mfrow=c(1,1))
par(xpd=TRUE) 
levels(data$sex)=c('Male','Female')
barplot(with(data,(prop.table(table(data$sex,data$clust),2))),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue3"),cex.lab=1.5,cex.main=2.1,main="Sex")
legend("topright", legend = levels(data$sex) , 
       col=c("aliceblue","lightblue3"), 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.8, horiz = FALSE, inset = c(0.0, 0.0))
# test 
with(data,chisq.test(sex,cluster,correct=TRUE,p=rep(1/length(sex),length(sex))))

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =levels(data$class_MUAC), pch=16, pt.cex=3, cex=1.5, bty='n',
       col = c("aliceblue","lightblue2","lightblue3","lightblue4"),horiz=TRUE)

par(mfrow=c(1,3))
barplot(with(data,(prop.table(table(data$class_MUAC,data$clust),2))),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),cex.lab=1.5,cex.main=2.1,main="Class MUAC")
legend("topright", legend = levels(data$class_MUAC) , 
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"), 
       bty = "n", pch=20 , pt.cex = 2, cex =1.5, horiz = FALSE, inset = c(0.0, 0.0))

barplot(with(data,(prop.table(table(data$class_WFL,data$clust),2))),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),cex.lab=1.5,cex.main=2.1,cex.main=2.1,main="Class WFL")
legend("topright", legend = levels(data$class_WFL) , 
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"), 
       bty = "n", pch=20 , pt.cex = 2, cex = 1.5, horiz = FALSE, inset = c(0.0, 0.0))

barplot(with(data,(prop.table(table(data$class_global,data$clust),2))),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),cex.lab=1.5,cex.main=2.1,main="Class global")
legend("topright", legend = levels(data$class_global) , 
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"), 
       bty = "n", pch=20 , pt.cex = 2, cex = 1.5, horiz = FALSE, inset = c(0.0, 0.0))

par(mfrow=c(2,3))
barplot(with(data,(prop.table(table(data$zlen_cat,data$clust),2))),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),cex.lab=1.5,cex.main=2.1,main="Length for age (zlen)")
barplot(with(data,(prop.table(table(data$zwei_cat,data$clust),2))),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),cex.lab=1.5,cex.main=2.1,main="Weight for age (zwei)")
barplot(with(data,(prop.table(table(data$zwfl_cat,data$clust),2))),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),cex.lab=1.5,cex.main=2.1,main="Weight for length (zwfl)")
barplot(with(data,(prop.table(table(data$zbmi_cat,data$clust),2))),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),cex.lab=1.5,cex.main=2.1,main="BMI for age (zbmi)")
barplot(with(data,(prop.table(table(data$zac_cat,data$clust),2))),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),cex.lab=1.5,cex.main=2.1,main="MUAC for age (zac)")
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend =levels(data$class_MUAC), pch=16, pt.cex=4, cex=1.8, bty='n',
       col = c("aliceblue","lightblue2","lightblue3","lightblue4"),horiz=FALSE)


df = data[,c('cluster','agemons')]
p1=ggplot(df, aes(x=agemons, color=cluster)) + ggtitle("Age") +
  geom_histogram(fill="white", alpha=0.6, position="identity",bins=15)

df = data[,c('cluster','height')]
p2=ggplot(df, aes(x=height, color=cluster)) + ggtitle("Height") +
  geom_histogram(fill="white", alpha=0.6, position="identity",bins=15)

df = data[,c('cluster','weight')]
p3=ggplot(df, aes(x=weight, color=cluster)) + ggtitle("Weight") +
  geom_histogram(fill="white", alpha=0.6, position="identity",bins=15)

grid.arrange(p1, p2, p3, ncol = 4, nrow = 2, layout_matrix= rbind(c(NA,1,1,NA),c(2,2,3,3)))

# z-scores
ggplot(data, aes(x = zlen)) +
  geom_histogram(fill = "lightblue3", colour = "lightblue4", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(data, aes(x = zwei)) +
  geom_histogram(fill = "lightblue3", colour = "lightblue4", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(data, aes(x = zwfl)) +
  geom_histogram(fill = "lightblue3", colour = "lightblue4", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(data, aes(x = zbmi)) +
  geom_histogram(fill = "lightblue3", colour = "lightblue4", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(data, aes(x = zac)) +
  geom_histogram(fill = "lightblue3", colour = "lightblue4", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)

par(mfrow=c(2,4))
for (i in c(c(44:84),c(86,c(90:94)))){
  hist(clust1[,i],col = rgb(0,1,1,alpha = 0.6),cex.lab=1.5,cex.main=2.1,main="Cluster1",xlab=paste("Coord",i-43),ylim=c(0,270),breaks=15)
  hist(clust2[,i],add=TRUE,col = rgb(1,1,0,alpha = 0.2),cex.lab=1.5,cex.main=2.1,main="Cluster2",xlab=paste("Coord",i-43),ylim=c(0,270),breaks=15)
}

# Decision tree
rpart.fit = rpart(cluster ~., method="class", data = data[,-c(c(1:83),87,99:109)]) # c(c(1:83),
# summary(rpart.fit)
par(mfrow=c(1,1))
rpart.plot(rpart.fit, digits = 2, fallen.leaves = TRUE,
           type = 5, extra=101, box.palette = "Blues")

df = data.frame(imp = rpart.fit$variable.importance)
df2 = df %>% 
  tibble::rownames_to_column() %>% 
  dplyr::rename("variable" = rowname) %>% 
  dplyr::arrange(imp) %>%
  dplyr::mutate(variable = forcats::fct_inorder(variable))

ggplot(df2) +
  geom_segment(aes(x = variable, y = 0, xend = variable, yend = imp), 
               size = 1, alpha = 0.7) +
  geom_point(aes(x = variable, y = imp, col = variable), 
             size = 2, show.legend = F) +
  coord_flip() +
  theme_bw()


### AGE UNDER AND ABOVE 17 MONTHS ####

Age1=data[agemons<=17,] 
Age2=data[agemons>17,]

# AGE<17 ##################
scaled_coords1=scale(Age1[,c(44:83)],center = TRUE,scale = TRUE)

k2_1 = kmeans(scaled_coords1, centers = 2, nstart = 25)
k3_1 = kmeans(scaled_coords1, centers = 3, nstart = 25)
k4_1 = kmeans(scaled_coords1, centers = 4, nstart = 25)
k5_1 = kmeans(scaled_coords1, centers = 5, nstart = 25)
# plots
p1_1 = fviz_cluster(k2_1, geom = "point", data = scaled_coords1) + ggtitle("k = 2")
p2_1 = fviz_cluster(k3_1, geom = "point",  data = scaled_coords1) + ggtitle("k = 3")
p3_1 = fviz_cluster(k4_1, geom = "point",  data = scaled_coords1) + ggtitle("k = 4")
p4_1 = fviz_cluster(k5_1, geom = "point",  data = scaled_coords1) + ggtitle("k = 5")
grid.arrange(p1_1, p2_1, p3_1, p4_1, nrow = 2) 

# k óptimo
fviz_nbclust(scaled_coords1, kmeans, method = "silhouette") # 2 clusters

Age1[,'cluster']=k2_1$cluster
Age1$cluster=as.factor(Age1$cluster)

# Mean shape of each cluster:
n1=nrow(clust1_1)
n2=nrow(clust2_1)
X_Proc_c1 = array(0,dim=c(20,2,n1))
X_Proc_c2 = array(0,dim=c(20,2,n2))
for (i in 1:n1){
  X_Proc_c1[,,i] = cbind(t(clust1_1[i,seq(44,82,2)]),t(clust1_1[i,seq(45,83,2)]))
}

for (i in 1:n2){
  X_Proc_c2[,,i] = cbind(t(clust2_1[i,seq(44,82,2)]),t(clust2_1[i,seq(45,83,2)]))
}



plotAllSpecimens(X_Proc_c1,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                    4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0.5,mean.cex=0.5,link.lwd=0.5))

plotAllSpecimens(X_Proc_c2,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                    4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0.5,mean.cex=0.5,link.lwd=0.5))


# Mean arm shape
C1=mshape(X_Proc_c1)
C2=mshape(X_Proc_c2)

ProcCoords_Age1=array(c(X_Proc_c1, X_Proc_c2), dim=c(20,2,n1+n2))

par(mfrow=c(1,2))
plotAllSpecimens(ProcCoords_Age1,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                          4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0,mean.cex=0.5,link.lwd=0.5))
lines(C1[joinline,],col="lightsteelblue3",lwd=2) # mean shape cluster 1 
title(main="Mean shape cluster 1",cex.main=1.7)


plotAllSpecimens(ProcCoords_Age1,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                          4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0,mean.cex=0.5,link.lwd=0.5))
lines(C2[joinline,],col="lightsteelblue3",lwd=2) # mean shape cluster 2 
title(main="Mean shape cluster 2",cex.main=1.7)



par(mfrow=c(1,1))
levels(Age1$sex)=c('Male','Female')
barplot(with(Age1,prop.table(table(Age1$sex,Age1$cluster),2)),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue3"),cex.lab=1.5,cex.main=2.1,main="sex")
legend("topright",legend = levels(Age1$sex),bty="n",pch=20,pt.cex=2,cex=1.2,horiz=FALSE,inset=c(-0.3,0.0),
       col=c("aliceblue","lightblue3"))

(table(clust1_1$sex))
(table(clust2_1$sex))

# test
with(Age1,chisq.test(sex,cluster,correct=TRUE,p=rep(1/length(sex),length(sex))))

par(mfrow=c(1,3))
barplot(with(Age1,prop.table(table(Age1$class_MUAC,Age1$cluster),2)),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),cex.lab=1.5,cex.main=2.1,main="Class MUAC")
legend("topright",legend = levels(Age1$class_MUAC),bty="n",pch=20,pt.cex=2,cex=1.2,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))
barplot(with(Age1,prop.table(table(Age1$class_WFL,Age1$cluster),2)),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),cex.lab=1.5,cex.main=2.1,main="Class WFL")
legend("topright",legend = levels(Age1$class_WFL),bty="n",pch=20,pt.cex=2,cex=1.2,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))
barplot(with(Age1,prop.table(table(Age1$class_global,Age1$cluster),2)),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),cex.lab=1.5,cex.main=2.1,main="Class global")
legend("topright",legend = levels(Age1$class_global),bty="n",pch=20,pt.cex=2,cex=1.2,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

par(mfrow=c(2,3))
barplot(with(Age1,prop.table(table(Age1$zlen_cat,Age1$cluster),2)),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),cex.lab=1.5,cex.main=2.1,main="Length for age (zlen)")
barplot(with(Age1,prop.table(table(Age1$zwei_cat,Age1$cluster),2)),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),cex.lab=1.5,cex.main=2.1,main="Weight for age (zwei)")
barplot(with(Age1,prop.table(table(Age1$zwfl_cat,Age1$cluster),2)),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),cex.lab=1.5,cex.main=2.1,main="Weight for length (zwfl)")
barplot(with(Age1,prop.table(table(Age1$zac_cat,Age1$cluster),2)),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),cex.lab=1.5,cex.main=2.1,main="MUAC for age (zac)")
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend =levels(data$class_MUAC), pch=16, pt.cex=4, cex=1.8, bty='n',
       col = c("aliceblue","lightblue2","lightblue3","lightblue4"),horiz=FALSE)

df = Age1[,c('cluster','agemons')]
p1=ggplot(df, aes(x=agemons, color=cluster)) + ggtitle("Age") +
  geom_histogram(fill="white", alpha=0.6, position="identity",bins=15)

df = Age1[,c('cluster','height')]
p2=ggplot(df, aes(x=height, color=cluster)) + ggtitle("Height") +
  geom_histogram(fill="white", alpha=0.6, position="identity",bins=15)

df = Age1[,c('cluster','weight')]
p3=ggplot(df, aes(x=weight, color=cluster)) + ggtitle("Weight") +
  geom_histogram(fill="white", alpha=0.6, position="identity",bins=15)

grid.arrange(p1, p2, p3, ncol = 4, nrow = 2, layout_matrix= rbind(c(NA,1,1,NA),c(2,2,3,3)))

ggplot(Age1, aes(x = muac)) +
  geom_histogram(fill = "lightblue3", colour = "lightblue4", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
# z-scores
ggplot(Age1, aes(x = zlen)) +
  geom_histogram(fill = "lightblue3", colour = "lightblue4", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age1, aes(x = zwei)) +
  geom_histogram(fill = "lightblue3", colour = "lightblue4", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age1, aes(x = zwfl)) +
  geom_histogram(fill = "lightblue3", colour = "lightblue4", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age1, aes(x = zbmi)) +
  geom_histogram(fill = "lightblue3", colour = "lightblue4", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age1, aes(x = zac)) +
  geom_histogram(fill = "lightblue3", colour = "lightblue4", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)


# Decision tree
par(mfrow=c(1,1))
rpart.fit = rpart(cluster ~., method="class", data = Age1[,-c(c(1:83),87,104:109)]) # [,-c(c(1:83))]
# summary(rpart.fit)
rpart.plot(rpart.fit, digits = 2, fallen.leaves = TRUE,
           type = 5, extra=101, box.palette = "Blues")

df = data.frame(imp = rpart.fit$variable.importance)
df2 = df %>% 
  tibble::rownames_to_column() %>% 
  dplyr::rename("variable" = rowname) %>% 
  dplyr::arrange(imp) %>%
  dplyr::mutate(variable = forcats::fct_inorder(variable))
ggplot(df2) +
  geom_segment(aes(x = variable, y = 0, xend = variable, yend = imp), 
               size = 1, alpha = 0.7) +
  geom_point(aes(x = variable, y = imp, col = variable), 
             size = 2, show.legend = F) +
  coord_flip() +
  theme_bw()


# AGE>=17 ##################
scaled_coords2=scale(Age2[,c(44:83)],center = TRUE,scale = TRUE)

k2_2 = kmeans(scaled_coords2, centers = 2, nstart = 25)
k3_2 = kmeans(scaled_coords2, centers = 3, nstart = 25)
k4_2 = kmeans(scaled_coords2, centers = 4, nstart = 25)
k5_2 = kmeans(scaled_coords2, centers = 5, nstart = 25)
# plots
p1_2 = fviz_cluster(k2_2, geom = "point", data = scaled_coords2) + ggtitle("k = 2")
p2_2 = fviz_cluster(k3_2, geom = "point",  data = scaled_coords2) + ggtitle("k = 3")
p3_2 = fviz_cluster(k4_2, geom = "point",  data = scaled_coords2) + ggtitle("k = 4")
p4_2 = fviz_cluster(k5_2, geom = "point",  data = scaled_coords2) + ggtitle("k = 5")
grid.arrange(p1_2, p2_2, p3_2, p4_2, nrow = 2) 

# k óptimo
fviz_nbclust(scaled_coords2, kmeans, method = "silhouette") # 2 clusters

# Mean shape of each cluster:
n1=nrow(clust1_2)
n2=nrow(clust2_2)
X_Proc_c1 = array(0,dim=c(20,2,n1))
X_Proc_c2 = array(0,dim=c(20,2,n2))
for (i in 1:n1){
  X_Proc_c1[,,i] = cbind(t(clust1_2[i,seq(44,82,2)]),t(clust1_2[i,seq(45,83,2)]))
}
for (i in 1:n2){
  X_Proc_c2[,,i] = cbind(t(clust2_2[i,seq(44,82,2)]),t(clust2_2[i,seq(45,83,2)]))
}


plotAllSpecimens(X_Proc_c1,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                    4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0.5,mean.cex=0.5,link.lwd=0.5))
plotAllSpecimens(X_Proc_c2,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                    4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0.5,mean.cex=0.5,link.lwd=0.5))


# Mean arm shape
C1=mshape(X_Proc_c1)
C2=mshape(X_Proc_c2)

ProcCoords_Age2=array(c(X_Proc_c1, X_Proc_c2), dim=c(20,2,n1+n2))

par(mfrow=c(1,2))
plotAllSpecimens(ProcCoords_Age2,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                          4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0,mean.cex=0.5,link.lwd=0.5))
lines(C1[joinline,],col="lightsteelblue3",lwd=2) # mean shape cluster 1 
title(main="Mean shape cluster 1",cex.main=1.7)


plotAllSpecimens(ProcCoords_Age2,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                          4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0,mean.cex=0.5,link.lwd=0.5))
lines(C2[joinline,],col="lightsteelblue3",lwd=2) # mean shape cluster 2 
title(main="Mean shape cluster 2",cex.main=1.7)


par(mfrow=c(1,1))
levels(Age2$sex)=c('Male','Female')
barplot(with(Age2,prop.table(table(Age2$sex,Age2$cluster),2)),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue3"),cex.lab=1.5,cex.main=2.1,main="Sex")
legend("topright",legend = levels(Age2$sex),bty="n",pch=20,pt.cex=2,cex=1.2,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue3"))

with(Age2,chisq.test(sex,cluster,correct=TRUE,p=rep(1/length(sex),length(sex))))

par(mfrow=c(1,3))
barplot(with(Age2,prop.table(table(Age2$class_MUAC,Age2$cluster),2)),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),cex.lab=1.5,cex.main=2.1,main="Class MUAC")
legend("topright",legend = levels(Age2$class_MUAC),bty="n",pch=20,pt.cex=2,cex=1.2,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age2,prop.table(table(Age2$class_WFL,Age2$cluster),2)),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),cex.lab=1.5,cex.main=2.1,main="Class WFL")
legend("topright",legend = levels(Age2$class_WFL),bty="n",pch=20,pt.cex=2,cex=1.2,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age2,prop.table(table(Age2$class_global,Age2$cluster),2)),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),cex.lab=1.5,cex.main=2.1,main="Class global")
legend("topright",legend = levels(Age2$class_global),bty="n",pch=20,pt.cex=2,cex=1.2,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

par(mfrow=c(2,3))
barplot(with(Age2,prop.table(table(Age2$zlen_cat,Age2$cluster),2)),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),cex.lab=1.5,cex.main=2.1,main="Length for age (zlen)")
barplot(with(Age2,prop.table(table(Age2$zwei_cat,Age2$cluster),2)),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),cex.lab=1.5,cex.main=2.1,main="Weight for age (zwei)")
barplot(with(Age2,prop.table(table(Age2$zwfl_cat,Age2$cluster),2)),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),cex.lab=1.5,cex.main=2.1,main="Weight for length (zwfl)")
barplot(with(Age2,prop.table(table(Age2$zbmi_cat,Age2$cluster),2)),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),cex.lab=1.5,cex.main=2.1,main="BMI for age (zbmi)")
barplot(with(Age2,prop.table(table(Age2$zac_cat,Age2$cluster),2)),beside=TRUE,
        xlab="Cluster", ylab="Frequency",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),cex.lab=1.5,cex.main=2.1,main="MUAC for age (zac)")
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend =levels(data$class_MUAC), pch=16, pt.cex=4, cex=1.8, bty='n',
       col = c("aliceblue","lightblue2","lightblue3","lightblue4"),horiz=FALSE)


df = Age2[,c('cluster','agemons')]
p1=ggplot(df, aes(x=agemons, color=cluster)) + #ggtitle("Age") +
  geom_histogram(fill="white", alpha=0.6, position="identity",bins=15)

df = Age2[,c('cluster','height')]
p2=ggplot(df, aes(x=height, color=cluster)) + ggtitle("Height") +
  geom_histogram(fill="white", alpha=0.6, position="identity",bins=15)

df = Age2[,c('cluster','weight')]
p3=ggplot(df, aes(x=weight, color=cluster)) + ggtitle("Weight") +
  geom_histogram(fill="white", alpha=0.6, position="identity",bins=15)

grid.arrange(p1, p2, p3, ncol = 4, nrow = 2, layout_matrix= rbind(c(NA,1,1,NA),c(2,2,3,3)))

# z-scores
ggplot(Age2, aes(x = zlen)) +
  geom_histogram(fill = "lightblue3", colour = "lightblue4", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age2, aes(x = zwei)) +
  geom_histogram(fill = "lightblue3", colour = "lightblue4", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age2, aes(x = zwfl)) +
  geom_histogram(fill = "lightblue3", colour = "lightblue4", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age2, aes(x = zbmi)) +
  geom_histogram(fill = "lightblue3", colour = "lightblue4", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age2, aes(x = zac)) +
  geom_histogram(fill = "lightblue3", colour = "lightblue4", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)

# Decision tree
par(mfrow=c(1,1))
rpart.fit = rpart(cluster ~., method="class", data = Age2[,-c(c(1:83),87,104:109)]) # [,-c(c(1:43),(44:83),87)]
# summary(rpart.fit)
rpart.plot(rpart.fit, digits = 2, fallen.leaves = TRUE,
           type = 3, extra=101, box.palette = "Blues")

df = data.frame(imp = rpart.fit$variable.importance)
df2 = df %>% 
  tibble::rownames_to_column() %>% 
  dplyr::rename("variable" = rowname) %>% 
  dplyr::arrange(imp) %>%
  dplyr::mutate(variable = forcats::fct_inorder(variable))
ggplot(df2) +
  geom_segment(aes(x = variable, y = 0, xend = variable, yend = imp), 
               size = 1, alpha = 0.7) +
  geom_point(aes(x = variable, y = imp, col = variable), 
             size = 2, show.legend = F) +
  coord_flip() +
  theme_bw()

####################################
### MAHALANOBIS DISTANCE MATRIX ####
####################################

library(pracma)
pseudoinv=pinv(cov(Proc_coords_R))
x0=as.matrix(Proc_coords_R)
mat = apply(x0, 1, function(i) mahalanobis(x0, i, cov = pseudoinv,inverted=TRUE))
mat=sqrt(mat)

########################
### PAM  CLUSTERING ####
########################

pam_2 = pam(mat,k=2,diss = TRUE) # best option
pam_3 = pam(mat,k=3,diss = TRUE)
pam_4 = pam(mat,k=4,diss = TRUE)
pam_5 = pam(mat,k=5,diss = TRUE)


par(mfrow=c(2,2))
# silhouettes
sil_pam_2 = silhouette(pam_2$cluster,mat)
plot(sil_pam_2,col='blue',border=NA)
sil_pam_3 = silhouette(pam_3$cluster,mat)
plot(sil_pam_3,col='blue',border=NA)
sil_pam_4 = silhouette(pam_4$cluster,mat)
plot(sil_pam_4,col='blue',border=NA)
sil_pam_5 = silhouette(pam_5$cluster,mat)
plot(sil_pam_5,col='blue',border=NA)

data[,'cluster']= pam_2$clustering 

round(prop.table(table(data$class_MUAC,data$cluster),2),3)
round(prop.table(table(data$class_WFL,data$cluster),2),3)
round(prop.table(table(data$class_global,data$cluster),2),3)

round(prop.table(table(data$zlen_cat,data$cluster),2),3)
round(prop.table(table(data$zwei_cat,data$cluster),2),3)
round(prop.table(table(data$zwfl_cat,data$cluster),2),3)
round(prop.table(table(data$zbmi_cat,data$cluster),2),3)
round(prop.table(table(data$zac_cat,data$cluster),2),3)


# datos del custer 1 vs cluster 2 por separado:
clust1_PAM=data[pam_2$clustering==1,]
clust2_PAM=data[pam_2$clustering==2,]

summary(clust1_PAM[,c(84:88)])
summary(clust2_PAM[,c(84:88)])

summary(clust1_PAM[,c(89:93)])
summary(clust2_PAM[,c(89:93)])

summary(clust1_PAM[,c(94:109)])
summary(clust2_PAM[,c(94:109)])

table(clust1_PAM$sex)
table(clust2_PAM$sex)

# sex
with(data,chisq.test(sex,cluster,correct=TRUE,p=rep(1/length(sex),length(sex))))
# class indicators
with(data,chisq.test(class_MUAC,cluster,correct=TRUE,p=rep(1/length(class_MUAC),length(class_MUAC))))
with(data,chisq.test(class_WFL,cluster,correct=TRUE,p=rep(1/length(class_WFL),length(class_WFL))))
with(data,chisq.test(class_global,cluster,correct=TRUE,p=rep(1/length(class_global),length(class_global))))
# z-scores
with(data,chisq.test(zwei_cat,cluster,correct=TRUE,p=rep(1/length(zwei_cat),length(zwei_cat))))
with(data,chisq.test(zlen_cat,cluster,correct=TRUE,p=rep(1/length(zlen_cat),length(zlen_cat))))
with(data,chisq.test(zwfl_cat,cluster,correct=TRUE,p=rep(1/length(zwfl_cat),length(zwfl_cat))))
with(data,chisq.test(zbmi_cat,cluster,correct=TRUE,p=rep(1/length(zbmi_cat),length(zbmi_cat))))
with(data,chisq.test(zac_cat,cluster,correct=TRUE,p=rep(1/length(zac_cat),length(zac_cat))))
# measure
with(data,chisq.test(measure,cluster,correct=TRUE,p=rep(1/length(measure),length(measure))))
# age
with(data,chisq.test(age,cluster,correct=TRUE,p=rep(1/length(age),length(age))))
# all p-values higher than 0.05!


C1=mshape(X_Proc_R$rotated[,,pam_2$cluster==1])
C2=mshape(X_Proc_R$rotated[,,pam_2$cluster==2])
joinline=c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1)
par(mfrow=c(1,2))
#,link.lty="dotted"
plotAllSpecimens(X_Proc_R$rotated,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                           4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0,mean.cex=0.5,link.lwd=0.5))
lines(C1[joinline,],col="lightsteelblue3",lwd=2) # mean shape cluster 1 
title(main="Mean shape cluster 1",cex.main=1.7)


plotAllSpecimens(X_Proc_R$rotated,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                           4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0,mean.cex=0.5,link.lwd=0.5))
lines(C2[joinline,],col="lightsteelblue3",lwd=2) # mean shape cluster 2 
title(main="Mean shape cluster 2",cex.main=1.7)


# Factor variables by cluster
(table(clust1_PAM$sex))
(table(clust2_PAM$sex))

(table(clust1_PAM$measure))
(table(clust2_PAM$measure))

(table(clust1_PAM$class_MUAC))
(table(clust2_PAM$class_MUAC))

(table(clust1_PAM$class_WFL))
(table(clust2_PAM$class_WFL))

(table(clust1_PAM$class_global))
(table(clust2_PAM$class_global))

(table(clust1_PAM$age))
(table(clust2_PAM$age))

(table(clust1_PAM$zlen_cat))
(table(clust2_PAM$zlen_cat))

(table(clust1_PAM$zwei_cat))
(table(clust2_PAM$zwei_cat))

(table(clust1_PAM$zwfl_cat))
(table(clust2_PAM$zwfl_cat))

(table(clust1_PAM$zbmi_cat))
(table(clust2_PAM$zbmi_cat))

(table(clust1_PAM$zac_cat))
(table(clust2_PAM$zac_cat))


