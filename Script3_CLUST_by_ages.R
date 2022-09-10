###########################
### CLUSTERING BY AGES ####
###########################

# age dataframes
Age1=data[agemons<=12,]
Age2=data[agemons>12&agemons<=24,]
Age3=data[agemons>24&agemons<=36,]
Age4=data[agemons>36&agemons<=48,]
Age5=data[agemons>48,]


# AGE<1 ##################
scaled_coords1=scale(Age1[,c(44:83)],center = TRUE,scale = TRUE)

# values of k
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

# optimal k
fviz_nbclust(scaled_coords1, kmeans, method = "silhouette") # 2 clusters

Age1[,'cluster']=k2_1$cluster
Age1$cluster=as.factor(Age1$cluster)


clust1_1=Age1[k2_1$cluster==1,]
clust2_1=Age1[k2_1$cluster==2,]

v = seq(44,83,by=2)
vv = seq(45,83,by=2)

par(mfrow=c(1,2))
coord_x_c1 = as.matrix(clust1_1[,v]) # x
coord_y_c1 = as.matrix(clust1_1[,vv]) # y
plot(coord_x_c1,coord_y_c1,main="Coords cluster 1")

coord_x_c2 = as.matrix(clust2_1[,v]) # x
coord_y_c2 = as.matrix(clust2_1[,vv]) # y
plot(coord_x_c2,coord_y_c2,main="Coords cluster 2")


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
                                                    4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),
                 plot.param = list(pt.cex=0.5,mean.cex=0.5,link.lwd=0.5))

plotAllSpecimens(X_Proc_c2,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                    4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),
                 plot.param = list(pt.cex=0.5,mean.cex=0.5,link.lwd=0.5))


# mean arm shape for each cluster

C1=mshape(X_Proc_c1)
C2=mshape(X_Proc_c2)

ProcCoords_Age1=array(c(X_Proc_c1, X_Proc_c2), dim=c(20,2,n1+n2))

par(mfrow=c(1,2))
plotAllSpecimens(ProcCoords_Age1,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                          4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0,mean.cex=0.5,link.lwd=0.5))
lines(C1[joinline,],col="lightsteelblue3",lwd=2) # mean shape cluster 1 
title(main="Mean shape cluster 1",cex.main=1.5)


plotAllSpecimens(ProcCoords_Age1,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                          4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0,mean.cex=0.5,link.lwd=0.5))
lines(C2[joinline,],col="lightsteelblue3",lwd=2) # mean shape cluster 1 
title(main="Mean shape cluster 2",cex.main=1.5)



par(mfrow=c(1,1))
par(xpd=TRUE) 

levels(Age1$sex)=c('Male','Female')
barplot(with(Age1,(table(Age1$sex,Age1$clust))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue3"),main="Sex")
legend("topright",legend = levels(Age1$sex),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue3"))

(table(clust1_1$sex))
(table(clust2_1$sex))

barplot(with(Age1,(table(Age1$class_MUAC,Age1$clust))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="Class MUAC")
legend("topright",legend = levels(Age1$class_MUAC),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age1,(table(Age1$class_WFL,Age1$clust))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"), main="Class WFL")
legend("topright",legend = levels(Age1$class_WFL),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age1,(table(Age1$class_global,Age1$clust))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="Class global")
legend("topright",legend = levels(Age1$class_global),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age1,(table(Age1$zlen_cat,Age1$clust))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZLEN")
legend("topright",legend = levels(Age1$zlen_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age1,(table(Age1$zwei_cat,Age1$clust))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZWEI")
legend("topright",legend = levels(Age1$zwei_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age1,(table(Age1$zbmi_cat,Age1$clust))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZBMI")
legend("topright",legend = levels(Age1$zbmi_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age1,(table(Age1$zwfl_cat,Age1$clust))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZWFL")
legend("topright",legend = levels(Age1$zwfl_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age1,(table(Age1$zac_cat,Age1$clust))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZAC")
legend("topright",legend = levels(Age1$zac_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))


# histograms of continuous variables
ggplot(Age1, aes(x = agemons)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age1, aes(x = height)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age1, aes(x = weight)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age1, aes(x = muac)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)

# z-scores
ggplot(Age1, aes(x = zlen)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age1, aes(x = zwei)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age1, aes(x = zwfl)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age1, aes(x = zbmi)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age1, aes(x = zac)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)


# Decision tree
par(mfrow=c(1,1))
rpart.fit = rpart(cluster ~., method="class", data = Age1[,-c(1:83,87)]) # -c(1:83,87,98:103)
# summary(rpart.fit)
rpart.plot(rpart.fit, digits = 3, fallen.leaves = TRUE,
           type = 3, extra=101)

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


# AGE<2 ##################
scaled_coords2=scale(Age2[,c(44:83)],center = TRUE,scale = TRUE)

# values of k
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

# optimal k
fviz_nbclust(scaled_coords2, kmeans, method = "silhouette") # 4 clusters -> ver si son los grupos nutricionales!


clust1_2=Age2[k4_2$cluster==1,]
clust2_2=Age2[k4_2$cluster==2,]
clust3_2=Age2[k4_2$cluster==3,]
clust4_2=Age2[k4_2$cluster==4,]


Age2[,'cluster']=k4_2$cluster 
Age2$cluster=as.factor(Age2$cluster)

round(prop.table(table(Age2$class_MUAC,Age2$cluster),2),3)
round(prop.table(table(Age2$class_WFL,Age2$cluster),2),3)
round(prop.table(table(Age2$class_global,Age2$cluster),2),3)


v = seq(44,83,by=2)
vv = seq(45,83,by=2)

par(mfrow=c(2,2))
coord_x_c1 = as.matrix(clust1_2[,v]) # x
coord_y_c1 = as.matrix(clust1_2[,vv]) # y
plot(coord_x_c1,coord_y_c1,main="Coords cluster 1")

coord_x_c2 = as.matrix(clust2_2[,v]) # x
coord_y_c2 = as.matrix(clust2_2[,vv]) # y
plot(coord_x_c2,coord_y_c2,main="Coords cluster 2")

coord_x_c3 = as.matrix(clust3_2[,v]) # x
coord_y_c3 = as.matrix(clust3_2[,vv]) # y
plot(coord_x_c3,coord_y_c3,main="Coords cluster 3")

coord_x_c4 = as.matrix(clust4_2[,v]) # x
coord_y_c4 = as.matrix(clust4_2[,vv]) # y
plot(coord_x_c4,coord_y_c4,main="Coords cluster 4")


# Mean shape of each cluster:
n1=nrow(clust1_2)
n2=nrow(clust2_2)
n3=nrow(clust3_2)
n4=nrow(clust4_2)
X_Proc_c1 = array(0,dim=c(20,2,n1))
X_Proc_c2 = array(0,dim=c(20,2,n2))
X_Proc_c3 = array(0,dim=c(20,2,n3))
X_Proc_c4 = array(0,dim=c(20,2,n4))
for (i in 1:n1){
  X_Proc_c1[,,i] = cbind(t(clust1_2[i,seq(44,82,2)]),t(clust1_2[i,seq(45,83,2)]))
}
for (i in 1:n2){
  X_Proc_c2[,,i] = cbind(t(clust2_2[i,seq(44,82,2)]),t(clust2_2[i,seq(45,83,2)]))
}
for (i in 1:n3){
  X_Proc_c3[,,i] = cbind(t(clust3_2[i,seq(44,82,2)]),t(clust3_2[i,seq(45,83,2)]))
}
for (i in 1:n4){
  X_Proc_c4[,,i] = cbind(t(clust4_2[i,seq(44,82,2)]),t(clust4_2[i,seq(45,83,2)]))
}


plotAllSpecimens(X_Proc_c1,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                    4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0.5,mean.cex=0.5,link.lwd=0.5))
plotAllSpecimens(X_Proc_c2,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                    4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0.5,mean.cex=0.5,link.lwd=0.5))
plotAllSpecimens(X_Proc_c3,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                    4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0.5,mean.cex=0.5,link.lwd=0.5))
plotAllSpecimens(X_Proc_c4,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                    4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0.5,mean.cex=0.5,link.lwd=0.5))



# mean arm shape for each cluster
C1=mshape(X_Proc_c1)
C2=mshape(X_Proc_c2)
C3=mshape(X_Proc_c3)
C4=mshape(X_Proc_c4)

ProcCoords_Age2=array(c(X_Proc_c1, X_Proc_c2,X_Proc_c3, X_Proc_c4), dim=c(20,2,n1+n2+n3+n4))
par(mfrow=c(2,2))
plotAllSpecimens(ProcCoords_Age2,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                          4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0,mean.cex=0.5,link.lwd=0.5))
lines(C1[joinline,],col="lightsteelblue3",lwd=2) # mean shape cluster 1 
title(main="Mean shape cluster 1",cex.main=1.5)


plotAllSpecimens(ProcCoords_Age2,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                          4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0,mean.cex=0.5,link.lwd=0.5))
lines(C2[joinline,],col="lightsteelblue3",lwd=2) # mean shape cluster 1 
title(main="Mean shape cluster 2",cex.main=1.5)

plotAllSpecimens(ProcCoords_Age2,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                          4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0,mean.cex=0.5,link.lwd=0.5))
lines(C3[joinline,],col="lightsteelblue3",lwd=2) # mean shape cluster 3
title(main="Mean shape cluster 3")

plotAllSpecimens(ProcCoords_Age2,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                          4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0,mean.cex=0.5,link.lwd=0.5))
lines(C4[joinline,],col="lightsteelblue3",lwd=2) # mean shape cluster 4
title(main="Mean shape cluster 4")




par(mfrow=c(1,1))
levels(Age2$sex)=c('Male','Female')
barplot(with(Age2,(table(Age2$sex,Age2$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue3"),main="Sex")
legend("topright",legend = levels(Age2$sex),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue3"))

par(mfrow=c(2,3))
barplot(with(Age2,(table(Age2$class_MUAC,Age2$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="Class MUAC")
legend("topright",legend = levels(Age2$class_MUAC),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age2,(table(Age2$class_WFL,Age2$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="Class WFL")
legend("topright",legend = levels(Age2$class_WFL),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age2,(table(Age2$class_global,Age2$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="Class global")
legend("topright",legend = levels(Age2$class_global),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age2,(table(Age2$zlen_cat,Age2$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZLEN")
legend("topright",legend = levels(Age2$zlen_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age2,(table(Age2$zwei_cat,Age2$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZWEI")
legend("topright",legend = levels(Age2$zwei_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age2,(table(Age2$zbmi_cat,Age2$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZBMI")
legend("topright",legend = levels(Age2$zbmi_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age2,(table(Age2$zwfl_cat,Age2$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZWFL")
legend("topright",legend = levels(Age2$zwfl_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age2,(table(Age2$zac_cat,Age2$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZAC")
legend("topright",legend = levels(Age2$zac_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

# histograms of continuous variables
ggplot(Age2, aes(x = agemons)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age2, aes(x = height)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age2, aes(x = weight)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age2, aes(x = muac)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)

# z-scores
ggplot(Age2, aes(x = zlen)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age2, aes(x = zwei)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age2, aes(x = zwfl)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age2, aes(x = zbmi)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age2, aes(x = zac)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)

# Decision tree
par(mfrow=c(1,1))
rpart.fit = rpart(cluster ~., method="class", data = Age2[,-c(1:83,87)]) # -c(1:83,87,98:103)
# summary(rpart.fit)
rpart.plot(rpart.fit, digits = 3, fallen.leaves = TRUE,
           type = 3, extra=101)

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


# AGE<3 ##################
scaled_coords3=scale(Age3[,c(44:83)],center = TRUE,scale = TRUE)

# values of k
k2_3 = kmeans(scaled_coords3, centers = 2, nstart = 25)
k3_3 = kmeans(scaled_coords3, centers = 3, nstart = 25)
k4_3 = kmeans(scaled_coords3, centers = 4, nstart = 25)
k5_3 = kmeans(scaled_coords3, centers = 5, nstart = 25)

# plots
p1_3 = fviz_cluster(k2_3, geom = "point", data = scaled_coords3) + ggtitle("k = 2")
p2_3 = fviz_cluster(k3_3, geom = "point",  data = scaled_coords3) + ggtitle("k = 3")
p3_3 = fviz_cluster(k4_3, geom = "point",  data = scaled_coords3) + ggtitle("k = 4")
p4_3 = fviz_cluster(k5_3, geom = "point",  data = scaled_coords3) + ggtitle("k = 5")
grid.arrange(p1_3, p2_3, p3_3, p4_3, nrow = 2) 

# optimal k
fviz_nbclust(scaled_coords3, kmeans, method = "silhouette") # 4 clusters -> ver si son los grupos nutricionales!

clust1_3=Age3[k4_3$cluster==1,]
clust2_3=Age3[k4_3$cluster==2,]
clust3_3=Age3[k4_3$cluster==3,]
clust4_3=Age3[k4_3$cluster==4,]

Age3[,'cluster']=k4_3$cluster
Age3$cluster=as.factor(Age3$cluster)


v = seq(44,83,by=2)
vv = seq(45,83,by=2)

par(mfrow=c(2,2))
coord_x_c1 = as.matrix(clust1_3[,v]) # x
coord_y_c1 = as.matrix(clust1_3[,vv]) # y
plot(coord_x_c1,coord_y_c1,main="Coords cluster 1")

coord_x_c2 = as.matrix(clust2_3[,v]) # x
coord_y_c2 = as.matrix(clust2_3[,vv]) # y
plot(coord_x_c2,coord_y_c2,main="Coords cluster 2")

coord_x_c3 = as.matrix(clust3_3[,v]) # x
coord_y_c3 = as.matrix(clust3_3[,vv]) # y
plot(coord_x_c3,coord_y_c3,main="Coords cluster 3")

coord_x_c4 = as.matrix(clust4_3[,v]) # x
coord_y_c4 = as.matrix(clust4_3[,vv]) # y
plot(coord_x_c4,coord_y_c4,main="Coords cluster 4")


# Mean shape of each cluster:
n1=nrow(clust1_3)
n2=nrow(clust2_3)
n3=nrow(clust3_3)
n4=nrow(clust4_3)
X_Proc_c1 = array(0,dim=c(20,2,n1))
X_Proc_c2 = array(0,dim=c(20,2,n2))
X_Proc_c3 = array(0,dim=c(20,2,n3))
X_Proc_c4 = array(0,dim=c(20,2,n4))
for (i in 1:n1){
  X_Proc_c1[,,i] = cbind(t(clust1_3[i,seq(44,82,2)]),t(clust1_3[i,seq(45,83,2)]))
}
for (i in 1:n2){
  X_Proc_c2[,,i] = cbind(t(clust2_3[i,seq(44,82,2)]),t(clust2_3[i,seq(45,83,2)]))
}
for (i in 1:n3){
  X_Proc_c3[,,i] = cbind(t(clust3_3[i,seq(44,82,2)]),t(clust3_3[i,seq(45,83,2)]))
}
for (i in 1:n4){
  X_Proc_c4[,,i] = cbind(t(clust4_3[i,seq(44,82,2)]),t(clust4_3[i,seq(45,83,2)]))
}
par(mfrow=c(2,2))

plotAllSpecimens(X_Proc_c1,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                    4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0.5,mean.cex=0.5,link.lwd=0.5))
plotAllSpecimens(X_Proc_c2,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                    4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0.5,mean.cex=0.5,link.lwd=0.5))
plotAllSpecimens(X_Proc_c3,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                    4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0.5,mean.cex=0.5,link.lwd=0.5))
plotAllSpecimens(X_Proc_c4,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                    4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0.5,mean.cex=0.5,link.lwd=0.5))


# mean arm shape for each cluster
C1=mshape(X_Proc_c1)
C2=mshape(X_Proc_c2)
C3=mshape(X_Proc_c3)
C4=mshape(X_Proc_c4)

ProcCoords_Age3=array(c(X_Proc_c1, X_Proc_c2,X_Proc_c3, X_Proc_c4), dim=c(20,2,n1+n2+n3+n4))
joinline=c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1)
par(mfrow=c(2,2))
plotAllSpecimens(ProcCoords_Age3,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                          4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0,mean.cex=0.5,link.lwd=0.5))
lines(C1[joinline,],col="lightsteelblue3",lwd=2) # mean shape cluster 1 
title(main="Mean shape cluster 1",cex.main=1.5)


plotAllSpecimens(ProcCoords_Age3,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                          4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0,mean.cex=0.5,link.lwd=0.5))
lines(C2[joinline,],col="lightsteelblue3",lwd=2) # mean shape cluster 1 
title(main="Mean shape cluster 2",cex.main=1.5)

plotAllSpecimens(ProcCoords_Age3,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                          4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0,mean.cex=0.5,link.lwd=0.5))
lines(C3[joinline,],col="lightsteelblue3",lwd=2) # mean shape cluster 3
title(main="Mean shape cluster 3")

plotAllSpecimens(ProcCoords_Age3,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                          4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0,mean.cex=0.5,link.lwd=0.5))
lines(C4[joinline,],col="lightsteelblue3",lwd=2) # mean shape cluster 4
title(main="Mean shape cluster 4")




par(mfrow=c(1,1))
levels(Age3$sex)=c('Male','Female')
barplot(with(Age3,(table(Age3$sex,Age3$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue3"),main="Sex")
legend("topright",legend = levels(Age3$sex),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue3"))

barplot(with(Age3,(table(Age3$class_MUAC,Age3$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="Class MUAC")
legend("topright",legend = levels(Age3$class_MUAC),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age3,(table(Age3$class_WFL,Age3$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="Class WFL")
legend("topright",legend = levels(Age3$class_WFL),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age3,(table(Age3$class_global,Age3$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="Class global")
legend("topright",legend = levels(Age3$class_global),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age3,(table(Age3$zlen_cat,Age3$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZLEN")
legend("topright",legend = levels(Age3$zlen_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age3,(table(Age3$zwei_cat,Age3$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZWEI")
legend("topright",legend = levels(Age3$zwei_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age3,(table(Age3$zbmi_cat,Age3$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZBMI")
legend("topright",legend = levels(Age3$zbmi_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age3,(table(Age3$zwfl_cat,Age3$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZWFL")
legend("topright",legend = levels(Age3$zwfl_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age3,(table(Age3$zac_cat,Age3$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZAC")
legend("topright",legend = levels(Age3$zac_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))


# histograms of continuous variables
ggplot(Age3, aes(x = agemons)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age3, aes(x = height)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age3, aes(x = weight)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age3, aes(x = muac)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)


# z-scores
ggplot(Age3, aes(x = zlen)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age3, aes(x = zwei)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age3, aes(x = zwfl)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age3, aes(x = zbmi)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age3, aes(x = zac)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)


# Decision tree
par(mfrow=c(1,1))
rpart.fit = rpart(cluster ~., method="class", data = Age3[,-c(c(1:83),87)]) # [,-c(c(1:83),87)]
# summary(rpart.fit)
rpart.plot(rpart.fit, digits = 3, fallen.leaves = TRUE,
           type = 3, extra=101)

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


# AGE<4 ##################
scaled_coords4=scale(Age4[,c(44:83)],center = TRUE,scale = TRUE)

# values of k
k2_4 = kmeans(scaled_coords4, centers = 2, nstart = 25)
k3_4 = kmeans(scaled_coords4, centers = 3, nstart = 25)
k4_4 = kmeans(scaled_coords4, centers = 4, nstart = 25)
k5_4 = kmeans(scaled_coords4, centers = 5, nstart = 25)

# plots
p1_4 = fviz_cluster(k2_4, geom = "point", data = scaled_coords4) + ggtitle("k = 2")
p2_4 = fviz_cluster(k3_4, geom = "point",  data = scaled_coords4) + ggtitle("k = 3")
p3_4 = fviz_cluster(k4_4, geom = "point",  data = scaled_coords4) + ggtitle("k = 4")
p4_4 = fviz_cluster(k5_4, geom = "point",  data = scaled_coords4) + ggtitle("k = 5")
grid.arrange(p1_4, p2_4, p3_4, p4_4, nrow = 2) 

# optimal k
fviz_nbclust(scaled_coords4, kmeans, method = "silhouette") # 2 clusters 

clust1_4=Age4[k2_4$cluster==1,]
clust2_4=Age4[k2_4$cluster==2,]

Age4[,'cluster']=k2_4$cluster
Age4$cluster=as.factor(Age4$cluster)


v = seq(44,83,by=2)
vv = seq(45,83,by=2)

par(mfrow=c(1,2))
coord_x_c1 = as.matrix(clust1_4[,v]) # x
coord_y_c1 = as.matrix(clust1_4[,vv]) # y
plot(coord_x_c1,coord_y_c1,main="Coords cluster 1")

coord_x_c2 = as.matrix(clust2_4[,v]) # x
coord_y_c2 = as.matrix(clust2_4[,vv]) # y
plot(coord_x_c2,coord_y_c2,main="Coords cluster 2")


# Mean shape of each cluster:
n1=nrow(clust1_4)
n2=nrow(clust2_4)
X_Proc_c1 = array(0,dim=c(20,2,n1))
X_Proc_c2 = array(0,dim=c(20,2,n2))
for (i in 1:n1){
  X_Proc_c1[,,i] = cbind(t(clust1_4[i,seq(44,82,2)]),t(clust1_4[i,seq(45,83,2)]))
}

for (i in 1:n2){
  X_Proc_c2[,,i] = cbind(t(clust2_4[i,seq(44,82,2)]),t(clust2_4[i,seq(45,83,2)]))
}


plotAllSpecimens(X_Proc_c1,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                    4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0.5,mean.cex=0.5,link.lwd=0.5))

plotAllSpecimens(X_Proc_c2,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                    4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0.5,mean.cex=0.5,link.lwd=0.5))


# mean arm shape for each cluster
C1=mshape(X_Proc_c1)
C2=mshape(X_Proc_c2)

ProcCoords_Age4=array(c(X_Proc_c1, X_Proc_c2), dim=c(20,2,n1+n2))

par(mfrow=c(1,2))
plotAllSpecimens(ProcCoords_Age4,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                          4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0,mean.cex=0.5,link.lwd=0.5))
lines(C1[joinline,],col="lightsteelblue3",lwd=2) # mean shape cluster 1 
title(main="Mean shape cluster 1",cex.main=1.5)


plotAllSpecimens(ProcCoords_Age4,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                          4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0,mean.cex=0.5,link.lwd=0.5))
lines(C2[joinline,],col="lightsteelblue3",lwd=2) # mean shape cluster 2 
title(main="Mean shape cluster 2",cex.main=1.5)




par(mfrow=c(1,1))
levels(Age4$sex)=c('Male','Female')
barplot(with(Age4,(table(Age4$sex,Age4$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue3"),main="Sex")
legend("topright",legend = levels(Age4$sex),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue3"))

barplot(with(Age4,(table(Age4$class_MUAC,Age4$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="Class MUAC")
legend("topright",legend = levels(Age4$class_MUAC),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age4,(table(Age4$class_WFL,Age4$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="Class WFL")
legend("topright",legend = levels(Age4$class_WFL),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age4,(table(Age4$class_global,Age4$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="Class global")
legend("topright",legend = levels(Age4$class_global),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age4,(table(Age4$zlen_cat,Age4$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZLEN")
legend("topright",legend = levels(Age4$zlen_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age4,(table(Age4$zwei_cat,Age4$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZWEI")
legend("topright",legend = levels(Age4$zwei_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age4,(table(Age4$zbmi_cat,Age4$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZBMI")
legend("topright",legend = levels(Age4$zbmi_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age4,(table(Age4$zwfl_cat,Age4$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZWFL")
legend("topright",legend = levels(Age4$zwfl_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age4,(table(Age4$zac_cat,Age4$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZAC")
legend("topright",legend = levels(Age4$zac_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))


# histograms of continuous variables
ggplot(Age4, aes(x = agemons)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age4, aes(x = height)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age4, aes(x = weight)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age4, aes(x = muac)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)

# z-scores
ggplot(Age4, aes(x = zlen)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age4, aes(x = zwei)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age4, aes(x = zwfl)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age4, aes(x = zbmi)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age4, aes(x = zac)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)

# Decision tree
par(mfrow=c(1,1))
rpart.fit = rpart(cluster ~., method="class", data = Age4[,-c(c(1:83),87)]) # [,-c(c(1:43),(44:83),87)]
# summary(rpart.fit)
rpart.plot(rpart.fit, digits = 3, fallen.leaves = TRUE,
           type = 3, extra=101)

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


# AGE<5 ##################
scaled_coords5=scale(Age5[,c(44:83)],center = TRUE,scale = TRUE)

# values of k
k2_5 = kmeans(scaled_coords5, centers = 2, nstart = 25)
k3_5 = kmeans(scaled_coords5, centers = 3, nstart = 25)
k4_5 = kmeans(scaled_coords5, centers = 4, nstart = 25)
k5_5 = kmeans(scaled_coords5, centers = 5, nstart = 25)

# plots
p1_5 = fviz_cluster(k2_5, geom = "point", data = scaled_coords5) + ggtitle("k = 2")
p2_5 = fviz_cluster(k3_5, geom = "point",  data = scaled_coords5) + ggtitle("k = 3")
p3_5 = fviz_cluster(k4_5, geom = "point",  data = scaled_coords5) + ggtitle("k = 4")
p4_5 = fviz_cluster(k5_5, geom = "point",  data = scaled_coords5) + ggtitle("k = 5")
grid.arrange(p1_5, p2_5, p3_5, p4_5, nrow = 2) 

# optimal k
fviz_nbclust(scaled_coords5, kmeans, method = "silhouette") # 2 clusters

Age5[,'cluster']=k2_5$cluster
Age5$cluster=as.factor(Age5$cluster)


clust1_5=Age5[k2_5$cluster==1,]
clust2_5=Age5[k2_5$cluster==2,]


v = seq(44,83,by=2)
vv = seq(45,83,by=2)

par(mfrow=c(1,2))
coord_x_c1 = as.matrix(clust1_5[,v]) # x
coord_y_c1 = as.matrix(clust1_5[,vv]) # y
plot(coord_x_c1,coord_y_c1,main="Coords cluster 1")

coord_x_c2 = as.matrix(clust2_5[,v]) # x
coord_y_c2 = as.matrix(clust2_5[,vv]) # y
plot(coord_x_c2,coord_y_c2,main="Coords cluster 2")


# Mean shape of each cluster:
n1=nrow(clust1_5)
n2=nrow(clust2_5)
X_Proc_c1 = array(0,dim=c(20,2,n1))
X_Proc_c2 = array(0,dim=c(20,2,n2))
for (i in 1:n1){
  X_Proc_c1[,,i] = cbind(t(clust1_5[i,seq(44,82,2)]),t(clust1_5[i,seq(45,83,2)]))
}

for (i in 1:n2){
  X_Proc_c2[,,i] = cbind(t(clust2_5[i,seq(44,82,2)]),t(clust2_5[i,seq(45,83,2)]))
}



plotAllSpecimens(X_Proc_c1,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                    4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0.5,mean.cex=0.5,link.lwd=0.5))
plotAllSpecimens(X_Proc_c2,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                    4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0.5,mean.cex=0.5,link.lwd=0.5))


# mean arm shape for each cluster
C1=mshape(X_Proc_c1)
C2=mshape(X_Proc_c2)

ProcCoords_Age5=array(c(X_Proc_c1, X_Proc_c2), dim=c(20,2,n1+n2))

par(mfrow=c(1,2))
plotAllSpecimens(ProcCoords_Age5,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                          4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0,mean.cex=0.5,link.lwd=0.5))
lines(C1[joinline,],col="lightsteelblue3",lwd=2) # mean shape cluster 1 
title(main="Mean shape cluster 1",cex.main=1.5)

plotAllSpecimens(ProcCoords_Age5,mean=TRUE,links=matrix(c(1,4,16,14,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,
                                                          4,14,16,12,18,10,8,6,20,3,19,5,7,9,17,11,13,15,1),ncol=2),plot.param = list(pt.cex=0,mean.cex=0.5,link.lwd=0.5))
lines(C2[joinline,],col="lightsteelblue3",lwd=2) # mean shape cluster 1 
title(main="Mean shape cluster 2",cex.main=1.5)




par(mfrow=c(1,1))
levels(Age5$sex)=c('Male','Female')
barplot(with(Age5,(table(Age5$sex,Age5$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue3"),main="Sex")
legend("topright",legend = levels(Age5$sex),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue3"))

with(Age5,chisq.test(sex,cluster,correct=TRUE,p=rep(1/length(sex),length(sex))))

(table(clust1_5$class_global))
(table(clust2_5$class_global))

barplot(with(Age5,(table(Age5$class_MUAC,Age5$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="Class MUAC")
legend("topright",legend = levels(Age5$class_MUAC),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age5,(table(Age5$class_WFL,Age5$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="Class WFL")
legend("topright",legend = levels(Age5$class_WFL),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age5,(table(Age5$class_global,Age5$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="Class global")
legend("topright",legend = levels(Age5$class_global),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age5,(table(Age5$zlen_cat,Age5$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZLEN")
legend("topright",legend = levels(Age5$zlen_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age5,(table(Age5$zwei_cat,Age5$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZWEI")
legend("topright",legend = levels(Age5$zwei_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age5,(table(Age5$zwfl_cat,Age5$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZWFL")
legend("topright",legend = levels(Age5$zwfl_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age5,(table(Age5$zbmi_cat,Age5$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZBMI")
legend("topright",legend = levels(Age5$zbmi_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

barplot(with(Age5,(table(Age5$zac_cat,Age5$cluster))),beside=TRUE,
        xlab="Cluster", ylab="Frecuencia",
        col=c("aliceblue","lightblue2","lightblue3","lightblue4"),main="ZAC")
legend("topright",legend = levels(Age5$zac_cat),bty="n",pch=20,pt.cex=2,cex=0.8,horiz=FALSE,inset=c(0.0,0.0),#inset=c(-0.2,0),
       col=c("aliceblue","lightblue2","lightblue3","lightblue4"))

# histograms of continuous variables
ggplot(Age5, aes(x = agemons)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age5, aes(x = height)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age5, aes(x = weight)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age5, aes(x = muac)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)

# z-scores
ggplot(Age5, aes(x = zlen)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age5, aes(x = zwei)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age5, aes(x = zwfl)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age5, aes(x = zbmi)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)
ggplot(Age5, aes(x = zac)) +
  geom_histogram(fill = "pink", colour = "brown", size = .3,bins=15) +
  scale_y_continuous(name = "Number of observations") +
  scale_x_continuous(name = "Cluster") +
  facet_wrap(~cluster)


# Decision tree
par(mfrow=c(1,1))
rpart.fit = rpart(cluster ~., method="class", data = Age5[,-c(c(1:83),87)]) # [,-c(c(1:43),(44:83),87)]
# summary(rpart.fit)
rpart.plot(rpart.fit, digits = 3, fallen.leaves = TRUE,
           type = 3, extra=101)

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
