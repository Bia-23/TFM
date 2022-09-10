### PRINCIPAL COMPONENT ANALYSIS ####

pca = prcomp(Proc_coords_R, scale=T)
pca
summary(pca) 
fviz_eig(pca,addlabels = TRUE)

par(mfrow=c(1,2))
plot(pca$x,main="PC1 vs PC2",cex.lab=1.3,cex.main=1.7)
plot(pca$x[,2:3],main="PC2 vs PC3",cex.lab=1.3,cex.main=1.7)

pairs(pca$x[,1:5],col="steelblue",pch=20,main="First five PCs")

# PCs in terms of other variables
pairs(pca$x[,1:5],col=measure,pch=20,main="First five PCs by measure")
pairs(pca$x[,1:5],col=sex,pch=20,main="First five PCs by sex")
pairs(pca$x[,1:5],col=age,pch=20,main="First five PCs by age groups")
pairs(pca$x[,1:5],col=class_global,pch=20,main="First five PCs by class global")

par(mfrow=c(2,4))
par(xpd=TRUE) 

#### PC1 vs PC2 ####
with(data,plot(pca$x,col=class_MUAC,main="PCs by class MUAC",cex.main=1.7))
with(data,plot(pca$x,col=class_WFL,main="PCs by class WFL",cex.main=1.7))
with(data,plot(pca$x,col=class_global,main="PCs by class global",cex.main=1.7))
with(data,plot(pca$x,col=age,main="PCs by age",cex.main=1.7))

with(data,plot(pca$x[,2:3],col=class_MUAC,main="PCs by class MUAC",cex.main=1.7))
with(data,plot(pca$x[,2:3],col=class_WFL,main="PCs by class WFL",cex.main=1.7))
with(data,plot(pca$x[,2:3],col=class_global,main="PCs by class global",cex.main=1.7))
with(data,plot(pca$x[,2:3],col=age,main="PCs by age",cex.main=1.7))


with(data,plot(pca$x,col=class_MUAC,main="PCs by class MUAC",cex.main=1.7))
legend("bottomleft", legend=levels(class_MUAC), pch=15, col=c(4,3,2,1)) 
with(data,plot(pca$x,col=class_MUAC_3,main="PCs by class MUAC - 3 groups"))
legend("bottomleft", legend=levels(class_MUAC_3), pch=15, col=c(3,2,1)) 
with(data,plot(pca$x,col=class_MUAC_2,main="PCs by class MUAC - 2 groups"))
legend("bottomleft", legend=levels(class_MUAC_2), pch=15, col=c(2,1)) 

with(data,plot(pca$x,col=class_WFL,main="PCs by class WFL"))
legend("bottomleft", legend=levels(class_WFL), pch=15, col=c(4,3,2,1)) 
with(data,plot(pca$x,col=class_WFL_3,main="PCs by class WFL - 3 groups"))
legend("bottomleft", legend=levels(class_WFL_3), pch=15, col=c(3,2,1)) 
with(data,plot(pca$x,col=class_WFL_2,main="PCs by class WFL - 2 groups"))
legend("bottomleft", legend=levels(class_WFL_2), pch=15, col=c(2,1)) 

with(data,plot(pca$x,col=class_global,main="PCs by class global"))
legend("bottomleft", legend=levels(class_global), pch=15, col=c(4,3,2,1)) 
with(data,plot(pca$x,col=class_global_3,main="PCs by class global - 3 groups"))
legend("bottomleft", legend=levels(class_global_3), pch=15, col=c(3,2,1)) 
with(data,plot(pca$x,col=class_global_2,main="PCs by class global - 2 groups"))
legend("bottomleft", legend=levels(class_global_2), pch=15, col=c(2,1)) 

par(mfrow=c(1,1))
with(data,plot(pca$x,col=sex,main="PCs by sex"))
legend("topright", legend=levels(sex), pch=15, col=unique(sex),inset=c(-0.3,0))

with(data,plot(pca$x,col=age,main="PCs by age")) # posibles grupos
legend("topright", legend=levels(age), pch=15, col=unique(age),inset=c(-0.3,0))

with(data,plot(pca$x,col=zac_cat,main="PCs by ZAC")) # posibles grupos
legend("topright", legend=levels(zac_cat), pch=15, col=unique(age),inset=c(-0.3,0))

with(data,plot(pca$x,col=zbmi_cat,main="PCs by ZBMI")) # posibles grupos
legend("topright", legend=levels(zbmi_cat), pch=15, col=unique(age),inset=c(-0.3,0))

with(data,plot(pca$x,col=zlen_cat,main="PCs by ZLEN"))
legend("topright", legend=levels(zlen_cat), pch=15, col=unique(age),inset=c(-0.3,0))

with(data,plot(pca$x,col=zwei_cat,main="PCs by ZWEI"))
legend("topright", legend=levels(zwei_cat), pch=15, col=unique(age),inset=c(-0.3,0))

with(data,plot(pca$x,col=zwfl_cat,main="PCs by ZWFL")) # posibles grupos
legend("topright", legend=levels(zwfl_cat), pch=15, col=unique(age),inset=c(-0.3,0))

with(data,plot(pca$x,col=measure,main="PCs by measure")) # posibles grupos?
legend("topright", legend=levels(measure), pch=15, col=unique(age),inset=c(-0.3,0))



#### PC2 vs PC3 ####
par(mfrow=c(1,3))
with(data,plot(pca$x[,2:3],col=c(4,3,2,1),main="PCs by class MUAC"))
legend("bottomleft", legend=levels(class_MUAC), pch=15, col=c(4,3,2,1)) 
with(data,plot(pca$x[,2:3],col=c(3,2,1),main="PCs by class MUAC - 3 groups"))
legend("bottomleft", legend=levels(class_MUAC_3), pch=15, col=c(3,2,1)) 
with(data,plot(pca$x[,2:3],col=c(2,1),main="PCs by class MUAC - 2 groups"))
legend("bottomleft", legend=levels(class_MUAC_2), pch=15, col=c(2,1)) 

with(data,plot(pca$x[,2:3],col=c(4,3,2,1),main="PCs by class WFL"))
legend("bottomleft", legend=levels(class_WFL), pch=15, col=c(4,3,2,1)) 
with(data,plot(pca$x[,2:3],col=c(3,2,1),main="PCs by class WFL - 3 groups"))
legend("bottomleft", legend=levels(class_WFL_3), pch=15, col=c(3,2,1)) 
with(data,plot(pca$x[,2:3],col=c(2,1),main="PCs by class WFL - 2 groups"))
legend("bottomleft", legend=levels(class_WFL_2), pch=15, col=c(2,1)) 

with(data,plot(pca$x[,2:3],col=c(4,3,2,1),main="PCs by class global"))
legend("bottomleft", legend=levels(class_global), pch=15, col=c(4,3,2,1)) 
with(data,plot(pca$x[,2:3],col=c(3,2,1),main="PCs by class global - 3 groups"))
legend("bottomleft", legend=levels(class_global_3), pch=15, col=c(3,2,1)) 
with(data,plot(pca$x[,2:3],col=c(2,1),main="PCs by class global - 2 groups"))
legend("bottomleft", legend=levels(class_global_2), pch=15, col=c(2,1)) 

par(mfrow=c(1,1))
with(data,plot(pca$x[,2:3],col=sex,main="PCs by sex"))
legend("topright", legend=levels(sex), pch=15, col=unique(sex),inset=c(-0.3,0))

with(data,plot(pca$x[,2:3],col=age,main="PCs by age")) # posibles grupos
legend("topright", legend=levels(age), pch=15, col=unique(age),inset=c(-0.3,0))

with(data,plot(pca$x[,2:3],col=zac_cat,main="PCs by ZAC")) # posibles grupos
legend("topright", legend=levels(zac_cat), pch=15, col=unique(zac_cat),inset=c(-0.3,0))

with(data,plot(pca$x[,2:3],col=zbmi_cat,main="PCs by ZBMI")) # posibles grupos
legend("topright", legend=levels(zbmi_cat), pch=15, col=unique(zbmi_cat),inset=c(-0.3,0))

with(data,plot(pca$x[,2:3],col=zlen_cat,main="PCs by ZLEN"))
legend("topright", legend=levels(zlen_cat), pch=15, col=unique(zlen_cat),inset=c(-0.3,0))

with(data,plot(pca$x[,2:3],col=zwei_cat,main="PCs by ZWEI"))
legend("topright", legend=levels(zwei_cat), pch=15, col=unique(zwei_cat),inset=c(-0.3,0))

with(data,plot(pca$x[,2:3],col=zwfl_cat,main="PCs by ZWFL")) # posibles grupos
legend("topright", legend=levels(zwfl_cat), pch=15, col=unique(zwfl_cat),inset=c(-0.3,0))

with(data,plot(pca$x[,2:3],col=measure,main="PCs by measure")) # posibles grupos?
legend("topright", legend=levels(measure), pch=15, col=unique(measure),inset=c(-0.3,0))


### 3D representation ####
pca3d(pca, group=class_global_3)



### MULTIDIMENSIONAL SCALING ####

mds = cmdscale(mat,eig=TRUE,k=2)
mds

x = mds$points[,1]
y = mds$points[,2]

plot(x, y, xlab="Dim.1", ylab="Dim.2", xlim=c(-4,2),ylim=c(-2,3),col=measure) # 2 clusters?
text(x, y, cex=.5)

mds.df = as.data.frame(mds$points) 
kmclusters = kmeans(mds.df, centers = 2, nstart = 1) # k-means clustering with 2 groups 
kmclusters = as.factor(kmclusters$cluster) # convert to a factor 
mds.df$groups = kmclusters # join to the existing data frame mds.df
mds.df=cbind(data,mds.df) 

ggscatter(mds.df, 
          x = "V1", 
          y = "V2", 
          color = "measure", # sex, age, ZLEN_cat...
          palette = "jco", 
          size = 1, xlab = "Principal component 1", ylab = "Principal component 2",
          # ellipse = TRUE,
          #font.label = c(20, "plain"),
          # ellipse.type = "convex",
          repel = TRUE)

ggscatter(mds.df, 
          x = "V1", 
          y = "V2", 
          color = "groups", 
          palette = "jco", 
          size = 1, xlab = "Principal component 1", ylab = "Principal component 2",
          font.label = c(20, "plain"),          ellipse = TRUE,
          ellipse.type = "convex",
          repel = TRUE)

# cluster 1:
ggscatter(mds.df[mds.df$groups==1,], 
          x = "V1", 
          y = "V2", 
          color = "measure", 
          palette = "jco", 
          size = 1, xlab = "Principal component 1", ylab = "Principal component 2",
          # ellipse = TRUE, 
          # ellipse.type = "convex", 
          repel = TRUE)

# cluster 2:
ggscatter(mds.df[mds.df$groups==2,], 
          x = "V1", 
          y = "V2", 
          color = "measure", 
          palette = "jco", 
          size = 1, xlab = "Principal component 1", ylab = "Principal component 2",
          # ellipse = TRUE, 
          # ellipse.type = "convex", 
          repel = TRUE)

# Visualizing a correlation matrix using Multidimensional Scaling:
res.cor = cor(Proc_coords_R, method = "spearman")
mds.cor = (1 - res.cor) %>%
  cmdscale() %>%
  as_tibble()
colnames(mds.cor) = c("Dim.1", "Dim.2")
ggscatter(mds.cor, x = "Dim.1", y = "Dim.2",
          size = 1,
          label = colnames(res.cor),
          repel = TRUE)



### CANONICAL DISCRIMINANT ANALYSIS ####

### FUNCTIONS ####
# Hotelling's test - Morphometrics with R
Hotellingsp=function(SSef, SSer, dfef, dfer, exact=F){ # podemos quitar el exact y el if de la funcion creo!
  p = qr(SSef+SSer)$rank
  k=dfef; w=dfer
  s=min(k,p)
  m=(w-p-1)/2
  t1=(abs(p-k)-1)/2
  Ht=sum(diag(SSef%*%ginv(SSer)))
  Fapprox=Ht*(2 * (s*m+1))/(s^2*(2*t1+s+1))
  ddfnum=s*(2*t1+s+1)
  ddfden=2*(s*m+1)
  pval= 1-pf(Fapprox, ddfnum, ddfden)
  if (exact){ # este if se puede quitar porque exact es por default FALSE
    b=(p+2*m)*(k+2*m)/((2*m+1)*(2*m-2))
    c1=(2+(p*k+2)/(b-1))/(2*m)
    Fapprox=((4+(p*k+2)/(b-1))/(p*k))*(Ht/c1)
    ddfnum=p*k
    ddfden=4+(p*k+2)/(b-1)}
  unlist(list("dfeffect"=dfef,"dferror"=dfer,"T2"=Ht,"Approx_F"=Fapprox,
              "df1"=ddfnum,"df2"=ddfden,"p"=pval))
}

radios = function(g,p,n,conf.level=0.95) {
  N = sum(n)
  F = qf(conf.level,p,N-g-p+1)
  sqrt(F*(N-g)*p/((N-g-p+1)*n))
}

### Analysis ####

# visualizing the high correlations:
pairs.panels(Proc_coords_R[,1:5])
pairs.panels(Proc_coords_R[,6:11])
pairs.panels(Proc_coords_R[,12:17])
pairs.panels(Proc_coords_R[,18:23])
pairs.panels(Proc_coords_R[,24:29])
pairs.panels(Proc_coords_R[,29:34])
pairs.panels(Proc_coords_R[,35:40])

# removing high correlations
df2 = cor(Proc_coords_R)
hc = findCorrelation(df2, cutoff=0.9)
hc = sort(hc)
reducedData=Proc_coords_R[,-c(hc)]
reducedData=Proc_coords_R[,c(1:4,6:8,11,13,15,22,23,25:27,31,34,36,38:40)]


### class global ####

mod=lm(as.matrix(cbind(data[,c(44:83)])) ~ class_global, data=data)
dfef= length(levels(class_global))-1
n=dim(data)[1]
dfer= n - length(levels(class_global))
SSef=(n-1)*var(mod$fitted.values)
SSer=(n-1)*var(mod$residuals)

Hotellingsp(SSef, SSer, dfef, dfer)

# linear model
mod = lm(as.matrix(reducedData) ~ class_global)
Anova(mod, test="Wilks") # Manova

# canonical discriminant analysis
can1 = candisc(mod, term="class_global")
heplot(can1,error.ellipse= FALSE)
plot(can1,conf=0.9,type="n",main="Class global")

can2=candisc(mod,term="class_global",ndim=1)
plot(can2)

scores = as.matrix(reducedData) %*% can1$coeffs.raw
plot(scores[,1],scores[,2], xlab="Canonical axis 1",ylab="Canonical axis 2")
medias=aggregate(reducedData,list(data$class_global),mean)
scores.medias = as.matrix(medias[,-1]) %*% as.matrix(can1$coeffs.raw)
text(scores.medias[,1],scores.medias[,2],1:5,pch=15,col="blue")

resumen = table(data$class_global)
g = length(resumen) # number of populations
p = dim(reducedData)[2] # number of variables
n = as.vector(resumen) # sample size in each population
r = radios(g,p,n,0.90)

# plots
par(mfrow=c(1,2))
plot(can1,conf=0.9,type="n",main="Class global")
plot.new()
plot.window(xlim=c(-32,-29),ylim=c(118,121))
axis(1)
axis(2)
box()
title(main="Confidence intervals", xlab="Can1", ylab="Can2")
text(scores.medias[,1],scores.medias[,2],labels=levels(data$class_global),col="blue")
symbols(scores.medias[,1],scores.medias[,2],circles=r,inches=FALSE,add=T,lwd=2,fg="blue")


### class_WFL #### 

mod=lm(as.matrix(cbind(data[,c(44:83)])) ~ class_WFL, data=data)
dfef= length(levels(class_WFL))-1
n=dim(data)[1]
dfer= n - length(levels(class_WFL))
SSef=(n-1)*var(mod$fitted.values)
SSer=(n-1)*var(mod$residuals)

Hotellingsp(SSef, SSer, dfef, dfer)

# linear model
mod = lm(as.matrix(reducedData) ~ class_WFL)
Anova(mod, test="Wilks") # Manova

# canonical discriminant analysis
can1 = candisc(mod, term="class_WFL")
heplot(can1,error.ellipse= FALSE)
plot(can1,conf=0.9,type="n",main="Class global")

can2=candisc(mod,term="class_WFL",ndim=1)
plot(can2)

scores = as.matrix(reducedData) %*% can1$coeffs.raw
plot(scores[,1],scores[,2], xlab="Canonical axis 1",ylab="Canonical axis 2")
medias=aggregate(reducedData,list(data$class_WFL),mean)
scores.medias = as.matrix(medias[,-1]) %*% as.matrix(can1$coeffs.raw)
text(scores.medias[,1],scores.medias[,2],1:5,pch=15,col="blue")

resumen = table(data$class_WFL)
g = length(resumen) # number of populations
p = dim(reducedData)[2] # number of variables
n = as.vector(resumen) # sample size in each population
r = radios(g,p,n,0.90)

# plots
par(mfrow=c(1,2))
plot(can1,conf=0.9,type="n",main="Class global")
plot.new()
plot.window(xlim=c(-32,-29),ylim=c(118,121))
axis(1)
axis(2)
box()
title(main="Confidence intervals", xlab="Can1", ylab="Can2")
text(scores.medias[,1],scores.medias[,2],labels=levels(data$class_WFL),col="blue")
symbols(scores.medias[,1],scores.medias[,2],circles=r,inches=FALSE,add=T,lwd=2,fg="blue")


### class_MUAC  ####

mod=lm(as.matrix(cbind(data[,c(44:83)])) ~ class_MUAC, data=data)
dfef= length(levels(class_MUAC))-1
n=dim(data)[1]
dfer= n - length(levels(class_MUAC))
SSef=(n-1)*var(mod$fitted.values)
SSer=(n-1)*var(mod$residuals)

Hotellingsp(SSef, SSer, dfef, dfer)

# linear model
mod = lm(as.matrix(reducedData) ~ class_MUAC)
Anova(mod, test="Wilks") # Manova

# canonical discriminant analysis
can1 = candisc(mod, term="class_MUAC")
heplot(can1,error.ellipse= FALSE)
plot(can1,conf=0.9,type="n",main="Class global")

can2=candisc(mod,term="class_MUAC",ndim=1)
plot(can2)

scores = as.matrix(reducedData) %*% can1$coeffs.raw
plot(scores[,1],scores[,2], xlab="Canonical axis 1",ylab="Canonical axis 2")
medias=aggregate(reducedData,list(data$class_MUAC),mean)
scores.medias = as.matrix(medias[,-1]) %*% as.matrix(can1$coeffs.raw)
text(scores.medias[,1],scores.medias[,2],1:5,pch=15,col="blue")

resumen = table(data$class_MUAC)
g = length(resumen) # number of populations
p = dim(reducedData)[2] # number of variables
n = as.vector(resumen) # sample size in each population
r = radios(g,p,n,0.90)

# plots
par(mfrow=c(1,2))
plot(can1,conf=0.9,type="n",main="Class global")
plot.new()
plot.window(xlim=c(-32,-29),ylim=c(118,121))
axis(1)
axis(2)
box()
title(main="Confidence intervals", xlab="Can1", ylab="Can2")
text(scores.medias[,1],scores.medias[,2],labels=levels(data$class_MUAC),col="blue")
symbols(scores.medias[,1],scores.medias[,2],circles=r,inches=FALSE,add=T,lwd=2,fg="blue")


### AGE  ####

mod=lm(as.matrix(cbind(data[,c(44:83)])) ~ age, data=data)
dfef= length(levels(age))-1
n=dim(data)[1]
dfer= n - length(levels(age))
SSef=(n-1)*var(mod$fitted.values)
SSer=(n-1)*var(mod$residuals)

Hotellingsp(SSef, SSer, dfef, dfer)

# linear model
mod = lm(as.matrix(reducedData) ~ age)
Anova(mod, test="Wilks") # Manova

# canonical discriminant analysis
can1 = candisc(mod, term="age")
heplot(can1,error.ellipse= FALSE)
plot(can1,conf=0.9,type="n",main="Class global")

can2=candisc(mod,term="age",ndim=1)
plot(can2)

scores = as.matrix(reducedData) %*% can1$coeffs.raw
plot(scores[,1],scores[,2], xlab="Canonical axis 1",ylab="Canonical axis 2")
medias=aggregate(reducedData,list(data$age),mean)
scores.medias = as.matrix(medias[,-1]) %*% as.matrix(can1$coeffs.raw)
text(scores.medias[,1],scores.medias[,2],1:5,pch=15,col="blue")

resumen = table(data$age)
g = length(resumen) # number of populations
p = dim(reducedData)[2] # number of variables
n = as.vector(resumen) # sample size in each population
r = radios(g,p,n,0.90)

# plots
par(mfrow=c(1,2))
plot(can1,conf=0.9,type="n",main="Class global")
plot.new()
plot.window(xlim=c(-32,-29),ylim=c(118,121))
axis(1)
axis(2)
box()
title(main="Confidence intervals", xlab="Can1", ylab="Can2")
text(scores.medias[,1],scores.medias[,2],labels=levels(data$age),col="blue")
symbols(scores.medias[,1],scores.medias[,2],circles=r,inches=FALSE,add=T,lwd=2,fg="blue")


### zwfl_cat ####

mod=lm(as.matrix(cbind(data[,c(44:83)])) ~ zwfl_cat, data=data)
dfef= length(levels(zwfl_cat))-1
n=dim(data)[1]
dfer= n - length(levels(zwfl_cat))
SSef=(n-1)*var(mod$fitted.values)
SSer=(n-1)*var(mod$residuals)

Hotellingsp(SSef, SSer, dfef, dfer)

# linear model
mod = lm(as.matrix(reducedData) ~ zwfl_cat)
Anova(mod, test="Wilks") # Manova

# canonical discriminant analysis
can1 = candisc(mod, term="zwfl_cat")
heplot(can1,error.ellipse= FALSE)
plot(can1,conf=0.9,type="n",main="Class global")

can2=candisc(mod,term="zwfl_cat",ndim=1)
plot(can2)

scores = as.matrix(reducedData) %*% can1$coeffs.raw
plot(scores[,1],scores[,2], xlab="Canonical axis 1",ylab="Canonical axis 2")
medias=aggregate(reducedData,list(data$zwfl_cat),mean)
scores.medias = as.matrix(medias[,-1]) %*% as.matrix(can1$coeffs.raw)
text(scores.medias[,1],scores.medias[,2],1:5,pch=15,col="blue")

resumen = table(data$zwfl_cat)
g = length(resumen) # number of populations
p = dim(reducedData)[2] # number of variables
n = as.vector(resumen) # sample size in each population
r = radios(g,p,n,0.90)

# plots
par(mfrow=c(1,2))
plot(can1,conf=0.9,type="n",main="Class global")
plot.new()
plot.window(xlim=c(-32,-29),ylim=c(118,121))
axis(1)
axis(2)
box()
title(main="Confidence intervals", xlab="Can1", ylab="Can2")
text(scores.medias[,1],scores.medias[,2],labels=levels(data$zwfl_cat),col="blue")
symbols(scores.medias[,1],scores.medias[,2],circles=r,inches=FALSE,add=T,lwd=2,fg="blue")


### zlen_cat ####

mod=lm(as.matrix(cbind(data[,c(44:83)])) ~ zlen_cat, data=data)
dfef= length(levels(zlen_cat))-1
n=dim(data)[1]
dfer= n - length(levels(zlen_cat))
SSef=(n-1)*var(mod$fitted.values)
SSer=(n-1)*var(mod$residuals)

Hotellingsp(SSef, SSer, dfef, dfer)

# linear model
mod = lm(as.matrix(reducedData) ~ zlen_cat)
Anova(mod, test="Wilks") # Manova

# canonical discriminant analysis
can1 = candisc(mod, term="zlen_cat")
heplot(can1,error.ellipse= FALSE)
plot(can1,conf=0.9,type="n",main="Class global")

can2=candisc(mod,term="zlen_cat",ndim=1)
plot(can2)

scores = as.matrix(reducedData) %*% can1$coeffs.raw
plot(scores[,1],scores[,2], xlab="Canonical axis 1",ylab="Canonical axis 2")
medias=aggregate(reducedData,list(data$zlen_cat),mean)
scores.medias = as.matrix(medias[,-1]) %*% as.matrix(can1$coeffs.raw)
text(scores.medias[,1],scores.medias[,2],1:5,pch=15,col="blue")

resumen = table(data$zlen_cat)
g = length(resumen) # number of populations
p = dim(reducedData)[2] # number of variables
n = as.vector(resumen) # sample size in each population
r = radios(g,p,n,0.90)

# plots
par(mfrow=c(1,2))
plot(can1,conf=0.9,type="n",main="Class global")
plot.new()
plot.window(xlim=c(-32,-29),ylim=c(118,121))
axis(1)
axis(2)
box()
title(main="Confidence intervals", xlab="Can1", ylab="Can2")
text(scores.medias[,1],scores.medias[,2],labels=levels(data$zlen_cat),col="blue")
symbols(scores.medias[,1],scores.medias[,2],circles=r,inches=FALSE,add=T,lwd=2,fg="blue")


### SEX  ####
mod=lm(as.matrix(cbind(data[,c(44:83)])) ~ sex, data=data)
dfef= length(levels(sex))-1
n=dim(data)[1]
dfer= n - length(levels(sex))
SSef=(n-1)*var(mod$fitted.values)
SSer=(n-1)*var(mod$residuals)

Hotellingsp(SSef, SSer, dfef, dfer)

model=lm(as.numeric(sex) ~ .,data=Proc_coords_R)
stepAIC(model)
mod = lm(as.matrix(reducedData) ~ sex)
Anova(mod, test="Wilks") # Manova

can1 = candisc(mod, term="sex")
can1 
summary(can1)
plot(can1,conf=0.9,type="n")

# Box's M-test for Homogeneity of Covariance Matrices:
box_m(reducedData,data$sex) 
box_m(reducedData,data$class_MUAC) 
box_m(reducedData,data$class_WFL) 
box_m(reducedData,data$class_global) 
box_m(reducedData,data$zwei_cat) 
box_m(reducedData,data$zlen_cat) 
box_m(reducedData,data$zwfl_cat) 
box_m(reducedData,data$zbmi_cat) 
box_m(reducedData,data$zac_cat) 
# all p-vaLues are < 0.05 -> cov matrices are different