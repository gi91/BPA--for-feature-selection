#########################################################
#Install package gRapHD which it is important for the BPA
#########################################################
#install.packages("remotes")
#install.packages("BiocManager")
#BiocManager::install("graph")
#remotes::install_version("gRapHD", version = "0.2.6")



library(igraph)
library(gRapHD)
library(gRbase)
library(gdata)
library(gRim)
library(gRbase)
library(gRain)
library(ggplot2)
library(mlbench)
library(GGally)
library(rlang)
library(tidyverse)
library("philentropy")
library(NMF)
library(np)
library(fastmit)




#########################################################
#Creation Dataset 
#########################################################
n <- 200
p<-40
sigma <- 1
set.seed(1235469744)
sim <- mlbench.friedman1(n, sd = sigma)
colnames(sim$x) <- c(paste("X", 1:5, sep = ""),paste("Noise", 1:5, sep = ""))
bogus <- matrix(rnorm(n * p), nrow = n)
colnames(bogus) <- paste("Noise", 5+(1:ncol(bogus)), sep = "")
x <- cbind(sim$x, bogus)
Y<- sim$y
x<-cbind(Y,x)


#########################################################
#Creation of he Dataset of interest
#########################################################
x<-as.data.frame(x)


#########################################################
#Step 0-creation of the HDGM
#########################################################
m <- minForest(x,stat = "AIC")
#Here you can USE LR, or AIC or BIC penalization: in this example we use AIC in order
# to introduce more noise to see how perform the Best Path Algorthim
print(m)
dev.off()
plot(m, numIter=2000 ,main="Min AIC forest",vert.radii = 0.05,cex.vert.label = 0.7)


k<-as.array(colnames(x))
#Adjacent MAtrix from HDGM

AM<-adjMat(m)
colnames(AM)<-c(k)
row.names(AM)<-c(k)
AM<-as.data.frame(AM)
nam_verx<-c(k)
numb_verx<-c(1:ncol(x))
nam_verx<-as.data.frame(nam_verx)
nam_verx<-cbind(nam_verx,numb_verx)
table(Degree(m))
vertices <- which(Degree(m) >= 1); vertices
str(x)


##################################################################
# Step 1: selection node of interest
##################################################################

neigh<- neighbourhood(m, orig = vertices[1],rad = 1001)
#Here we put vertices[1] because our variable of inerest is the node 1 
hlv<- c(1,neigh$v[-1, 1])
hlc <- rep(c( "green","orange"), c(1, length(neigh$v[-1, 1])))
dev.off()
plot(m, vert = neigh$v[,1], numIter=2000,vert.hl = hlv, vert.radii = 0.05,vert.label=TRUE ,col.hl = hlc, cex.vert.label = 0.7)
#We cheek with this plot if our variable selection strs from the node of interest,- The green one



#This procedure needs only for prety plot
AM<-adjMat(m)
colnames(AM)<-c(1:ncol(x))
row.names(AM)<-c(1:ncol(x))
AM<-as.data.frame(AM)

upperTriangle(AM, diag=FALSE,byrow = TRUE)<-0

AM[AM == 0] <- NA
AM<-as.data.frame(AM)
AM$from<-c(1:ncol(x))

dim(AM)



#transform 
edges<- AM%>% 
  gather(key="to", value="value", -(ncol(x)+1)) %>%
  na.omit()

dim(edges)
rownames(edges)<-c(1:dim(edges)[1])

nodes<-c(1:ncol(x))
nodes<-as.data.frame(nodes)
colnames(nodes)<-"id"
label<-c(colnames(x))
nodes<-cbind(nodes,label)
nodes<-as.data.frame(nodes)


nodes$color.background<-c("gold",rep("chartreuse1",time=50))

edges$value<-NULL
edges$color<-"black"



net<- graph.data.frame(edges, nodes, directed=F)
ggnet2(net)
colnames(x)


a<-ggnet2(net, node.size = 13, node.color =c(nodes$color.background), 
          alpha = 3,label = colnames(x), label.size = 3,label.color ="black",fontface = "bold")


a




#Dataset of interest if we are in the forest case
sel<-neigh$v[,1]
X<-x[,sel]
dim(X)
X<-data.frame(X)

#dataset of support for the different sub ses
p<-as.data.frame(neigh$v)
colnames(p)<-c("numb_verx","dist_1")
l<- merge(p,nam_verx,by="numb_verx")
l<-as.data.frame(l);dim(l)
rownames(l)<-c(1:dim(l)[1])
l<-l[order(l$dist_1),]






l<-l[-1,]
dim(l)
l$Order<-c(1:dim(l)[1])
dataf = matrix(0, nrow = nrow(x), ncol = dim(l)[1]);dim(dataf)
dataf<-as.data.frame(dataf)
l

P<-X[,-1]
colnames(P)
y<-X$Y
str(P)
j<-ncol(P)
##################################################################
# Step 2 :Estimation of f(Y|X)
##################################################################
#Algoritmo
for(i in 1:j){
  T<-P[,l$Order==i]
  set.seed(123456)
  # Train the model
  bw <-npcdensbw(xdat=T, ydat=y) 
  fhat <-fitted(npcdens(bws=bw))
  dataf[,i]<-c(fhat)
}

colnames(dataf)<-c(colnames(P))


##################################################################
# BEST PATH STEP
##################################################################
j<-max(l$dist)
datakL = matrix(NA, nrow = j, ncol = 1);dim(datakL)
datakL<-as.data.frame(datakL)
rownames(datakL)<-paste("path",1:j,sep="")
colnames(datakL)<-"EC"



#Algoritmo
for(i in 1:j){
  set.seed(123456)
  T<-dataf[,l$dist<=i]
  T<-as.data.frame(T)
  set.seed(123)
  # Train the model
  fy<-rowSums(T)/ncol(T)
  r<-ncol(T)
  KLM<-matrix(NA,nrow=r,ncol=1)
  for(o in 1:r){
    set.seed(123)
    h<-rbind(fy,T[,o])
    a1<-KL(h, est.prob = "empirical")
    h<-rbind(T[,o],fy)
    a2<-KL(h, est.prob = "empirical")
    a<-0.5*(a1+a2)
    KLM[o,]<-c(a)
  }
  EC<-mean(KLM[,1])
  datakL[i,]<-c(EC)
}
##################################################################
# Step 3:  pick as best-path-step the one with the highest EC
##################################################################
which.max(datakL$EC)
datakL
P<-P[,l$dist<=which.max(datakL$EC)]

##################################################################
# Step 4: implement the Kraskov test
##################################################################
datp= matrix(NA, nrow = ncol(P), ncol = 1);dim(datp)
datp<-as.data.frame(datp)
rownames(datp)<-paste(colnames(P))
colnames(datp)<-"p"


for(i in 1:ncol(P)){
  a<-mi.test(P[,i],y,distance = FALSE, num.permutations = 99,
             seed = 1233)
  datp[i,]<-c(a$p.value)
  
}


datp
##################################################################
# Linear Case
##################################################################
#Dataset of interest 

sel<-neigh$v[,1]
X<-x[,sel]
dim(X)
X<-data.frame(X)



#dataset of support for the different sub ses
p<-as.data.frame(neigh$v)
colnames(p)<-c("numb_verx","dist")
l<- merge(p,nam_verx,by="numb_verx")
l<-as.data.frame(l);dim(l)
l<-l[order(l$dist),]
dim(l)
rownames(l)<-c(1:ncol(X))

max(l$dist)
j<-max(l$dist)


#train.control <- trainControl(method = "cv", number = 15)
dataf = matrix(, nrow = j, ncol = 3);dim(dataf)
dataf<-as.data.frame(dataf)
colnames(dataf)<-c("R2","ADjR2","Rcoef")
rownames(dataf)<-paste("path",1:j,sep="")
#Algoritmo
for(i in 1:j){
  T<-X[,l$dist<=i]
  set.seed(123)
  # Train the model
  k<- lm( Y~., data = T)
  d<-summary(k)
  Rcoef<-(d$adj.r.squared/(1-d$adj.r.squared))
  dataf[i,]<-c(d$r.squared,d$adj.r.squared,Rcoef)
}
dataf


X<-X[,l$dist<=which.max(dataf$Rcoef)]
colnames(X)
OLS1<-lm(Y~.,data =X)
summary(OLS1)
colnames(X)
xg<-X[,c(1:5)]

OLS<-lm(Y~.,data =xg)
summary(OLS)


#Please if you use this code cite: Feature selection based on the best-path algorithm in high dimensional graphical models
#@article{riso2023feature,
#title={Feature selection based on the best-path algorithm in high dimensional graphical models},
#author={Riso, Luigi and Zoia, Maria G and Nava, Consuelo R},
#journal={Information Sciences},
#volume={649},
#pages={119601},
#year={2023},
#publisher={Elsevier}
#}