tifpath <- ""
setwd(tifpath)
library(tiff) 
library(raster)
library(rgdal)
library(MODIS)
library(gtools)

names<-list.files()
namesN<-names
data1 <- list()
for (i in 1:6) {
  data1[[namesN[i]]] <- readTIFF(namesN[i])
  data1[[namesN[i]]]<-data1[[namesN[i]]][,,1]
  data1[[namesN[i]]]<-data1[[namesN[i]]][1:496, 497:1520]
}

datacol<-ncol(data1[[namesN[1]]])
datarow<-nrow(data1[[namesN[1]]])
Guess<-matrix(c(0.075,0.124,0.075,0.124,0.204,0.124,0.075,0.124,0.075),nrow=3,ncol=3)
dataN<-list()
for(i in 1:6){
  dataN[[namesN[i]]]<-data1[[namesN[i]]]
}

for(i in 1:6){
  for(t in 2:(datarow-1)){
    for(k in 2:(datacol-1)){
      Block<-data1[[namesN[i]]][(t-1):(t+1),(k-1):(k+1)]
      B<-Block*Guess
      BSum<-sum(B)
      dataN[[namesN[i]]][t,k]<-BSum
    }
  }
}
Dis<-function(a,b){
  D<-sqrt((a[1]-b[1])^2+(a[2]-b[2])^2)
  return(D)
}
dataO<-dataN
for(i in 1:6){
  totsu<-otsu(dataN[[namesN[i]]])
  dataO[[namesN[i]]][dataO[[namesN[i]]]>totsu]<-1
  dataO[[namesN[i]]][dataO[[namesN[i]]]<totsu]<-0
}
AveMean<-matrix(c(1/9,1/9,1/9,1/9,1/9,1/9,1/9,1/9,1/9),nrow=3,ncol=3)
for(i in 1:6){
  for(t in 2:(datarow-1)){
    for(k in 2:(datacol-1)){
      Block<-dataO[[namesN[i]]][(t-1):(t+1),(k-1):(k+1)]
      B<-Block*AveMean
      BSum<-sum(B)
      dataO[[namesN[i]]][t,k]<-BSum
    }
  }
}
for(i in 1:6){
  totsu2<-5/9
  dataO[[namesN[i]]][dataO[[namesN[i]]]>totsu2]<-1
  dataO[[namesN[i]]][dataO[[namesN[i]]]<totsu2]<-0
}
for(i in 1:6){
  for(t in 2:(datarow-1)){
    for(k in 2:(datacol-1)){
      Block<-dataO[[namesN[i]]][(t-1):(t+1),(k-1):(k+1)]
      B<-Block*AveMean
      BSum<-sum(B)
      dataO[[namesN[i]]][t,k]<-BSum
    }
  }
}
for(i in 1:6){
  totsu2<-5/9
  dataO[[namesN[i]]][dataO[[namesN[i]]]>totsu2]<-1
  dataO[[namesN[i]]][dataO[[namesN[i]]]<totsu2]<-0
}


a = list()
a[[1]]=c(283, 249)
a[[2]]=c(631, 218)
a[[3]]=c(335, 254)
a[[4]]=c(642, 214)
a[[5]]=c(383, 249)
a[[6]]=c(645, 223)
a[[7]]=c(390, 256)
a[[8]]=c(659, 228)
a[[9]]=c(407, 245)
a[[10]]=c(688, 233)
a[[11]]=c(399, 231)
a[[12]]=c(719, 242)

for(i in 1:12){
  a[[i]][1]<-a[[i]][1]-512
}
a11=a
for(i in 1:12){
  a[[i]][1]<-a11[[i]][2]
  a[[i]][2]<-a11[[i]][1]
}

Edge<-list()
for(i in 1:6){
  Line<-list()
  for(k in 1:datacol){
    Line[[k]]<-c(which(dataO[[namesN[i]]][,k]>0)[1],k-512)
    Line[[k]][is.na(Line[[k]])]<-0
  } 
  Edge[[i]]<-Line
}

CupV=0
RimV=0
for(i in 1:6){
  for(k in a[[2*i-1]][2]:a[[2*i]][2]){
    kline<-(a[[2*i]][1]-a[[2*i-1]][1])/(a[[2*i]][2]-a[[2*i-1]][2])
    yline<-kline*(k-a[[2*i-1]][2])+a[[2*i-1]][1]
    if((yline-Edge[[i]][[k+512]][1])>0){
      RimV=RimV+(yline-Edge[[i]][[k+512]][1])/50*1/36*2*pi*abs(k)/36/6
    }
    else{
      CupV=CupV+(Edge[[i]][[k+512]][1]-yline)/50*1/36*2*pi*abs(k)/36/6   
    }
  }
}
CupV
RimV