#load package
library(ncdf4)
library(chron)
library(lattice)
library(RColorBrewer)
library(plotrix)
library(forecast)
library(FNN)
library(gaussDiff)
library(GoFKernel)
library(geoR)
library(ggplot2)

setwd("")
cord.hist<-readRDS("CORDEX4_historical.RDA")
cord.areg<-readRDS("cord.areg.RDS")
merra.hist<-readRDS("merra_hist.RDS")

year<-c(1980:2005)
years<-rep(year,times=c(365,365,365,365,365,365,365,365,365,365,365,365,365,365,365,365,365,365,365,365,365,365,365,365,365,365))
luna<-seq(60,9497,by=1460)

rm.m<-array(0,c(35,34,9490))
fit.m<-array(0,c(35,34,9490))
for (i in 1:34){
  for (j in 1:35){
    ob<-ts(merra.hist[j,i,-luna],start=c(1980,1),end=c(2005,365),fr=365)
    X=fourier(ob,K=2)
    mod.hr=tslm(ob~X+years) 
    resi<-as.numeric(mod.hr$residuals)
    fit<-as.numeric(mod.hr$fitted.values)
    fit.m[j,i,]<-fit
    rm.m[j,i,]<-resi
  }
}

rm.c<-array(0,c(35,34,9490))
fit.c<-array(0,c(35,34,9490))
for (i in 1:34){
  for (j in 1:35){
    ob<-ts(cord.hist[j,i,-luna],start=c(1980,1),end=c(2005,365),fr=365)
    X=fourier(ob,K=2)
    mod.hr=tslm(ob~X+years) 
    resi<-as.numeric(mod.hr$residuals)
    rm.c[j,i,]<-resi
    fit<-as.numeric(mod.hr$fitted.values)
    fit.c[j,i,]<-fit
  }
}

rm.a<-array(0,c(35,34,9490))
fit.a<-array(0,c(35,34,9490))
for (i in 1:34){
  for (j in 1:35){
    ob<-ts(cord.areg[j,i,-luna],start=c(1980,1),end=c(2005,365),fr=365)
    X=fourier(ob,K=2)
    mod.hr=tslm(ob~X+years) 
    resi<-as.numeric(mod.hr$residuals)
    rm.a[j,i,]<-resi
    fit<-as.numeric(mod.hr$fitted.values)
    fit.a[j,i,]<-fit
  }
}

newrm.merra.hist<-array(0,c(35,34,9490))
newrm.cord.hist<-array(0,c(35,34,9490))
newrm.cord.areg<-array(0,c(35,34,9490))
for(i in 1:35){
  for (j in 1:34){
    ar<-arima(rm.m[i,j,],order=c(2,0,0))
    newrm.merra.hist[i,j,]<-ar$residuals
    ar<-arima(rm.c[i,j,],order=c(2,0,0))
    newrm.cord.hist[i,j,]<-ar$residuals
    ar<-arima(rm.a[i,j,],order=c(2,0,0))
    newrm.cord.areg[i,j,]<-ar$residuals
  }
}

saveRDS(newrm.cord.hist,file="newrm.cord.hist.RDS")
saveRDS(newrm.cord.areg,file="newrm.cord.areg.RDS")
saveRDS(newrm.merra.hist,file="newrm.merra.hist.RDS")