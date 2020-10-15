library(forecast)
library(FNN)
library(gaussDiff)
library(GoFKernel)
library(geoR)
library(ggplot2)

pred.m.fit<-readRDS("pred.RDS")
loc<-readRDS("/loc.RDS")
lac<-readRDS("lac.RDS")
turb<-readRDS("Optimal_DataFrame-Borders50km-16GW-run4.RDS")
turb<-turb[!is.na(turb[,3]),]
merra_mask<-readRDS("merra_mask.RDS")
cord.row<-readRDS("merra.hist.row.RDS")


mask.row<-c(rbind(merra_mask[1:35,]))

matern.cord<-cbind(loc,lac,cord.row)

vg <- variog(coords=matern.m[,1:2],data=matern.cord[,3:9492])
h  <- vg$u
vg$v <- rowMeans(vg$v)

cord.fit<-variofit(vg,kappa=0.5,fix.kappa = T)

krig.cord<-matrix(0,nrow=41071,ncol=20)
for(d in 1:20){
  day1<-cbind(loc,lac,cord.row[,d])
  d1<-as.geodata(day1,data.col = 3)   
  gr=cbind(turb$Lon,turb$Lat)
  kc= krige.control(type="ok",obj.mod = cord.fit)
  sk= krige.conv(d1,krige=kc,loc=gr)
  krig.merra[,d]<-sk$predict
}

matern.m<-cbind(loc,lac,pred.m.fit)

vg.m <- variog(coords=matern.m[,1:2],data=matern.m[,3:9492])
h.m  <- vg$u
vg.m$v <- rowMeans(vg$v)

m.fit<-variofit(vg,kappa=0.5,fix.kappa = T)


krig.merra<-matrix(0,nrow=41071,ncol=20)
for(d in 1:20){
  day1<-cbind(loc,lac,pred.m.fit[,d])
  d1<-as.geodata(day1,data.col = 3)   
  gr=cbind(turb$Lon,turb$Lat)
  kc= krige.control(type="ok",obj.mod = m.fit)
  sk= krige.conv(d1,krige=kc,loc=gr)
  krig[,d]<-sk$predict
}


turb.full<-readRDS("Optimal_DataFrame-Borders50km-16GW-run4.RDS")
turb$mw<-ifelse(turb$OptTurb==1|turb$OptTurb==2|turb$OptTurb==3|turb$OptTurb==4|turb$OptTurb==5|
                  turb$OptTurb==6,2,
                ifelse(turb$OptTurb==7|turb$OptTurb==8|turb$OptTurb==9|turb$OptTurb==10|turb$OptTurb==11|
                         turb$OptTurb==12,3.45,
                       ifelse(turb$OptTurb==13,2.75,
                              ifelse(turb$OptTurb==14|turb$OptTurb==15|turb$OptTurb==16|turb$OptTurb==17|turb$OptTurb==18,
                                     3.4,
                                     ifelse(turb$OptTurb==19|turb$OptTurb==20|turb$OptTurb==21|turb$OptTurb==22,2.75,
                                            ifelse(turb$OptTurb==23|turb$OptTurb==24|turb$OptTurb==25,2.5,
                                                   ifelse(turb$OptTurb==26|turb$OptTurb==27|turb$OptTurb==28|turb$OptTurb==29|turb$OptTurb==30|turb$OptTurb==31,3.3,3)))))))

turb$number<-floor(turb$PHI/turb$mw)


turb$height<-ifelse(turb$OptTurb==1,75,
                    ifelse(turb$OptTurb==2,80,
                           ifelse(turb$OptTurb==3,95,
                                  ifelse(turb$OptTurb==4,110,
                                         ifelse(turb$OptTurb==5,120,
                                                ifelse(turb$OptTurb==6,125,
                                                       ifelse(turb$OptTurb==7,87,
                                                              ifelse(turb$OptTurb==8,117,
                                                                     ifelse(turb$OptTurb==9,137,
                                                                            ifelse(turb$OptTurb==10,147,
                                                                                   ifelse(turb$OptTurb==11,149,
                                                                                          ifelse(turb$OptTurb==12,166,
                                                                                                 ifelse(turb$OptTurb==13,90,
                                                                                                        ifelse(turb$OptTurb==14,85,
                                                                                                               ifelse(turb$OptTurb==15,110,
                                                                                                                      ifelse(turb$OptTurb==16,131.4,
                                                                                                                             ifelse(turb$OptTurb==17,134,
                                                                                                                                    ifelse(turb$OptTurb==18,164.5,
                                                                                                                                           ifelse(turb$OptTurb==19,75,
                                                                                                                                                  ifelse(turb$OptTurb==20,85,
                                                                                                                                                         ifelse(turb$OptTurb==21,98.3,
                                                                                                                                                                ifelse(turb$OptTurb==22,123.5,
                                                                                                                                                                       ifelse(turb$OptTurb==23,75,
                                                                                                                                                                              ifelse(turb$OptTurb==24,80,
                                                                                                                                                                                     ifelse(turb$OptTurb==25,100,
                                                                                                                                                                                            ifelse(turb$OptTurb==26,84,
                                                                                                                                                                                                   ifelse(turb$OptTurb==27,106,
                                                                                                                                                                                                          ifelse(turb$OptTurb==28,112,
                                                                                                                                                                                                                 ifelse(turb$OptTurb==29,114,
                                                                                                                                                                                                                        ifelse(turb$OptTurb==30,120,
                                                                                                                                                                                                                               ifelse(turb$OptTurb==31,134,100)))))))))))))))))))))))))))))))


typ1<-c(1:41071)[turb$OptTurb==26|turb$OptTurb==27|turb$OptTurb==28|turb$OptTurb==29|turb$OptTurb==30|turb$OptTurb==31]

typ2<-c(1:41071)[turb$OptTurb==14|turb$OptTurb==15|turb$OptTurb==16|turb$OptTurb==17|turb$OptTurb==18]
typ3<-c(1:41071)[turb$OptTurb==19|turb$OptTurb==20|turb$OptTurb==21|turb$OptTurb==22]
typ4<-c(1:41071)[turb$OptTurb==23|turb$OptTurb==24|turb$OptTurb==25]
typ5<-c(1:41071)[turb$OptTurb==13]


alpha<-readRDS("alpha_day_100m_2016.RDS")
sig<-readRDS("sigma_day_100m_2016.RDS")
sd<-readRDS("stdev_day_100m_2016.RDS")
sa_lat<-readRDS("LAT_GRID.RDS")
sa_lon<-readRDS("LON_GRID.RDS")
sa_mask<-readRDS("saudi_mask-run4.RDS")

lat<-readRDS("lat.RDS")

lon<-readRDS("lon.RDS")

#alpha.row<-matrix(0,ncol=366,nrow = )
#for(i in 1:366){
  alpha.row[,i]<-cbind(c(alpha[,,i]))
#}

vara<-array(0,c(549,499,366))
for(d in 1:366){
  vara[,,d]<-sqrt(sd[,,d]^2)+(sig[,,d]^2)
}

day<-c(1:366,1:365,1:365,1:365,1:366,1:365,1:365,1:365,1:366,1:365,1:365,1:365,1:366,1:365,1:365,1:365,1:366,1:365,1:365,1:365,1:366,1:365,1:365,1:365,1:366,1:365)
diff<-matrix(0,nrow=41071,ncol=9497)
expfuture<-matrix(0,nrow=41071,ncol=9497)
for(i in 1:9497){
  alpha.row<-c(cbind(alpha[,,day[i]]))
  sd.row<-c(cbind(vara[,,day[i]]))
  al<-alpha.row[!is.na(turb.full$OptTurb)]
  sdp<-sd.row[!is.na(turb.full$OptTurb)]
  alp<-c()
  for(j in 1:41071){
    a<-rnorm(1,mean=al[j],sd=sdp[j])
    alp<-c(alp,a)
  }
  krig<-krig120[,i]
  expkrig<-krig*((turb$height/10)^alp)
  expfuture[,i]<-expkrig
  power<-c(1:41071)
  power[typ1]<-ifelse(expkrig[typ1]<3.5,0,
                      ifelse(expkrig[typ1]>=3.5&expkrig[typ1]<4,106,
                             ifelse(expkrig[typ1]>=4&expkrig[typ1]<4.5,197,
                                    ifelse(expkrig[typ1]>=4.5&expkrig[typ1]<5,311,
                                           ifelse(expkrig[typ1]>=5&expkrig[typ1]<5.5,447,
                                                  ifelse(expkrig[typ1]>=5.5&expkrig[typ1]<6,610,
                                                         ifelse(expkrig[typ1]>=6&expkrig[typ1]<6.5,804,
                                                                ifelse(expkrig[typ1]>=6.5&expkrig[typ1]<7,1032,
                                                                       ifelse(expkrig[typ1]>=7&expkrig[typ1]<7.5,1298,
                                                                              ifelse(expkrig[typ1]>=7.5&expkrig[typ1]<8,1601,
                                                                                     ifelse(expkrig[typ1]>=8&expkrig[typ1]<8.5,1936,
                                                                                            ifelse(expkrig[typ1]>=8.5&expkrig[typ1]<9,2292, 
                                                                                                   ifelse(expkrig[typ1]>=9&expkrig[typ1]<9.5,2635,
                                                                                                          ifelse(expkrig[typ1]>=9.5&expkrig[typ1]<10,2901,
                                                                                                                 ifelse(expkrig[typ1]>=10&expkrig[typ1]<10.5,3091,
                                                                                                                        ifelse(expkrig[typ1]>=10.5&expkrig[typ1]<11,3215,
                                                                                                                               ifelse(expkrig[typ1]>=11&expkrig[typ1]<11.5,3281,3300)))))))))))))))))
  
  
  power[typ2]<-ifelse(expkrig[typ2]<3.5,0,
                      ifelse(expkrig[typ2]>=3.5&expkrig[typ2]<4,140,
                             ifelse(expkrig[typ2]>=4&expkrig[typ2]<4.5,234,
                                    ifelse(expkrig[typ2]>=4.5&expkrig[typ2]<5,362,
                                           ifelse(expkrig[typ2]>=5&expkrig[typ2]<5.5,497,
                                                  ifelse(expkrig[typ2]>=5.5&expkrig[typ2]<6,686,
                                                         ifelse(expkrig[typ2]>=6&expkrig[typ2]<6.5,880,
                                                                ifelse(expkrig[typ2]>=6.5&expkrig[typ2]<7,1142,
                                                                       ifelse(expkrig[typ2]>=7&expkrig[typ2]<7.5,1407,
                                                                              ifelse(expkrig[typ2]>=7.5&expkrig[typ2]<8,1734,
                                                                                     ifelse(expkrig[typ2]>=8&expkrig[typ2]<8.5,2060,
                                                                                            ifelse(expkrig[typ2]>=8.5&expkrig[typ2]<9,2386,
                                                                                                   ifelse(expkrig[typ2]>=9&expkrig[typ2]<9.5,2709,
                                                                                                          ifelse(expkrig[typ2]>=9.5&expkrig[typ2]<10,2935,
                                                                                                                 ifelse(expkrig[typ2]>=10&expkrig[typ2]<10.5,3156,
                                                                                                                        ifelse(expkrig[typ2]>=10.5&expkrig[typ2]<11,3280,
                                                                                                                               ifelse(expkrig[typ2]>=11&expkrig[typ2]<11.5,3368,
                                                                                                                                      ifelse(expkrig[typ2]>=11.5&expkrig[typ2]<12,3415,
                                                                                                                                             ifelse(expkrig[typ2]>=12&expkrig[typ2]<12.5,3441,3500)))))))))))))))))))
  
  power[typ3]<-ifelse(expkrig[typ3]<3.5,0,
                      ifelse(expkrig[typ3]>=3.5&expkrig[typ3]<4,68,
                             ifelse(expkrig[typ3]>=4&expkrig[typ3]<4.5,114,
                                    ifelse(expkrig[typ3]>=4.5&expkrig[typ3]<5,177,
                                           ifelse(expkrig[typ3]>=5&expkrig[typ3]<5.5,243,
                                                  ifelse(expkrig[typ3]>=5.5&expkrig[typ3]<6,347,
                                                         ifelse(expkrig[typ3]>=6&expkrig[typ3]<6.5,452,
                                                                ifelse(expkrig[typ3]>=6.5&expkrig[typ3]<7,595,
                                                                       ifelse(expkrig[typ3]>=7&expkrig[typ3]<7.5,738,
                                                                              ifelse(expkrig[typ3]>=7.5&expkrig[typ3]<8,907,
                                                                                     ifelse(expkrig[typ3]>=8&expkrig[typ3]<8.5,1076,
                                                                                            ifelse(expkrig[typ3]>=8.5&expkrig[typ3]<9,1307,
                                                                                                   ifelse(expkrig[typ3]>=9.5&expkrig[typ3]<10,1786,
                                                                                                          ifelse(expkrig[typ3]>=10&expkrig[typ3]<10.5,2033,
                                                                                                                 ifelse(expkrig[typ3]>=10.5&expkrig[typ3]<11,2219,
                                                                                                                        ifelse(expkrig[typ3]>=11&expkrig[typ3]<11.5,2405,
                                                                                                                               ifelse(expkrig[typ3]>=11.5&expkrig[typ3]<12,2535,
                                                                                                                                      ifelse(expkrig[typ3]>=12&expkrig[typ3]<12.5,2633,
                                                                                                                                             ifelse(expkrig[typ3]>=12.5&expkrig[typ3]<13,2710,2750)))))))))))))))))))
  
  
  
  power[typ4]<-ifelse(expkrig[typ4]<3.5,0,
                      ifelse(expkrig[typ4]>=3.5&expkrig[typ4]<4,50,
                             ifelse(expkrig[typ4]>=4&expkrig[typ4]<4.5,136,
                                    ifelse(expkrig[typ4]>=4.5&expkrig[typ4]<5,221,
                                           ifelse(expkrig[typ4]>=5&expkrig[typ4]<5.5,326,
                                                  ifelse(expkrig[typ4]>=5.5&expkrig[typ4]<6,431,
                                                         ifelse(expkrig[typ4]>=6&expkrig[typ4]<6.5,576,
                                                                ifelse(expkrig[typ4]>=6.5&expkrig[typ4]<7,720,
                                                                       ifelse(expkrig[typ4]>=7&expkrig[typ4]<7.5,911,
                                                                              ifelse(expkrig[typ4]>=7.5&expkrig[typ4]<8,1102,
                                                                                     ifelse(expkrig[typ4]>=8&expkrig[typ4]<8.5,1339,
                                                                                            ifelse(expkrig[typ4]>=8.5&expkrig[typ4]<9,1575,
                                                                                                   ifelse(expkrig[typ4]>=9&expkrig[typ4]<9.5,1575,
                                                                                                          ifelse(expkrig[typ4]>=9.5&expkrig[typ4]<10,1797,
                                                                                                                 ifelse(expkrig[typ4]>=10&expkrig[typ4]<10.5,2019,
                                                                                                                        ifelse(expkrig[typ4]>=10.5&expkrig[typ4]<11,2162,
                                                                                                                               ifelse(expkrig[typ4]>=11&expkrig[typ4]<11.5,2304,
                                                                                                                                      ifelse(expkrig[typ4]>=11.5&expkrig[typ4]<12,2390,
                                                                                                                                             ifelse(expkrig[typ4]>=12&expkrig[typ4]<12.5,2458,2500)))))))))))))))))))
  power[typ5]<-ifelse(expkrig[typ5]<3,0,ifelse(expkrig[typ5]>=3&expkrig[typ5]<3.5,89,
                                               ifelse(expkrig[typ5]>=3.5&expkrig[typ5]<4,171,
                                                      ifelse(expkrig[typ5]>=4&expkrig[typ5]<4.5,269,
                                                             ifelse(expkrig[typ5]>=4.5&expkrig[typ5]<5,389,
                                                                    ifelse(expkrig[typ5]>=5&expkrig[typ5]<5.5,533,
                                                                           ifelse(expkrig[typ5]>=5.5&expkrig[typ5]<6,704,
                                                                                  ifelse(expkrig[typ5]>=6&expkrig[typ5]<6.5,906,
                                                                                         ifelse(expkrig[typ5]>=6.5&expkrig[typ5]<7,1136,
                                                                                                ifelse(expkrig[typ5]>=7&expkrig[typ5]<7.5,1400,
                                                                                                       ifelse(expkrig[typ5]>=7.5&expkrig[typ5]<8,1674,
                                                                                                              ifelse(expkrig[typ5]>=8&expkrig[typ5]<8.5,1945,
                                                                                                                     ifelse(expkrig[typ5]>=8.5&expkrig[typ5]<9,2173,
                                                                                                                            ifelse(expkrig[typ5]>=9&expkrig[typ5]<9.5,2373,
                                                                                                                                   ifelse(expkrig[typ5]>=9.5&expkrig[typ5]<10,2518,
                                                                                                                                          ifelse(expkrig[typ5]>=10&expkrig[typ5]<10.5,2619,
                                                                                                                                                 ifelse(expkrig[typ5]>=10.5&expkrig[typ5]<11,2696,
                                                                                                                                                        ifelse(expkrig[typ5]>=11&expkrig[typ5]<11.5,2739,
                                                                                                                                                               ifelse(expkrig[typ5]>=11.5&expkrig[typ5]<12,2766,2780)))))))))))))))))))
  
  histwind<-merra_krig120[,i]
  exphist<-histwind*((turb$height/10)^alp)
  
  powerhist<-c(1:41071)
  powerhist[typ1]<-ifelse(exphist[typ1]<3.5,0,
                          ifelse(exphist[typ1]>=3.5&exphist[typ1]<4,106,
                                 ifelse(exphist[typ1]>=4&exphist[typ1]<4.5,197,
                                        ifelse(exphist[typ1]>=4.5&exphist[typ1]<5,311,
                                               ifelse(exphist[typ1]>=5&exphist[typ1]<5.5,447,
                                                      ifelse(exphist[typ1]>=5.5&exphist[typ1]<6,610,
                                                             ifelse(exphist[typ1]>=6&exphist[typ1]<6.5,804,
                                                                    ifelse(exphist[typ1]>=6.5&exphist[typ1]<7,1032,
                                                                           ifelse(exphist[typ1]>=7&exphist[typ1]<7.5,1298,
                                                                                  ifelse(exphist[typ1]>=7.5&exphist[typ1]<8,1601,
                                                                                         ifelse(exphist[typ1]>=8&exphist[typ1]<8.5,1936,
                                                                                                ifelse(exphist[typ1]>=8.5&exphist[typ1]<9,2292,
                                                                                                       ifelse(exphist[typ1]>=9&exphist[typ1]<9.5,2635,
                                                                                                              ifelse(exphist[typ1]>=9.5&exphist[typ1]<10,2901,
                                                                                                                     ifelse(exphist[typ1]>=10&exphist[typ1]<10.5,3091,
                                                                                                                            ifelse(exphist[typ1]>=10.5&exphist[typ1]<11,3215,
                                                                                                                                   ifelse(exphist[typ1]>=11&exphist[typ1]<11.5,3281,3300)))))))))))))))))
  
  
  
  powerhist[typ2]<-ifelse(exphist[typ2]<3.5,0,
                          ifelse(exphist[typ2]>=3.5&exphist[typ2]<4,140,
                                 ifelse(exphist[typ2]>=4&exphist[typ2]<4.5,234,
                                        ifelse(exphist[typ2]>=4.5&exphist[typ2]<5,362,
                                               ifelse(exphist[typ2]>=5&exphist[typ2]<5.5,497,
                                                      ifelse(exphist[typ2]>=5.5&exphist[typ2]<6,686,
                                                             ifelse(exphist[typ2]>=6&exphist[typ2]<6.5,880,
                                                                    ifelse(exphist[typ2]>=6.5&exphist[typ2]<7,1142,
                                                                           ifelse(exphist[typ2]>=7&exphist[typ2]<7.5,1407,
                                                                                  ifelse(exphist[typ2]>=7.5&exphist[typ2]<8,1734,
                                                                                         ifelse(exphist[typ2]>=8&exphist[typ2]<8.5,2060,
                                                                                                ifelse(exphist[typ2]>=8.5&exphist[typ2]<9,2386,
                                                                                                       ifelse(exphist[typ2]>=9&exphist[typ2]<9.5,2709,
                                                                                                              ifelse(exphist[typ2]>=9.5&exphist[typ2]<10,2935,
                                                                                                                     ifelse(exphist[typ2]>=10&exphist[typ2]<10.5,3156,
                                                                                                                            ifelse(exphist[typ2]>=10.5&exphist[typ2]<11,3280,
                                                                                                                                   ifelse(exphist[typ2]>=11&exphist[typ2]<11.5,3368,
                                                                                                                                          ifelse(exphist[typ2]>=11.5&exphist[typ2]<12,3415,
                                                                                                                                                 ifelse(exphist[typ2]>=12&exphist[typ2]<12.5,3441,3500)))))))))))))))))))
  
  
  
  powerhist[typ3]<-ifelse(exphist[typ3]<3.5,0,
                          ifelse(exphist[typ3]>=3.5&exphist[typ3]<4,68,
                                 ifelse(exphist[typ3]>=4&exphist[typ3]<4.5,114,
                                        ifelse(exphist[typ3]>=4.5&exphist[typ3]<5,177,
                                               ifelse(exphist[typ3]>=5&exphist[typ3]<5.5,243,
                                                      ifelse(exphist[typ3]>=5.5&exphist[typ3]<6,347,
                                                             ifelse(exphist[typ3]>=6&exphist[typ3]<6.5,452,
                                                                    ifelse(exphist[typ3]>=6.5&exphist[typ3]<7,595,
                                                                           ifelse(exphist[typ3]>=7&exphist[typ3]<7.5,738,
                                                                                  ifelse(exphist[typ3]>=7.5&exphist[typ3]<8,907,
                                                                                         ifelse(exphist[typ3]>=8&exphist[typ3]<8.5,1076,
                                                                                                ifelse(exphist[typ3]>=8.5&exphist[typ3]<9,1307,
                                                                                                       ifelse(exphist[typ3]>=9&exphist[typ3]<9.5,1538,
                                                                                                              ifelse(exphist[typ3]>=9.5&exphist[typ3]<10,1786,
                                                                                                                     ifelse(exphist[typ3]>=10&exphist[typ3]<10.5,2033,
                                                                                                                            ifelse(exphist[typ3]>=10.5&exphist[typ3]<11,2219,
                                                                                                                                   ifelse(exphist[typ3]>=11&exphist[typ3]<11.5,2405,
                                                                                                                                          ifelse(exphist[typ3]>=11.5&exphist[typ3]<12,2535,
                                                                                                                                                 ifelse(exphist[typ3]>=12&exphist[typ3]<12.5,2633,
                                                                                                                                                        ifelse(exphist[typ3]>=12.5&exphist[typ3]<13,2710,2750))))))))))))))))))))
  
  
  powerhist[typ4]<-ifelse(exphist[typ4]<3.5,0,
                          ifelse(exphist[typ4]>=3.5&exphist[typ4]<4,50,
                                 ifelse(exphist[typ4]>=4&exphist[typ4]<4.5,136,
                                        ifelse(exphist[typ4]>=4.5&exphist[typ4]<5,221,
                                               ifelse(exphist[typ4]>=5&exphist[typ4]<5.5,326,
                                                      ifelse(exphist[typ4]>=5.5&exphist[typ4]<6,431,
                                                             ifelse(exphist[typ4]>=6&exphist[typ4]<6.5,576,
                                                                    ifelse(exphist[typ4]>=6.5&exphist[typ4]<7,720,
                                                                           ifelse(exphist[typ4]>=7&exphist[typ4]<7.5,911,
                                                                                  ifelse(exphist[typ4]>=7.5&exphist[typ4]<8,1102,
                                                                                         ifelse(exphist[typ4]>=8&exphist[typ4]<8.5,1339,
                                                                                                ifelse(exphist[typ4]>=8.5&exphist[typ4]<9,1575,
                                                                                                       ifelse(exphist[typ4]>=9&exphist[typ4]<9.5,1575,
                                                                                                              ifelse(exphist[typ4]>=9.5&exphist[typ4]<10,1797,
                                                                                                                     ifelse(exphist[typ4]>=10&exphist[typ4]<10.5,2019,
                                                                                                                            ifelse(exphist[typ4]>=10.5&exphist[typ4]<11,2162,
                                                                                                                                   ifelse(exphist[typ4]>=11&exphist[typ4]<11.5,2304,
                                                                                                                                          ifelse(exphist[typ4]>=11.5&exphist[typ4]<12,2390,
                                                                                                                                                 ifelse(exphist[typ4]>=12&exphist[typ4]<12.5,2458,2500)))))))))))))))))))
  
  
  powerhist[typ5]<-ifelse(exphist[typ5]<3,0,ifelse(exphist[typ5]>=3&exphist[typ5]<3.5,89,
                                                   ifelse(exphist[typ5]>=3.5&exphist[typ5]<4,171,
                                                          ifelse(exphist[typ5]>=4&exphist[typ5]<4.5,269,
                                                                 ifelse(exphist[typ5]>=4.5&exphist[typ5]<5,389,
                                                                        ifelse(exphist[typ5]>=5&exphist[typ5]<5.5,533,
                                                                               ifelse(exphist[typ5]>=5.5&exphist[typ5]<6,704,
                                                                                      ifelse(exphist[typ5]>=6&exphist[typ5]<6.5,906,
                                                                                             ifelse(exphist[typ5]>=6.5&exphist[typ5]<7,1136,
                                                                                                    ifelse(exphist[typ5]>=7&exphist[typ5]<7.5,1400,
                                                                                                           ifelse(exphist[typ5]>=7.5&exphist[typ5]<8,1674,
                                                                                                                  ifelse(exphist[typ5]>=8&exphist[typ5]<8.5,1945,
                                                                                                                         ifelse(exphist[typ5]>=8.5&exphist[typ5]<9,2173,
                                                                                                                                ifelse(exphist[typ5]>=9&exphist[typ5]<9.5,2373,
                                                                                                                                       ifelse(exphist[typ5]>=9.5&exphist[typ5]<10,2518,
                                                                                                                                              ifelse(exphist[typ5]>=10&exphist[typ5]<10.5,2619,
                                                                                                                                                     ifelse(exphist[typ5]>=10.5&exphist[typ5]<11,2696,
                                                                                                                                                            ifelse(exphist[typ5]>=11&exphist[typ5]<11.5,2739,
                                                                                                                                                                   ifelse(exphist[typ5]>=11.5&exphist[typ5]<12,2766,2780)))))))))))))))))))
  
  total<-turb$number*power
  totalhist<-turb$number*powerhist
  diff[,i]<-total-totalhist
}

diff.mean<-rowMeans(diff)
diff.sd<-c()
for(i in 1:41071){
  diff.sd<-c(diff.sd,sd(diff[i,]))
}


i=1
tt<-c()
for(j in 1:273951){
  if(is.na(turb.full$OptTurb[j])){
    tt[j]<-NA
  }
  else{
    tt[j]<-diff.mean[i]
    i=i+1
  }
}

power.mat<-matrix(c(tt),ncol=499,nrow=549)
power.mat<-ifelse(is.na(power.mat),0,power.mat)
agreg_power<-matrix(0,nrow=35,ncol=34)
for (i in 213:246){
  a<-ifelse(sa_lat>lat[i]&sa_lat<lat[i+1],1,0)
  for(j in 344:378){
    b<-ifelse(sa_lon>lon[j]&sa_lon<lon[j+1],1,0)
    wrf_run1<-a*b*power.mat
    agreg_power[j-343,i-212]<-mean(wrf_run1[wrf_run1!=0])
  }
}

power.row<-c(rbind(agreg_power))
power.row<-ifelse(mask.row==1&is.na(power.row),0,power.row)
power.dollar<-ifelse(power.row<=6000,power.row*0.045,power.row*0.075)
cor.c<-data.frame(loc,lac,power.dollar)
colnames(cor.c)<-c("loc","lac","cor")
p<-ggplot(cor.c,aes(x=loc,y=lac,fill=cor))
p+geom_tile() +
  scale_fill_gradient2(low="red",high="blue",midpoint = 0,mid = "grey95")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.8),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),axis.text =element_text(colour="black",size=12),panel.background = element_rect(fill = "white"),
        legend.text=element_text(size=13))

