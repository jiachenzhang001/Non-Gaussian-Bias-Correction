setwd("")
newrm.merra.hist<-readRDS("newrm.merra.hist.RDS")
newrm.cord.hist<-readRDS("newrm.cord.hist.RDS")
newrm.cord.areg<-readRDS("newrm.cord.areg.RDS")
merra_mask<-readRDS("merra_mask.RDS")
mask<-merra_mask[1:35,]
mask.row<-c(rbind(mask))
m.k<-readRDS("mk.RDS")

loc<-rep(seq(34.375,55.625,by=0.625),times=34)
lac<-rep(seq(16,32.5,by=0.5),each=35)
loc<-ifelse(mask.row==0,NA,loc)
lac<-ifelse(mask.row==0,NA,lac)
loc<-as.numeric(na.omit(loc))
lac<-as.numeric(na.omit(lac))


mask.row<-c(rbind(mask))
m.row<-matrix(0,nrow=1190,ncol=9490)
for (l in 1:9490){
  m.row[,l]<-cbind(c(newrm.merra.hist[,,l]))
}

cord.row<-matrix(0,nrow=1190,ncol=9490)
for (l in 1:9490){
  cord.row[,l]<-cbind(c(newrm.cord.hist[1:35,1:34,l]))
}
for(i in 1:9490){
  m.row[,i]<-ifelse(mask.row==0,NA,m.row[,i])
}
for(i in 1:9490){
  cord.row[,i]<-ifelse(mask.row==0,NA,cord.row[,i])
}

a.row<-matrix(0,nrow=1190,ncol=9490)
for (l in 1:9490){
  a.row[,l]<-cbind(c(newrm.cord.areg[,,l]))
}
for(i in 1:9490){
  a.row[,i]<-ifelse(mask.row==0,NA,a.row[,i])
}

t<-!is.na(m.row[,1])
m.row<-m.row[t,]
cord.row<-cord.row[t,]
a.row<-a.row[t,]

pred.m<-matrix(0,nrow=614,ncol=4745)
for (i in 1:614){
  mu<-mean(m.row[i,1:4745])
  osd<-sd(m.row[i,1:4745])
  wsd<-sd(cord.row[i,1:4745])
  wmu<-mean(cord.row[i,1:4745])
  
  wbc<-cord.row[i,4746:9490]+mu-wmu
  pred.m[i,]<-wbc
  
}

KL.dist(t(pred.m),t(m.row[,4746:9490]),k=70)


pred.mv<-matrix(0,nrow=614,ncol=4745)
for (i in 1:614){
  mu<-mean(m.row[i,1:4745])
  osd<-sd(m.row[i,1:4745])
  wsd<-sd(cord.row[i,1:4745])
  wmu<-mean(cord.row[i,1:4745])
  
  wbc<-mu+(osd/wsd)*(cord.row[i,4746:9490]-wmu)
  pred.mv[i,]<-wbc
  
}
KL.dist(t(pred.mv),t(m.row[,4746:9490]),k=70)

#matern covariance function
matern.m<-cbind(loc,lac,m.row)

vg <- variog(coords=matern.m[,1:2],data=matern.m[,3:9492])
h  <- vg$u
vg$v <- rowMeans(vg$v)

m.fit<-variofit(vg,kappa=0.5,fix.kappa = T)
rho<-m.fit$cov.pars[2]

matern.wrf<-cbind(loc,lac,cord.row)

vg.wrf <- variog(coords=matern.wrf[,1:2],data=matern.wrf[,3:9492])
vg.wrf$v<-rowMeans(vg.wrf$v)
wrf.fit<-variofit(vg.wrf,kappa=0.5,fix.kappa = T)
rho.wrf<-wrf.fit$cov.pars[2]

ob<-m.row[,1:4745]
sim<-cord.row[,1:4745]

mu<-rowMeans(m.row[,1:4745])
wmu<-rowMeans(cord.row[,1:4745])

cov.mat<-matrix(0,nrow=614,ncol=614)
cov.wrf.mat<-matrix(0,nrow=614,ncol=614)
for(n in 1:614){
  d<-sqrt((lac-lac[n])^2+(loc-loc[n])^2)
  cov.mat[,n]<-matern(d, rho, 0.5)
  cov.wrf.mat[,n]<-matern(d,rho.wrf,0.5)
}


w<-solve(chol(cov.wrf.mat))
o<-chol(cov.mat)

sigma<-o%*%w
correct<-sigma%*%(cord.row[,4746:9490]-wmu)
p.row<-mu+correct

KL.dist(t(p.row),t(m.row[,4746:9490]),k=70)


#one transformation parameter
lambda<-seq(0.1,1.5,by=0.05)


kld.lamb<-c()  
for(g in 1:29){
  lam.m<-matrix(0,nrow=614,ncol=9490)
  for(i in 1:614){
    lam.m[i,]<-ifelse(m.row[i,]>=0,((m.row[i,]+1)^lambda[g]-1)/lambda[g],-((-m.row[i,]+1)^(2-lambda[g])-1)/(2-lambda[g]))
  }
  matern.m<-cbind(loc,lac,lam.m)
  
  vg <- variog(coords=matern.m[,1:2],data=matern.m[,3:9492])
  vg$v<-rowMeans(vg$v)
  m.fit<-variofit(vg,kappa=0.5,fix.kappa = T)
  rho<-m.fit$cov.pars[2]

  for(gw in 1:29){
    lam.cord<-matrix(0,nrow=614,ncol=9490)
    for(i in 1:614){
      lam.cord[i,]<-ifelse(cord.row[i,]>=0,((cord.row[i,]+1)^lambda[gw]-1)/lambda[gw],-((-cord.row[i,]+1)^(2-lambda[gw])-1)/(2-lambda[gw]))
    }  
    matern.wrf<-cbind(loc,lac,lam.cord)
    vg.wrf <- variog(coords=matern.wrf[,1:2],data=matern.wrf[,3:9492])

    vg.wrf$v<-rowMeans(vg.wrf$v)
    wrf.fit<-variofit(vg.wrf,kappa=0.5,fix.kappa = T)
    rho.wrf<-wrf.fit$cov.pars[2]
    
    ob<-lam.m[,1:4745]
    sim<-lam.cord[,1:4745]
    
    mu<-rowMeans(lam.m[,1:4745])
    wmu<-rowMeans(lam.cord[,1:4745])
    
    cov.mat<-matrix(0,nrow=614,ncol=614)
    cov.wrf.mat<-matrix(0,nrow=614,ncol=614)
    for(n in 1:614){
      d<-sqrt((lac-lac[n])^2+(loc-loc[n])^2)
      cov.mat[,n]<-matern(d, rho, 0.5)
      cov.wrf.mat[,n]<-matern(d,rho.wrf,0.5)
    }
    
    w<-solve(chol(cov.wrf.mat))
    o<-chol(cov.mat)
    
    sigma<-o%*%w
    correct<-sigma%*%(lam.cord[,4746:9490]-wmu)
    p.row<-mu+correct
    
    m.func<-function(x){
      y<-ifelse(x>=0,((x+1)^lambda[g]-1)/lambda[g],-((-x+1)^(2-lambda[g])-1)/(2-lambda[g]))
      return(y)
    }
    m.inv<-inverse(m.func)
    for(r in 1:614){
      p.row[r,]<-m.inv(p.row[r,])
    
    k<-KL.dist(t(p.row),t(m.row[,4746:9490]),k=70)
    kld.lamb<-c(kld.lamb,k[70])
    
  }
  
}

#clusterwise transformation parameter
gam.m<-c()
gam.wrf<-c()

for(c in c(1:20)){
  gam.count<-c()
  gam.min<-c()
  m.mat<-m.row[m.k==c,]
  wrf.mat<-cord.row[m.k==c,]      
  la.c<-lac[m.k==c]
  lo.c<-loc[m.k==c]
  for(g in 1:29){
    kldwf<-c()
    mr<-nrow(m.mat)
    lam.m<-matrix(0,nrow=mr,ncol=9490)
    for(i in 1:mr){
      lam.m[i,]<-ifelse(m.mat[i,]>=0,((m.mat[i,]+1)^lambda[g]-1)/lambda[g],-((-m.mat[i,]+1)^(2-lambda[g])-1)/(2-lambda[g]))
    }
    matern.m<-cbind(lo.c,la.c,lam.m)
    
    vg <- variog(coords=matern.m[,1:2],data=matern.m[,3:9492])
    vg$v  <- rowMeans(vg$v)
    m.fit<-variofit(vg,kappa=0.5,fix.kappa = T)
    rho<-m.fit$cov.pars[2]
    for(gw in 1:29){
      lam.wrf<-matrix(0,nrow=mr,ncol=9490)
      for(j in 1:mr){
        lam.wrf[j,]<-ifelse(wrf.mat[j,]>=0,((wrf.mat[j,]+1)^lambda[gw]-1)/lambda[gw],-((-wrf.mat[j,]+1)^(2-lambda[gw])-1)/(2-lambda[gw]))
      }  
      matern.wrf<-cbind(lo.c,la.c,lam.wrf)
      vg.wrf <- variog(coords=matern.wrf[,1:2],data=matern.wrf[,3:9492])
      vg.wrf$v  <- rowMeans(vg.wrf$v)
      
      fit.wrf <- variofit(vg.wrf,kappa=0.5,fix.kappa = T)
      rho.wrf<-fit.wrf$cov.pars[2]
      
      ob<-lam.m[,1:4745]
      sim<-lam.wrf[,1:4745]
      
      mu<-rowMeans(lam.m[,1:4745])
      wmu<-rowMeans(lam.wrf[,1:4745])
      
      cov.mat<-matrix(0,nrow=mr,ncol=mr)
      cov.wrf.mat<-matrix(0,nrow=mr,ncol=mr)
      for(n in 1:mr){
        d<-sqrt((la.c-la.c[n])^2+(lo.c-lo.c[n])^2)
        cov.mat[,n]<-matern(d,rho,0.5)
        cov.wrf.mat[,n]<-matern(d,rho.wrf,0.5)
      }
      w<-solve(chol(cov.wrf.mat))
      o<-chol(cov.mat)
      
      sigma<-o%*%w
      correct<-sigma%*%(lam.wrf[,4746:9490]-wmu)
      p<-mu+correct
      
      m.func<-function(x){
        y<-ifelse(x>=0,((x+1)^lambda[g]-1)/lambda[g],-((-x+1)^(2-lambda[g])-1)/(2-lambda[g]))
        return(y)
      }
      m.inv<-inverse(m.func)
      for(r in 1:614){
        p[r,]<-m.inv(p[r,])
      }
      
      k<-KL.dist(t(p),t(m.row[,4746:9490]),k=70)
      kldwf<-c(kldwf,k[70])
      
    }
    gam.count<-c(gam.count,which.min(kldwf))
    gam.min<-c(gam.min,kldwf[which.min(kldwf)])
  }
  gam.m<-c(gam.m,lambda[which.min(gam.min)])
  gam.wrf<-c(gam.wrf,lambda[gam.count[which.min(gam.min)]])
}



#nonstationary matern

fit<-NSconvo_fit(coords = cbind(loc,lac),data=m.row ,fit.radius = 5, N.mc = 4)
covmat<-fit$Cov.mat



fit<-NSconvo_fit(coords = cbind(loc,lac),data=cord.row ,fit.radius = 5, N.mc = 4)
cordmat<-fit$Cov.mat



ob<-m.row[,1:4745]
sim<-cord.row[,1:4745]

mu<-rowMeans(m.row[,1:4745])
wmu<-rowMeans(cord.row[,1:4745])


w<-solve(chol(cordmat))
o<-chol(covmat)
sigma<-o%*%w
correct<-sigma%*%(cord.row[,4746:9490]-wmu)
p.row<-mu+correct


KL.dist(t(p.row),t(m.row[,4746:9490]),k=70)

#nonstationary one parameter trasfoamtion
kld.lamb<-c()  
for(g in 1:29){
  lam.m<-matrix(0,nrow=614,ncol=9490)
  for(i in 1:614){
    lam.m[i,]<-ifelse(m.row[i,]>=0,((m.row[i,]+1)^lambda[g]-1)/lambda[g],-((-m.row[i,]+1)^(2-lambda[g])-1)/(2-lambda[g]))
  }

  fit<-NSconvo_fit(coords = cbind(loc,lac),data=lam.m ,fit.radius = 5, N.mc = 4)
  covmat<-fit$Cov.mat
  
  for(gw in 1:29){
    lam.cord<-matrix(0,nrow=614,ncol=9490)
    for(i in 1:614){
      lam.cord[i,]<-ifelse(cord.row[i,]>=0,((cord.row[i,]+1)^lambda[gw]-1)/lambda[gw],-((-cord.row[i,]+1)^(2-lambda[gw])-1)/(2-lambda[gw]))
    }  
  
    fit<-NSconvo_fit(coords = cbind(loc,lac),data=lam.cord ,fit.radius = 5, N.mc = 4)
    cordmat<-fit$Cov.mat
    

    
    ob<-lam.m[,1:4745]
    sim<-lam.cord[,1:4745]
    
    mu<-rowMeans(lam.m[,1:4745])
    wmu<-rowMeans(lam.cord[,1:4745])
    
    w<-solve(chol(cordmat))
    o<-chol(covmat)
    
    sigma<-o%*%w
    correct<-sigma%*%(cord.row[,4746:9490]-wmu)
    p.row<-mu+correct
    
    m.func<-function(x){
      y<-ifelse(x>=0,((x+1)^lambda[g]-1)/lambda[g],-((-x+1)^(2-lambda[g])-1)/(2-lambda[g]))
      return(y)
    }
    m.inv<-inverse(m.func)
    for(r in 1:614){
      p.row[r,]<-m.inv(p.row[r,])
    }
    
    k<-KL.dist(t(p.row),t(m.row[,4746:9490]),k=70)
    kld.lamb<-c(kld.lamb,k[70])
    
  }
  
}

which.min(kld.lamb)



#cluster-wise
gam.m<-c()
gam.cord<-c()

for(c in c(1:20)){
  gam.count<-c()
  gam.min<-c()
  m.mat<-m.row[m.k==c,]
  wrf.mat<-cord.row[m.k==c,]      
  la.c<-lac[m.k==c]
  lo.c<-loc[m.k==c]
  for(g in 1:29){
    kldwf<-c()
    mr<-nrow(m.mat)
    lam.m<-matrix(0,nrow=mr,ncol=9490)
    for(i in 1:mr){
      lam.m[i,]<-ifelse(m.mat[i,]>=0,((m.mat[i,]+1)^lambda[g]-1)/lambda[g],-((-m.mat[i,]+1)^(2-lambda[g])-1)/(2-lambda[g]))
    }
  
    fit<-NSconvo_fit(coords = cbind(loc,lac),data=lam.m ,fit.radius = 5, N.mc = 4)
    covmat<-fit$Cov.mat
    
    mean.cov<-matrix(0,nrow=614,ncol=614)
    for(i in 1:614){
      for(j in 1:614){
        mean.cov[i,j]<-mean(covmat[i,j,])
      }
    }
    for(gw in 1:29){
      lam.wrf<-matrix(0,nrow=mr,ncol=9490)
      for(j in 1:mr){
        lam.wrf[j,]<-ifelse(wrf.mat[j,]>=0,((wrf.mat[j,]+1)^lambda[gw]-1)/lambda[gw],-((-wrf.mat[j,]+1)^(2-lambda[gw])-1)/(2-lambda[gw]))
      }  
     
      fit<-NSconvo_fit(coords = cbind(loc,lac),data=lam.cord ,fit.radius = 5, N.mc = 4)
      cordmat[,,i]<-fit$Cov.mat
     
   
      
      ob<-lam.m[,1:4745]
      sim<-lam.cord[,1:4745]
      
      mu<-rowMeans(lam.m[,1:4745])
      wmu<-rowMeans(lam.cord[,1:4745])
      

      w<-solve(chol(cordmat))
      o<-chol(covmat)
      
      sigma<-o%*%w
      correct<-sigma%*%(lam.wrf[,4746:9490]-wmu)
      p<-mu+correct
      
      m.func<-function(x){
        y<-ifelse(x>=0,((x+1)^lambda[g]-1)/lambda[g],-((-x+1)^(2-lambda[g])-1)/(2-lambda[g]))
        return(y)
      }
      m.inv<-inverse(m.func)
      for(r in 1:614){
        p[r,]<-m.inv(p[r,])
      }
      
      k<-KL.dist(t(p),t(m.row[,4746:9490]),k=70)
      kldwf<-c(kldwf,k[70])
      
    }
    gam.count<-c(gam.count,which.min(kldwf))
    gam.min<-c(gam.min,kldwf[which.min(kldwf)])
    k1<-c(k1,kldwf[which.min(kldwf)])
  }
  gam.m<-c(gam.m,lambda[which.min(gam.min)])
  gam.cord<-c(gam.cord,lambda[gam.count[which.min(gam.min)]])
  k2<-c(k2,k1[which.min(k1)])
}


yj.m<-matrix(0,nrow=614,ncol=4745)
yj.cord<-matrix(0,nrow=614,ncol=4745)
for(i in 1:20){
  num<-sum(m.k==i)
  rownum<-c(1:614)[m.k==i]
  for(j in 1:num){
    yj.m[rownum[j],]<-ifelse(m.row[rownum[j],]>=0,((m.row[rownum[j],]+1)^gam.m[i]-1)/gam.m[i],-((-m.row[rownum[j],]+1)^(2-gam.m[i])-1)/(2-gam.m[i]))
    yj.cord[rownum[j],]<-ifelse(cord.row[rownum[j],]>=0,((cord.row[rownum[j],]+1)^gam.cord[i]-1)/gam.cord[i],-((-cord.row[rownum[j],]+1)^(2-gam.cord[i])-1)/(2-gam.cord[i]))
    }
}


fit<-NSconvo_fit(coords = cbind(loc,lac),data=yj.m ,fit.radius = 5, N.mc = 4)
yjcov.m<-fit$Cov.mat



fit<-NSconvo_fit(coords = cbind(loc,lac),data=yj.cord ,fit.radius = 5, N.mc = 4)
yjcov.cord<-fit$Cov.mat



ob<-yj.m[,1:4745]
sim<-yj.cord[,1:4745]

mu<-rowMeans(yj.m[,1:4745])
wmu<-rowMeans(yj.cord[,1:4745])


w<-solve(chol(yjcov.cord))
o<-chol(yjcov.m)
sigma<-o%*%w
correct<-sigma%*%(yj.cord[,4746:9490]-wmu)
p.row<-mu+correct

pred<-c()
for(c in 1:20){
  m.func<-function(x){
    y<-ifelse(x>=0,((x+1)^gam.m[c]-1)/gam.m[c],-((-x+1)^(2-gam.m[c])-1)/(2-gam.m[c]))
    return(y)
  }
  m.inv<-inverse(m.func)
  mat<-p.row[m.k==c,]
  nr<-dim(mat)[1]
  for(i in 1:nr){
    pred<-rbind(pred,m.inv(mat[nr,]))
  }
}
KL.dist(t(pred),t(m.row[,4746:9490]),k=70)
