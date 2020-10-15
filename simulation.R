library(convoSPAT)
library(FNN)
library(geoR)

#cross-validation

m.row<-readRDS("m.row.RDS")
cord.row<-readRDS("cord.row.RDS")

kld.m<-c()
kld.mv<-c()
for(k in 1:1000){
  sam<-sample(c(1:614),300,replace=F)
  pred.m<-c()
  sim.m<-m.row[sam,]
  sim.cord<-cord.row[sam,]
  for (i in 1:300){
    mu<-mean(sim.m[i,1:4745])
    osd<-sd(sim.m[i,1:4745])
    wsd<-sd(sim.cord[i,1:4745])
    wmu<-mean(sim.cord[i,1:4745])
    
    wbc<-sim.cord[i,4746:9490]+mu-wmu
    wbc.mv<-mu+(osd/wsd)*(sim.cord[i,7301:9490]-wmu)
    pred.m<-rbind(pred.m,wbc)
    pred.mv<-rbind(pred.mv,wbc)
  }
  kl<-KL.dist(t(pred.m),t(sim.m[,4746:9490]),k=70)
  klmv<-KL.dist(t(pred.mv),t(sim.m[,4746:9490]),k=70)
  kld.m<-c(kld.m,kl[70])
  kld.mv<-c(kld.mv,klmv[70])
}

kld.mn<-c()
for(k in 1:1000){
  sam<-sample(614,size=300,replace=F)
  pred.m<-c()
  sim.m<-m.row[sam,]
  sim.wrf<-wrf.row[sam,]
  lo.sim<-lo[sam]
  la.sim<-la[sam]
  
  fit<-NSconvo_fit(coords = cbind(lo.sim,la.sim),data=sim.m ,fit.radius = 5, N.mc = 4)
  covmat<-fit$Cov.mat
    
  fit<-NSconvo_fit(coords = cbind(lo.sim,la.sim),data=sim.cord ,fit.radius = 5, N.mc = 4)
  cordmat<-fit$Cov.mat
      
  ob<-sim.m[,1:4745]
  sim<-sim.cord[,1:4745]
      
  mu<-rowMeans(sim.m[,1:4745])
  wmu<-rowMeans(sim.cord[,1:4745])
      
      
  w<-solve(chol(cordmat))
  o<-chol(covmat)
      
  sigma<-o%*%w
  correct<-sigma%*%(sim.cord[,4746:9490]-wmu)
  p.row<-mu+correct
      
  k<-KL.dist(t(p.row),t(sim.m[,4746:9490]),k=70)
  kld.mn<-c(kld.mn,kl[70])
 
}

gam.m<-c()
gam.wrf<-c()
kld.t1<-c()
for(k in 1:1000){
  sam<-sample(614,size=300,replace=F)
  pred.m<-c()
  sim.m<-m.row[sam,]
  sim.wrf<-wrf.row[sam,]
  lo.sim<-lo[sam]
  la.sim<-la[sam]
  
  gam.count<-c()
  gam.min<-c()
  for(g in 1:29){
    kldwf<-c()
    lam.m<-matrix(0,nrow=rn,ncol=9490)
    for(i in 1:rn){
      lam.m[i,]<-ifelse(sim.m[i,]>=0,((sim.m[i,]+1)^lambda[g]-1)/lambda[g],-((-sim.m[i,]+1)^(2-lambda[g])-1)/(2-lambda[g]))
    }

    fit<-NSconvo_fit(coords = cbind(lo.clust,la.clust),data=lam.m ,fit.radius = 5, N.mc = 4)
    covmat<-fit$Cov.mat
    

    for(gw in 1:29){
      lam.cord<-matrix(0,nrow=rn,ncol=9490)
      for(i in 1:rn){
        lam.cord[i,]<-ifelse(sim.cord[i,]>=0,((sim.cord[i,]+1)^lambda[gw]-1)/lambda[gw],-((-sim.cord[i,]+1)^(2-lambda[gw])-1)/(2-lambda[gw]))
      }  
  
      fit<-NSconvo_fit(coords = cbind(lo.clust,la.clust),data=lam.cord ,fit.radius = 5, N.mc = 4)
      cordmat<-fit$Cov.mat
      
  
      ob<-lam.m[,1:4745]
      sim<-lam.cord[,1:4745]
      
      mu<-rowMeans(lam.m[,1:4745])
      wmu<-rowMeans(lam.cord[,1:4745])
      
      
      w<-solve(chol(cordmat))
      o<-chol(covmat)
      
      sigma<-o%*%w
      correct<-sigma%*%(lam.cord[,4746:9490]-wmu)
      p.row<-mu+correct
      
      k<-KL.dist(t(p.row),t(lam.m[,4746:9490]),k=70)
      kldwf<-c(kldwf,kl[70])
    }
    gam.count<-c(gam.count,which.min(kldwf))
    gam.min<-c(gam.min,kldwf[which.min(kldwf)])
  }
  gam.m<-c(gam.m,gam[which.min(gam.min)])
  gam.cord<-c(gam.cord,gam[gam.count[which.min(gam.min)]])
  kld.t1<-c(kld.t1,gam.min[which.min(gam.min)])
}

kld.tc<-c()
for(k in 1:1000){
  lo.sim<-c()
  la.sim<-c()
  gam.m<-c()
  gam.cord<-c()
  mk.clust<-c()
  mat.m<-c()
  mat.cord<-c()
  for(j in 1:20){
    pos<-c(1:614)[m.k==j]
    sam<-sample(cluster.size[j],size=sam.num[j],replace=F)
    sim.m<-m.row[pos[sam],]
    mat.m<-rbind(mat.m,sim.m)
    
    sim.cord<-cord.row[pos[sam],]
    mat.cord<-rbind(mat.cord,sim.cord)
    
    lo.sim<-c(lo.sim,lo[pos[sam]])
    la.sim<-c(la.sim,la[pos[sam]])
    la.clust<-la[pos[sam]]
    lo.clust<-lo[pos[sam]]
    
    mk.clust<-c(mk.clust,m.k[pos[sam]])
    
    gam.count<-c()
    gam.min<-c()
    rn<-nrow(sim.m)
    for(g in 1:29){
      kldwf<-c()
      lam.m<-matrix(0,nrow=rn,ncol=9490)
      for(i in 1:rn){
        lam.m[i,]<-ifelse(sim.m[i,]>=0,((sim.m[i,]+1)^lambda[g]-1)/lambda[g],-((-sim.m[i,]+1)^(2-lambda[g])-1)/(2-lambda[g]))
      }
      covmat<-array(0,c(rn,rn,4745))
      
      fit<-NSconvo_fit(coords = cbind(lo.clust,la.clust),data=lam.m ,fit.radius = 5, N.mc = 4)
      covmat<-fit$Cov.mat
 
      
      for(gw in 1:29){
        lam.cord<-matrix(0,nrow=rn,ncol=9490)
        for(i in 1:rn){
          lam.cord[i,]<-ifelse(sim.cord[i,]>=0,((sim.cord[i,]+1)^lambda[gw]-1)/lambda[gw],-((-sim.cord[i,]+1)^(2-lambda[gw])-1)/(2-lambda[gw]))
        }  
        cordmat<-array(0,c(rn,rn,4745))
        
        fit<-NSconvo_fit(coords = cbind(lo.clust,la.clust),data=lam.cord ,fit.radius = 5, N.mc = 4)
        cordmat<-fit$Cov.mat

        
        ob<-lam.m[,1:4745]
        sim<-lam.cord[,1:4745]
        
        mu<-rowMeans(lam.m[,1:4745])
        wmu<-rowMeans(lam.cord[,1:4745])
        
        w<-solve(chol(cordmat))
        o<-chol(covmat)
        
        sigma<-o%*%w
        correct<-sigma%*%(lam.cord[,4746:9490]-wmu)
        p.row<-mu+correct
        
        k<-KL.dist(t(p.row),t(lam.m[,4746:9490]),k=70)
        kldwf<-c(kldwf,kl[70])
      }
      gam.count<-c(gam.count,which.min(kldwf))
      gam.min<-c(gam.min,kldwf[which.min(kldwf)])
    }
    gam.m<-c(gam.m,gam[which.min(gam.min)])
    gam.cord<-c(gam.cord,gam[gam.count[which.min(gam.min)]])
  }
  m.yj<-c()
  cord.yj<-c()
  
  for(c in 1:20){
    m.mat<-mat.m[c(1:301)[mk.clust==c],]
    cord.mat<-mat.cord[c(1:301)[mk.clust==c],]      
    mr<-nrow(m.mat)
    for(p in 1:mr){
      m.mat[p,]<-ifelse(m.mat[p,]>=0,((m.row[p,]+1)^gam.m[c]-1)/gam.m[c],-((-m.mat[p,]+1)^(2-gam.m[c])-1)/(2-gam.m[c]))
    }
    m.yj<-rbind(m.yj,m.mat)
    
    for(l in 1:mr){
      cord.mat[l,]<-ifelse(cord.mat[l,]>=0,((cord.mat[l,]+1)^gam.cord[c]-1)/gam.cord[c],-((-cord.mat[i,]+1)^(2-gam.cord[c])-1)/(2-gam.cord[c]))
    }  
    cord.yj<-rbind(cord.yj,cord.mat)
  }

  fit<-NSconvo_fit(coords = cbind(lo.sim,la.sim),data=m.yj ,fit.radius = 5, N.mc = 4)
  covmat<-fit$Cov.mat
  


  fit<-NSconvo_fit(coords = cbind(lo.sim,la.sim),data=cord.yj ,fit.radius = 5, N.mc = 4)
  cordmat<-fit$Cov.mat


  
  ob<-m.yj[,1:4745]
  sim<-cord.yj[,1:4745]
  
  mu<-rowMeans(m.yj[,1:4745])
  wmu<-rowMeans(cord.yj[,1:4745])
  
  
  w<-solve(chol(cordmat))
  o<-chol(covmat)
  
  sigma<-o%*%w
  correct<-sigma%*%(ob-wmu)
  p<-mu+correct
  kld<-KL.dist(p,cord.yj[,4746:9490],k=70)
  kld.tc<-c(kld.tc,kld[70])
}


#simulation

lat<-rep(seq(from=0,by=0.2,length.out = 10),times=20)
lon<-rep(seq(from=0,by=0.2,length.out = 20),each=10)
cls<-c(rep(rep(1:2,each=5),times=5),rep(rep(3:4,each=5),times=5),rep(rep(5:6,each=5),times=5),rep(rep(7:8,each=5),times=5))


skt.array<-array(0,dim=c(200,100,1000))
glg.array<-array(0,dim=c(200,100,1000))
for(i in 1:1000){

  skt<-matrix(0,ncol=100,nrow=200)
  glg<-matrix(0,ncol=100,nrow=200)
  for(t in 1:100){
    ur=grf(200,cov.pars = c(1.5,0.2),cov.model = "exponential",grid = cbind(lon,lat))
    eta=grf(200,cov.pars = c(1.5,0.5),cov.model = "exponential",grid = cbind(lon,lat))
    rgam=rgamma(200,4,4)
    ys=(0.5*abs(ur)+eta$data)/sqrt(rgam)
    skt[,t]<-ys
    
    xi=rlnorm(200,meanlog = 0,sdlog = 0.7)
    eps=grf(200,cov.pars = c(1.5,0.5),cov.model = "exponential",grid = cbind(lon,lat))
    g=(eta$data/sqrt(xi))+eps
    glg[,t]<-g
  }
  
  skt.array[,,i]<-skt
  glg.array[,,i]<-glg
  
}



kld.m<-c()
kld.mv<-c()
for(i in 1:1000){
  for (j in 1:200){
    mu<-mean(skt.array[j,1:50,i])
    osd<-sd(skt.array[j,1:50,i])
    wsd<-sd(glg.array[j,1:50,i])
    wmu<-mean(glg.array[j,1:50,i])
    
    wbc<-glg.array[j,51:100,i]+mu-wmu
    wbc.mv<-mu+(osd/wsd)*(glg.array[j,51:100,i]-wmu)
    pred.m<-rbind(pred.m,wbc)
    pred.mv<-rbind(pred.mv,wbc)
  }
  kl<-KL.dist(t(pred.m),t(skt.array[j,51:100,i]),k=7)
  klmv<-KL.dist(t(pred.mv),t(skt.array[j,51:100,i]),k=7)
  kld.m<-c(kld.m,kl[7])
  kld.mv<-c(kld.mv,klmv[7])
  
  
}




kld.mn<-c()
for(i in 1:1000){
  
    fit<-NSconvo_fit(coords = cbind(lo.sim,la.sim),data=skt.array[,1:50,i] ,fit.radius = 5, N.mc = 4)
    covmat<-fit$Cov.mat
    
    fit<-NSconvo_fit(coords = cbind(lo.sim,la.sim),data=glg.array[,1:50,i] ,fit.radius = 5, N.mc = 4)
    cordmat<-fit$Cov.mat
    
    ob<-skt.array[,1:50,i]
    sim<-glg.array[,1:50,i]
    
    mu<-rowMeans(skt.array[,1:50,i])
    wmu<-rowMeans(glg.array[,1:50,i])
    
    
    w<-solve(chol(cordmat))
    o<-chol(covmat)
    
    sigma<-o%*%w
    correct<-sigma%*%(glg.array[,51:100,i]-wmu)
    p.row<-mu+correct
    
    k<-KL.dist(t(p.row),t(skt.array[,51:100,i]),k=7)
    kld.mn<-c(kld.mn,kl[7])
  
}



kld.t1<-c()
for(i in 1:1000){
  gam.count<-c()
  gam.min<-c()
  for(g in 1:29){
    kldwf<-c()
    lam.m<-matrix(0,nrow=200,ncol=100)
    for(t in 1:200){
      lam.m[t,]<-ifelse(skt.array[t,,i]>=0,((skt.array[t,,i]+1)^lambda[g]-1)/lambda[g],-((-skt.array[t,,i]+1)^(2-lambda[g])-1)/(2-lambda[g]))
    }
    
    fit<-NSconvo_fit(coords = cbind(lon,lat),data=lam.m ,fit.radius = 5, N.mc = 4)
    covmat<-fit$Cov.mat
    
    
    for(gw in 1:29){
      lam.cord<-matrix(0,nrow=200,ncol=100)
      for(p in 1:200){
        lam.cord[p,]<-ifelse(glg.array[p,,i]>=0,((glg.array[p,,i]+1)^lambda[gw]-1)/lambda[gw],-((-glg.array[p,,i]+1)^(2-lambda[gw])-1)/(2-lambda[gw]))
      }  
      
      fit<-NSconvo_fit(coords = cbind(lon,lat),data=lam.cord ,fit.radius = 5, N.mc = 4)
      cordmat<-fit$Cov.mat
      
      
      ob<-lam.m[,1:50]
      sim<-lam.cord[,1:50]
      
      mu<-rowMeans(lam.m[,1:50])
      wmu<-rowMeans(lam.cord[,1:50])
      
      
      w<-solve(chol(cordmat))
      o<-chol(covmat)
      
      sigma<-o%*%w
      correct<-sigma%*%(lam.cord[,51:100]-wmu)
      p.row<-mu+correct
      
      k<-KL.dist(t(p.row),t(lam.m[,51:100]),k=7)
      kldwf<-c(kldwf,kl[7])
    }
    gam.count<-c(gam.count,which.min(kldwf))
    gam.min<-c(gam.min,kldwf[which.min(kldwf)])
  }
  gam.m<-c(gam.m,gam[which.min(gam.min)])
  gam.cord<-c(gam.cord,gam[gam.count[which.min(gam.min)]])
  kld.t1<-c(kld.t1,gam.min[which.min(gam.min)])
 
  
}


kld.tc<-c()
for(k in 1:1000){
  lo.sim<-c()
  la.sim<-c()
  gam.m<-c()
  gam.cord<-c()
  mk.clust<-c()
  mat.m<-c()
  mat.cord<-c()
  for(j in 1:8){
    pos<-c(1:200)[cls==j]
    sam<-sample(cluster.size[j],size=sam.num[j],replace=F)
    sim.m<-m.row[pos[sam],]
    mat.m<-rbind(mat.m,sim.m)
    
    sim.cord<-cord.row[pos[sam],]
    mat.cord<-rbind(mat.cord,sim.cord)
    
    lo.sim<-c(lo.sim,lo[pos[sam]])
    la.sim<-c(la.sim,la[pos[sam]])
    la.clust<-la[pos[sam]]
    lo.clust<-lo[pos[sam]]
    
    mk.clust<-c(mk.clust,m.k[pos[sam]])
    
    gam.count<-c()
    gam.min<-c()
    rn<-nrow(sim.m)
    for(g in 1:29){
      kldwf<-c()
      lam.m<-matrix(0,nrow=rn,ncol=9490)
      for(i in 1:rn){
        lam.m[i,]<-ifelse(sim.m[i,]>=0,((sim.m[i,]+1)^lambda[g]-1)/lambda[g],-((-sim.m[i,]+1)^(2-lambda[g])-1)/(2-lambda[g]))
      }
      covmat<-array(0,c(rn,rn,4745))
      
      fit<-NSconvo_fit(coords = cbind(lo.clust,la.clust),data=lam.m ,fit.radius = 5, N.mc = 4)
      covmat<-fit$Cov.mat
      
      
      for(gw in 1:29){
        lam.cord<-matrix(0,nrow=rn,ncol=9490)
        for(i in 1:rn){
          lam.cord[i,]<-ifelse(sim.cord[i,]>=0,((sim.cord[i,]+1)^lambda[gw]-1)/lambda[gw],-((-sim.cord[i,]+1)^(2-lambda[gw])-1)/(2-lambda[gw]))
        }  
        cordmat<-array(0,c(rn,rn,4745))
        
        fit<-NSconvo_fit(coords = cbind(lo.clust,la.clust),data=lam.cord ,fit.radius = 5, N.mc = 4)
        cordmat<-fit$Cov.mat
        
        
        ob<-lam.m[,1:4745]
        sim<-lam.cord[,1:4745]
        
        mu<-rowMeans(lam.m[,1:4745])
        wmu<-rowMeans(lam.cord[,1:4745])
        
        w<-solve(chol(cordmat))
        o<-chol(covmat)
        
        sigma<-o%*%w
        correct<-sigma%*%(lam.cord[,4746:9490]-wmu)
        p.row<-mu+correct
        
        k<-KL.dist(t(p.row),t(lam.m[,4746:9490]),k=70)
        kldwf<-c(kldwf,kl[70])
      }
      gam.count<-c(gam.count,which.min(kldwf))
      gam.min<-c(gam.min,kldwf[which.min(kldwf)])
    }
    gam.m<-c(gam.m,gam[which.min(gam.min)])
    gam.cord<-c(gam.cord,gam[gam.count[which.min(gam.min)]])
  }
  m.yj<-c()
  cord.yj<-c()
  
  for(c in 1:20){
    m.mat<-mat.m[c(1:301)[mk.clust==c],]
    cord.mat<-mat.cord[c(1:301)[mk.clust==c],]      
    mr<-nrow(m.mat)
    for(p in 1:mr){
      m.mat[p,]<-ifelse(m.mat[p,]>=0,((m.row[p,]+1)^gam.m[c]-1)/gam.m[c],-((-m.mat[p,]+1)^(2-gam.m[c])-1)/(2-gam.m[c]))
    }
    m.yj<-rbind(m.yj,m.mat)
    
    for(l in 1:mr){
      cord.mat[l,]<-ifelse(cord.mat[l,]>=0,((cord.mat[l,]+1)^gam.cord[c]-1)/gam.cord[c],-((-cord.mat[i,]+1)^(2-gam.cord[c])-1)/(2-gam.cord[c]))
    }  
    cord.yj<-rbind(cord.yj,cord.mat)
  }
  
  fit<-NSconvo_fit(coords = cbind(lo.sim,la.sim),data=m.yj ,fit.radius = 5, N.mc = 4)
  covmat<-fit$Cov.mat
  
  
  
  fit<-NSconvo_fit(coords = cbind(lo.sim,la.sim),data=cord.yj ,fit.radius = 5, N.mc = 4)
  cordmat<-fit$Cov.mat
  
  
  
  ob<-m.yj[,1:4745]
  sim<-cord.yj[,1:4745]
  
  mu<-rowMeans(m.yj[,1:4745])
  wmu<-rowMeans(cord.yj[,1:4745])
  
  
  w<-solve(chol(cordmat))
  o<-chol(covmat)
  
  sigma<-o%*%w
  correct<-sigma%*%(ob-wmu)
  p<-mu+correct
  kld<-KL.dist(p,cord.yj[,4746:9490],k=70)
  kld.tc<-c(kld.tc,kld[70])
}

