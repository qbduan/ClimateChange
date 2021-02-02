
QR_Spatial_2<-function(y,s,X,L=4,
                     mn_range=-2, sd_range=1, iniallist = list(),
                     a.var=0.1, b.var=0.1,
                     sd.mean=1000,
                     tau=seq(0.01,0.99,0.01),
                     iters=10000,burn=1000,update=10){
  
  
  library(emulator)
  library(fields)
  
  ns<-ncol(y)
  ntau<-length(tau)
  nt<-nrow(y)
  p<-dim(X)[3]
  d<-as.matrix(dist(s,upper=T,diag=T))
  z<-y
  
  ################################################:
  #######     UPDATE INITIAL VALUES    ###########:
  ################################################:
  
  med<-matrix(0,p,ns)
  theta<-array(0,c(p,ns,L))
  
  if(length(iniallist)==0) {
    for (j in 1:ns) {
      fff <- lm(y[, j] ~ X[, j, ] - 1)
      med[, j] <- fff$coef
      theta[1, j, ] <- sd(fff$residuals)
    }
    range <- rangez <- exp(mn_range)
    ratioz <- 0.5
    sdy<-sd(as.vector(y))
    MHmed<-rep(sdy/10,p)
    MHtheta<-matrix(sdy/10,p,L)
    
    tau_med<-rep(1,p)
    tau_theta<-matrix(1,p,L)
    mn_med<-rowMeans(med)
    mn_theta<-apply(theta,c(1,3),mean)
    
    P<-solve(covexp_cpp(d,range,1))      
    Pz<-solve(covexp_cpp(d,rangez,ratioz))      
    
    curll<-rep(0,ns)
    for(j in 1:ns){
      curll[j]<-sum(dSplitNorm_cpp(y[,j],X[,j,],med[,j],theta[,j,],log=T))
      z[,j]<-qnorm(pSplitNorm_cpp(y[,j],X[,j,],med[,j],theta[,j,]))
    }
    canz<-z
    
  }
  
  # added by QB: use inital values
  if(length(iniallist)>0) { 
    len <- length(iniallist$range)
    
    range <- iniallist$range[len]
    rangez <- iniallist$range_res[len]
    ratioz <- iniallist$ratio_res[len]
    med <- iniallist$beta[len, , ]
    theta <- iniallist$theta[len, , , ]
    
    mn_med <- iniallist$mn_med
    mn_theta <- iniallist$mn_theta
    tau_med <- iniallist$tau_med
    tau_theta <- iniallist$tau_theta
    P <-iniallist$P
    Pz <- iniallist$Pz
    z <- iniallist$z
    canz <- iniallist$canz
    curll <-iniallist$curll
    
  # MHmed and MHtheta will put here.
    burn <- 100 # if initallist is given, then less burn
  }
  
  
  ################################################:
  ########     KEEP TRACK OF STUFF    ############:
  ################################################:
  keep.range<-keep.rangez<-keep.ratioz<-rep(0,iters)
  keep.med<-array(0,c(iters,p,ns))
  keep.theta<-array(0,c(iters,p,ns,L))
  Q1<-Q2<-Q3<-array(0,c(ntau,p,ns))
  dimnames(Q1)<-dimnames(Q2)<-dimnames(Q3)<-list(
    paste("tau=",round(tau,3),sep=""),
    paste("X",1:p,sep=""),
    paste("Site #",1:ns,sep=""))
  
  attmed<-accmed<-rep(0,p)
  atttheta<-acctheta<-matrix(0,p,L)
  
  
  start.time<-proc.time()
  for(i in 1:iters){
    ################################################:
    ###########     UPDATE MEDIAN    ###############:
    ################################################:
    
    a <-UPDATEMEDIAN(y,X,med,mn_med,tau_med,canz,curll,MHmed,theta,attmed,accmed,P,Pz,z)
    med <- a$med
    attme <- a$attmed
    accmed <- a$accmed
    canz <- a$canz
    curll <- a$curll
    z <- a$z
    rm(a)
    
    #hyperparameters:
    
    for(j in 1:p){
      SSS<-quad.form(P,med[j,]-mn_med[j])
      tau_med[j]<-rgamma(1,ns/2+a.var,SSS/2+b.var)
      
      P1<-rowSums(P)
      VVV<-tau_med[j]*sum(P1)+1/sd.mean^2
      MMM<-tau_med[j]*sum(P1*med[j,])
      mn_med[j]<-rnorm(1,MMM/VVV,1/sqrt(VVV))
    }
    
    
    
    ################################################:
    ###########     UPDATE THETA    ################:
    ################################################:
    b<- UPDATETHETA(y,X,med,mn_theta,tau_theta,canz,curll,MHtheta,theta,atttheta,acctheta,P,Pz,z);
    theta <- b$theta
    atttheta <- b$atttheta
    acctheta <- b$acctheta
    canz <- b$canz
    curll <- b$curll
    z <- b$z
    rm(b)
    #hyperparameters:
    for(j in 1:p){for(l in 1:L){
      P1<-rowSums(P)
      P2<-sum(P1)
      
      RRR<-theta[j,,l]-mn_theta[j,l]
      SSS<-quad.form(P,RRR)
      tau_theta[j,l]<-rgamma(1,ns/2+a.var,SSS/2+b.var)
      
      VVV<-tau_theta[j,l]*P2+1/sd.mean^2
      MMM<-tau_theta[j,l]*sum(P1*theta[j,,l])
      mn_theta[j,l]<-rnorm(1,MMM/VVV,1/sqrt(VVV))
    }}
    
    
    canrange<-exp(rnorm(1,log(range),0.05))
    canP<-solve(covexp_cpp(d,canrange,1))
    
    SSo<-SSn<-0
    for(j in 1:p){for(l in 1:(L+1)){
      if(l<L+1){
        yyy<-theta[j,,l]-mn_theta[j,l]
        ttt<-tau_theta[j,l]
        SSn<-SSn-0.5*ttt*quad.form(canP,yyy)
        SSo<-SSo-0.5*ttt*quad.form(P,yyy)
      }
      if(l==L+1){
        yyy<-med[j,]-mn_med[j]
        ttt<-tau_med[j] 
        SSn<-SSn-0.5*ttt*quad.form(canP,yyy)
        SSo<-SSo-0.5*ttt*quad.form(P,yyy)
      }
    }}
    
    R<-dnorm(log(canrange),mn_range,sd_range,log=T)-
      dnorm(log(range),mn_range,sd_range,log=T)+
      0.5*p*(L+1)*determinant(canP)$modulus-
      0.5*p*(L+1)*determinant(P)$modulus+
      SSn-SSo
    if(!is.na(exp(R))){if(runif(1)<exp(R)){
      range<-canrange;
      P<-canP
    }}
    
    ###########################################
    ######    RESIDUAL CORRELATION      #######:
    ###########################################
    
    curL<-0.5*nt*determinant(Pz)$modulus
    for(t in 1:nt){
      curL<-curL-0.5*quad.form(Pz,z[t,])
    }
    
    canratioz<-rnorm(1,ratioz,0.025)
    if(canratioz>0 & canratioz<1){
      canPz<-solve(covexp_cpp(d,rangez,canratioz))
      canL<-0.5*nt*determinant(canPz)$modulus
      for(t in 1:nt){
        canL<-canL-0.5*quad.form(canPz,z[t,])
      }
      R<-exp(canL-curL)
      if(!is.na(R)){if(runif(1)<R){
        ratioz<-canratioz;
        Pz<-canPz;
        curL<-canL
      }}
    }
    
    canrangez<-exp(rnorm(1,log(rangez),0.025))
    canPz<-solve(covexp_cpp(d,rangez,canratioz))
    canL<-0.5*nt*determinant(canPz)$modulus
    for(t in 1:nt){
      canL<-canL-0.5*quad.form(canPz,z[t,])
    }
    R<-canL-curL+
      dnorm(log(canrangez),mn_range,sd_range,log=T)-
      dnorm(log(rangez),mn_range,sd_range,log=T)
    R<-exp(R)
    if(!is.na(R)){if(runif(1)<R){
      rangez<-canrangez;
      Pz<-canPz;
      curL<-canL
    }}
    
    
    ###########################################
    ######  UPDATE MH CANDIDATE DIST    #######:
    ###########################################
    
    for(j in 1:p){if(i<burn/2 & attmed[j]>50){
      if(accmed[j]/attmed[j]<0.3){MHmed[j]<-MHmed[j]*0.8}
      if(accmed[j]/attmed[j]>0.6){MHmed[j]<-MHmed[j]*1.2}
      accmed[j]<-attmed[j]<-0
    }}
    for(j in 1:p){for(l in 1:L){
      if(i<burn/2 & atttheta[j,l]>50){
        if(acctheta[j,l]/atttheta[j,l]<0.3){MHtheta[j,l]<-MHtheta[j,l]*0.8}
        if(acctheta[j,l]/atttheta[j,l]>0.6){MHtheta[j,l]<-MHtheta[j,l]*1.2}
        acctheta[j,l]<-atttheta[j,l]<-0
      }
    }}
    MHmed[MHmed<sdy/1000]<-sdy/1000
    MHtheta[MHtheta<sdy/1000]<-sdy/1000
    
    
    
    ###########################################
    ######      KEEP TRACK OF OUTPUT    #######:
    ###########################################
    
    keep.range[i]<-range
    keep.rangez[i]<-rangez
    keep.ratioz[i]<-ratioz
    keep.theta[i,,,]<-theta
    keep.med[i,,]<-med
    
    
    if(i>burn){
      nnn<-iters-burn
      X0<-matrix(0,1,p)
      X0[,1]<-1
      
      for(k in 1:ns){
        qqq0<-qSplitNorm_cpp(tau,X0,med[,k],theta[,k,])
        for(j in 1:p){
          XXX<-matrix(0,1,p)
          XXX[,c(1,j)]<-1
          qqq<-qSplitNorm_cpp(tau,XXX,med[,k],theta[,k,])
          if(j>1){qqq<-qqq-qqq0}
          Q1[,j,k]<-Q1[,j,k]+qqq/nnn
          Q2[,j,k]<-Q2[,j,k]+qqq*qqq/nnn
          Q3[,j,k]<-Q3[,j,k]+(qqq>0)/nnn
        }
      }
    }
    
    ####  Display current iteration:
    if(i%%update==0){
      print(paste("Done with",i,"of",iters,"samples"))
      if(i==update){
        now<-proc.time()
        now<-as.numeric(now[3]-start.time[3])
        print("Approximate time left (min)")
        print((iters*now/update)/60)
        print("Approximate time left (hours)")
        print((iters*now/update)/(60*60))
      }
      if(i%%100==0){
	intermlist <- list(
    	range=keep.range,
    	range_res=keep.rangez,
    	ratio_res=keep.ratioz,
    	qfx.mn=Q1,
    	qfx.var=Q2-Q1^2,
    	qfx.prob.g.0=Q3,
    	theta=keep.theta,
    	beta=keep.med,
    	MHmed = MHmed, 
    	MHtheta = MHtheta,
    	mn_med = mn_med,
    	mn_theta = mn_theta,
    	tau_med = tau_med,
    	tau_theta = tau_theta,
    	P = P,
    	Pz = Pz,
    	z = z,
    	canz = canz,
    	curll =curll
    	
  	)
	saveRDS(intermlist,file ="intermlist.rds")
	rm(intermlist)
	}
    }
    
  }
  
  
  list(
    range=keep.range,
    range_res=keep.rangez,
    ratio_res=keep.ratioz,
    qfx.mn=Q1,
    qfx.var=Q2-Q1^2,
    qfx.prob.g.0=Q3,
    theta=keep.theta,
    beta=keep.med,
    MHmed = MHmed, 
    MHtheta = MHtheta,
    mn_med = mn_med,
    mn_theta = mn_theta,
    tau_med = tau_med,
    tau_theta = tau_theta,
    P = P,
    Pz = Pz,
    z = z,
    canz = canz,
    curll =curll
  )}


############################################################################:
#############     OTHER FUNCTIONS USED IN THE LIKELIHOOD        ############:
############################################################################:


