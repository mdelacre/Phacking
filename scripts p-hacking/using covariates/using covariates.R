library(MASS)

for (package in "MASS") {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}

# Generate DV and covariates

correlDV <- function(n, mu, sds,nb_covar,corr_cov=.5,corr_dvcov=0) {
  mu=rep(mu,(nb_covar+1))
  Sigma <- matrix(corr_cov, nrow = (nb_covar+1), ncol = (nb_covar+1))
  Sigma[,1]=corr_dvcov # corrélation entre la VD et chaque covariable
  Sigma[1,]=corr_dvcov # corrélation entre la VD et chaque covariable
  diag(Sigma) <- 1
  rawvars <- mvrnorm(n=n, mu=mu, Sigma=Sigma)*sds
  return(rawvars)
}

covariates_moresig <- function(n=100,nb_covar=3,corr_cov=.5,corr_dvcov=0,    # arguments relatifs aux DV corrélées
                               mu=0,sds=2,                           # arguments relatifs aux éch générés
                               alpha=.05, graph=TRUE, nSims=100000){   # arguments relatifs au p-hacking
  
  out <- matrix(0, nrow=nSims, ncol=2)
  
  for(i in 1:nSims) { 
    g1<-correlDV(n=n,mu=mu,sds=sds,nb_covar=nb_covar, corr_cov=corr_cov,corr_dvcov=corr_dvcov)
    g2<-correlDV(n=n,mu=mu,sds=sds,nb_covar=nb_covar, corr_cov=corr_cov,corr_dvcov=corr_dvcov)
    quant_var <- rbind(g1,g2)
    fact <- rep(1:2,each=n)
    
    # réaliser diverses ancova en condidérant une seule covariable à chaque fois
    test<-1:nb_covar
    for (j in 1:nb_covar) {test[j]=summary(lm(quant_var[,1] ~ quant_var[,j+1] + fact))$coefficients[3,4]}
    out[i,1]=test[1]
    out[i,2]=min(test)
  }
  initialsig_p=out[,1][out[,1]<.05]
  hackedsig_p=out[,2][out[,2]<.05]
  
  freq.sigfirst<-length(initialsig_p)/nSims
  freq.sigafter<-length(hackedsig_p)/nSims
  
  par(mfrow=c(2,1))
  A=hist(initialsig_p,plot=F,breaks=5)
  B=hist(hackedsig_p,plot=F,breaks=5)  
  A$counts=A$counts/sum(A$counts)
  B$counts=B$counts/sum(B$counts)
  plot(A,freq=T,ylab="Relative Frequency",xlim=c(0,alpha),main=paste("p curve when one automatically selects the first covariate ","\n","frequency of p-values under",alpha,"=",round(freq.sigfirst,3)),col="pink")
  plot(B,freq=T,ylab="Relative Frequency",xlim=c(0,alpha),main=paste("p curve","\n","selection among ",nb_covar,"covariates with r =",corr_cov,"\n","frequency of p-values under",alpha,"=",round(freq.sigafter,3)),col="lightblue")

}

covariates_moresig(n=100,nb_covar=3,corr_cov=.5,corr_dvcov=.3,mu=0,sds=2,alpha=.05, graph=TRUE, nSims=10000)
covariates_moresig(n=100,nb_covar=3,corr_cov=0,corr_dvcov=.3,mu=0,sds=2,alpha=.05, graph=TRUE, nSims=10000)
covariates_moresig(n=100,nb_covar=5,corr_cov=.1,corr_dvcov=.5,mu=0,sds=2,alpha=.05, graph=TRUE, nSims=10000)

#####  2)   select a covariate only if she gives significant resuls, otherwise looking for another one ######
###########################################################################################################

covariates_firstSig <- function(n=100,nb_covar=3,corr_cov=.5,corr_dvcov=0,    # arguments relatifs aux DV corrélées
                                mu=0,sds=2,                           # arguments relatifs aux éch générés
                                alpha=.05, graph=TRUE, nSims=100000){   # arguments relatifs au p-hacking
  
  out <- matrix(0, nrow=nSims, ncol=3)
  
  for(i in 1:nSims) {
    g1<-correlDV(n=n,mu=mu,sds=sds,nb_covar=nb_covar, corr_cov=corr_cov,corr_dvcov=corr_dvcov)
    g2<-correlDV(n=n,mu=mu,sds=sds,nb_covar=nb_covar, corr_cov=corr_cov,corr_dvcov=corr_dvcov)
    quant_var <- rbind(g1,g2)
    fact <- rep(1:2,each=n)
    
    # réaliser un t.test sur chaque DV
    
    initial.p=summary(lm(quant_var[,1] ~ quant_var[,2] + fact))$coefficients[3,4]
    current_p=initial.p
    explored_dv <- 1
    while(current_p > alpha & explored_dv < nb_covar) {    # add subjects as long as p is not significant
      explored_dv <- explored_dv + 1
      current_p=summary(lm(quant_var[,1] ~ quant_var[,(explored_dv+1)] + fact))$coefficients[3,4]
    }
    out[i,1]=initial.p
    out[i,2]=current_p
    out[i,3]=explored_dv
  }
  initialsig_p=out[,1][out[,1]< alpha]
  hackedsig_p=out[,2][out[,2]< alpha]
  count=out[,3]#[out[,1]<.05] 
  
  freq.sigfirst<-length(initialsig_p)/nSims
  freq.sigafter<-length(hackedsig_p)/nSims
  
  par(mfrow=c(2,1),mar = rep(4, 4))
  A=hist(initialsig_p,plot=F,breaks=5)   
  B=hist(hackedsig_p,plot=F,breaks=5)   
  A$counts=A$counts/sum(A$counts)
  B$counts=B$counts/sum(B$counts)
  plot(A,freq=T,xlim=c(0,alpha),ylab="Relative Frequency",xlab="p-values",main=paste("p curve when one automatically selects the first DV ","\n","frequency of p-values under",alpha,"=",round(freq.sigfirst,2)),col="pink")
  plot(B,freq=T,xlim=c(0,alpha),ylab="Relative Frequency",xlab= "p-values",main=paste("p curve when selection among ",nb_DV,"DV with r =",corr,"\n","average explored dv=",round(mean(count),2),"\n","frequency of p-values under",alpha,"=",round(freq.sigafter,2)),col="lightblue")
}
 
covariates_firstsig(n=100,nb_covar=3,corr_cov=.5,corr_dvcov=.3,mu=0,sds=2,alpha=.05, graph=TRUE, nSims=1000)
covariates_moresig(n=100,nb_covar=3,corr_cov=0,corr_dvcov=.3,mu=0,sds=2,alpha=.05, graph=TRUE, nSims=10000)
covariates_moresig(n=100,nb_covar=5,corr_cov=.1,corr_dvcov=.5,mu=0,sds=2,alpha=.05, graph=TRUE, nSims=10000)


