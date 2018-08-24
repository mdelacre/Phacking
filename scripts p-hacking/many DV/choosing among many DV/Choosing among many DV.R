## Require the script "Generate x DV that are correlated.R"

# Correlation pour covariates
# rem.: si je veux générer de l'hétéroscédasticité, il suffit de multiplier les colonnes par le sd désiré.
# car par défaut, la sd de tous les groupes vaut 1.
library(MASS)
correlDV <- function(n,sds,nb_DV,corr) {
  
  mu=rep(0,nb_DV)
  Sigma <- matrix(corr, nrow = nb_DV, ncol = nb_DV)
  diag(Sigma) <- 1
  rawvars <- mvrnorm(n=n, mu=mu, Sigma=Sigma)*sds
  return(rawvars)
}

#######          1)   select the more significant result among many DV             ########
############################################################################################
 
phacking_moreSig <- function(n=20, nb_DV=5,corr=0.5, sds=2,          # arguments relatifs aux éch générés
                             alpha=.05, graph=TRUE, nSims=100000){   # arguments relatifs au p-hacking

  out <- matrix(0, nrow=nSims, ncol=2)

  for(i in 1:nSims) { 
    g1<-correlDV(n=n,corr=corr,nb_DV=nb_DV,sds=sds)
    g2<-correlDV(n=n,corr=corr,nb_DV=nb_DV,sds=sds)
    DVs <- rbind(g1,g2)
    fact <- rep(1:2,each=n)

    # réaliser un t.test sur chaque DV 
    test<-1:nb_DV
    for (j in 1:nb_DV) {test[j]=t.test(x=DVs[,j][fact==1],y=DVs[,j][fact==2])$p.value}

    out[i,1]=test[1]
    out[i,2]=min(test)
  }
  
  initialsig_p=out[,1][out[,1]< alpha]
  hackedsig_p=out[,2][out[,2]< alpha]
  
  freq.sigfirst<-length(initialsig_p)/nSims
  freq.sigafter<-length(hackedsig_p)/nSims

  if (graph==TRUE){
  par(mfrow=c(2,1))
  A=hist(initialsig_p,plot=F,breaks=5)
  B=hist(hackedsig_p,plot=F,breaks=5)   
  A$counts=A$counts/sum(A$counts)
  B$counts=B$counts/sum(B$counts)
  plot(A,freq=T,ylab="Relative Frequency",xlim=c(0,alpha),main=paste("p curve when one automatically selects the first DV ","\n","frequency of p-values under",alpha,"=",round(freq.sigfirst,3)),col="pink")
  plot(B,freq=T,ylab="Relative Frequency",,xlim=c(0,alpha),main=paste("p curve","\n","selection among ",nb_DV,"DV with r =",corr,"\n","frequency of p-values under",alpha,"=",round(freq.sigafter,3)),col="lightblue")

      }
} 

phacking_moreSig(n=20,sd=1,corr=.5, nb_DV=5,alpha=.05,nSims=1000)

library(aprof)


#####  2)   select a DV only if she is significant, otherwise looking for another one ######
############################################################################################

phacking_firstSig <- function(n=20, corr=.5, nb_DV=5, sds=2,     
                             alpha=.05, graph=TRUE, nSims=100000){   # arguments relatifs au p-hacking
    
  out <- matrix(0, nrow=nSims, ncol=3)

  for(i in 1:nSims) {
    g1<-correlDV(n=n,corr=corr,nb_DV=nb_DV,sds=sds)
    g2<-correlDV(n=n,corr=corr,nb_DV=nb_DV,sds=sds)
    DVs <- rbind(g1,g2)
    
    fact <- rep(1:2,each=n)
    
    # réaliser un t.test sur chaque DV

    initial.p=t.test(x=DVs[,1][fact==1],y=DVs[,1][fact==2])$p.value
    current_p=initial.p
    explored_dv <- 1
      while(current_p > alpha & explored_dv < nb_DV) {    # add subjects as long as p is not significant
      explored_dv <- explored_dv + 1
      current_p=t.test(x=DVs[,explored_dv][fact==1],y=DVs[,explored_dv][fact==2])$p.value
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

  par(mfrow=c(2,1))
    A=hist(initialsig_p,plot=F,breaks=5)  
    B=hist(hackedsig_p,plot=F,breaks=5)   
    A$counts=A$counts/sum(A$counts)
    B$counts=B$counts/sum(B$counts)
    plot(A,freq=T,xlim=c(0,alpha),ylab="Relative Frequency",xlab="p-values",main=paste("p curve when one automatically selects the first DV ","\n","frequency of p-values under",alpha,"=",round(freq.sigfirst,2)),col="pink")
    plot(B,freq=T,xlim=c(0,alpha),ylab="Relative Frequency",xlab= "p-values",main=paste("p curve when selection among ",nb_DV,"DV with r =",corr,"\n","average explored dv=",round(mean(count),2),"\n","frequency of p-values under",alpha,"=",round(freq.sigafter,2)),col="lightblue")
}
 
phacking_moreSig(n=20,corr=0.2,nb_DV=2,sds=2,nSims=100000)
phacking_moreSig(n=20,corr=0.2,nb_DV=4,sds=2,nSims=100000)
phacking_moreSig(n=20,corr=0.2,nb_DV=8,sds=2,nSims=100000)
phacking_moreSig(n=20,corr=0.2,nb_DV=32,sds=2,nSims=100000)
phacking_moreSig(n=20,corr=0.5,nb_DV=2,sds=2,nSims=100000)
phacking_moreSig(n=20,corr=0.5,nb_DV=4,sds=2,nSims=100000)
phacking_moreSig(n=20,corr=0.5,nb_DV=8,sds=2,nSims=100000)
phacking_moreSig(n=20,corr=0.5,nb_DV=16,sds=2,nSims=100000)
phacking_moreSig(n=20,corr=0.8,nb_DV=2,sds=2,nSims=100000)
phacking_moreSig(n=20,corr=0.8,nb_DV=4,sds=2,nSims=100000)
phacking_moreSig(n=20,corr=0.8,nb_DV=8,sds=2,nSims=100000)
phacking_moreSig(n=20,corr=0.8,nb_DV=16,sds=2,nSims=100000)

phacking_firstSig(n=20,corr=0.2,nb_DV=2,sds=2,nSims=100)
phacking_firstSig(n=20,corr=0.2,nb_DV=4,sds=2,nSims=100000)
phacking_firstSig(n=20,corr=0.2,nb_DV=8,sds=2,nSims=100000)
phacking_firstSig(n=20,corr=0.2,nb_DV=16,sds=2,nSims=100000)
phacking_firstSig(n=20,corr=0.5,nb_DV=2,sds=2,nSims=100000)
phacking_firstSig(n=20,corr=0.5,nb_DV=4,sds=2,nSims=100000)
phacking_firstSig(n=20,corr=0.5,nb_DV=8,sds=2,nSims=100000)
phacking_firstSig(n=20,corr=0.5,nb_DV=16,sds=2,nSims=100000)
phacking_firstSig(n=20,corr=0.8,nb_DV=2,sds=2,nSims=100000)
phacking_firstSig(n=20,corr=0.8,nb_DV=4,sds=2,nSims=100000)
phacking_firstSig(n=20,corr=0.8,nb_DV=8,sds=2,nSims=100000)
phacking_firstSig(n=20,corr=0.8,nb_DV=16,sds=2,nSims=100000)

