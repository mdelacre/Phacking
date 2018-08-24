## K DV, select it only if all are significant
# compute the type 1 error rate

# see Schimmack, 2012 and Ulrich et al. 2014

## Require the script "Generate x DV that are correlated.R"

correlDV <- function(n,sds,nb_DV,corr) {
  
  mu=rep(0,nb_DV)
  Sigma <- matrix(corr, nrow = nb_DV, ncol = nb_DV)
  diag(Sigma) <- 1
  rawvars <- mvrnorm(n=n, mu=mu, Sigma=Sigma)*sds
  return(rawvars)
}


########          1)   select the more significant result among many DV             ########
############################################################################################
 
conservative_res <- function(n=20, correlation=0.5, nb_DV=5,sds=2,             # arguments relatifs aux éch générés
                             expected_Nsig=5, alpha=.05, graph=TRUE, nSims=100000){   # arguments relatifs au p-hacking
  
  out <- list()

  for(i in 1:nSims) {
    g1<-correlDV(n=n,corr=corr,nb_DV=nb_DV,sds=sds)
    g2<-correlDV(n=n,corr=corr,nb_DV=nb_DV,sds=sds)
    DVs <- rbind(g1,g2)
    fact <- rep(1:2,each=n)
    
    # réaliser un t.test sur chaque DV
    test<-1:nb_DV
    for (j in 1:nb_DV) {test[j]=t.test(x=DVs[,j][fact==1],y=DVs[,j][fact==2])$p.value}
    p_values=NULL
    if (sum(test<.05)>=expected_Nsig) {results=test[test<.05]} else {results=NULL}
    out[[i]]=c(p_values,results) 
  }
  
  p_values=unlist(out)
  freq.sig<-length(p_values)/(nSims*nb_DV)

  A=hist(p_values,plot=F,breaks=5)
  A$counts=A$counts/sum(A$counts)
  plot(A,freq=T,ylab="Relative Frequency",xlim=c(0,.05),main=paste("p curve when a min of",expected_Nsig,"significant results are required in order to report sig results","\n","frequency of p-values under",alpha,"=",round(freq.sig,3)),col="pink")

}

conservative_res(n=20, corr=0.5, nb_DV=2,sds=2,expected_Nsig=2, alpha=.05, graph=TRUE, nSims=100000)
