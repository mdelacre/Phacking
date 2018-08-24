## Require the script "Generate x DV that are correlated.R"

correlDV <- function(n,sds,nb_DV,corr) {
  
  mu=rep(0,nb_DV)
  Sigma <- matrix(corr, nrow = nb_DV, ncol = nb_DV)
  diag(Sigma) <- 1
  rawvars <- mvrnorm(n=n, mu=mu, Sigma=Sigma)*sds
  return(rawvars)
}

##### 1) probability that at least one DV is significant, among 2 DV and its average ######
###########################################################################################
# @Simmons et al.

phacking_freqSig <- function(n=20, corr=.5,sds=2,             # arguments relatifs aux éch générés
                             alpha=.05, graph=TRUE, nSims=100000){   # arguments relatifs au p-hacking
  
  out <- matrix(0, nrow=nSims, ncol=5)
  colnames(out)=c("DV1","DV2","Combi","at least one sig","Lowest_p")
  for(i in 1:nSims) {
    g1<-correlDV(n=n,corr=corr,nb_DV=2,sds=sds)
    g2<-correlDV(n=n,corr=corr,nb_DV=2,sds=sds)
    A <- rbind(g1,g2)
    B<- (A[,1]+A[,2])/2
    DVs<-cbind(A,B)
    fact <- rep(1:2,each=n)
    
    # réaliser un t.test sur chaque DV
    test<-1:3
    for (j in 1:3) {test[j]=t.test(x=DVs[,j][fact==1],y=DVs[,j][fact==2])$p.value}
    out[i,1:3]=test
    out[i,4]=any(test<alpha)
    out[i,5]= which(test==min(test))
  }
  
  alpha=length(out[,4][out[,4]==1])/nSims
  
  allres<-out[,5]
  sigres<-subset(out[,5], out[,4]==1)
  
  par(mfrow=c(2,1))
  barplot(table(allres),main="Among all results, how often does the smallest value fall in each categories?",names.arg=c("A","B","C"),col="pink")
  barplot(table(sigres),main="Among significant results, how often does the smallest value fall in each categories?",names.arg=c("A","B","C"),col="pink")
  
  return(alpha)
}

##### 2) probability that at least one DV is significant, among 2 DV and its weighted average ######
####################################################################################################
# same script than 1) except that scores are weighted by "-1" if t stat is negative, by "+1" if t stat is positive
# in order to adjust the influence of each variable to maximize the mean difference between groups
# as suggested by Ulrich and Miller (2014)

phackingweighted_freqSig <- function(n=20, corr=.5, nb_DV=2,sds=2,              # arguments relatifs aux éch générés
                             alpha=.05, graph=TRUE, nSims=100000){   # arguments relatifs au p-hacking

  out <- matrix(0, nrow=nSims, ncol=(nb_DV+3))
  colnames(out)=c("DV1","DV2","Combi","at least one sig","Lowest_p")

    for(i in 1:nSims) {

    DVs <- matrix(0, nrow=2*n, ncol=nb_DV+1)
    names=expand.grid(DV="DV",numb=1:nb_DV)
    colnames(DVs)=c(paste0(names$DV,names$numb),"Compo")

    DVs[(1:n),1:nb_DV] <- correlDV(n=n,corr=corr,nb_DV=nb_DV,sds=sds)
    DVs[(n+1):(2*n),1:nb_DV] <- correlDV(n=n,corr=corr,nb_DV=nb_DV,sds=sds)
    fact <- rep(1:2,each=n)
    # calculer une statistique t pour chaque vd et identifier le signe de la statistique (en vue de la pondération)
    tstat<-1:(nb_DV)
    w<-1:(nb_DV)
    for (j in 1:(nb_DV)) { 
      tstat[j]=t.test(x=DVs[,j][fact==1],y=DVs[,j][fact==2])$statistic
      if (tstat[j]<0){w[j]=-1}
      else {w[j]=1}
      }
    weighted_ddb<-DVs[,1:nb_DV]*matrix(rep(w,each=2*n),2*n,nb_DV)
    DVs[,nb_DV+1]<-apply(weighted_ddb[,1:nb_DV],1,sum)/nb_DV
    
    # réaliser un t.test sur chaque DV (+ la combi) et extraire la p-valeur
    test<-1:(nb_DV+1)
    for (j in 1:(nb_DV+1)) {test[j]=t.test(x=DVs[,j][fact==1],y=DVs[,j][fact==2])$p.value}
    out[i,1:(nb_DV+1)]=test
    out[i,(nb_DV+2)]=any(test<alpha)
    out[i,(nb_DV+3)]= which(test==min(test))  
   }

  alpha=length(out[,(nb_DV+2)][out[,(nb_DV+2)]==1])/nSims
  
  allres<-out[,5]
  sigres<-subset(out[,5], out[,4]==1)
   
  par(mfrow=c(2,1))
  barplot(table(allres),main="Among all results, how often does the smallest value fall in each categories?",names.arg=c("A","B","C"),col="pink")
  barplot(table(sigres),main="Among significant results, how often does the smallest value fall in each categories?",names.arg=c("A","B","C"),col="pink")
  
  return(alpha)
  
}

#####                       combine the m more significant p-values                      ######
###############################################################################################

# Ulrich et al.

phacking_combimoreSig <- function(n=20, corr=.5, nb_DV=4,M=2,sds=2,                # arguments relatifs aux éch générés
                             alpha=.05, graph=TRUE, nSims=100000){    # arguments relatifs au p-hacking
  
  combi_p=rep(0,nSims)
  smallest_p=rep(0,nSims)
  firstDV_p=rep(0,nSims)

  
  for(i in 1:nSims) {
    DVs <- matrix(0, nrow=2*n, ncol=nb_DV)
    names=expand.grid(DV="DV",numb=1:nb_DV)
    colnames(DVs)=paste0(names$DV,names$numb)
    DVs[(1:n),1:nb_DV] <- correlDV(n=n,corr=corr,nb_DV=nb_DV,sds=sds)
    DVs[(n+1):(2*n),1:nb_DV] <- correlDV(n=n,corr=corr,nb_DV=nb_DV,sds=sds)
    fact <- rep(1:2,each=n)

    # calculer une statistique t pour chaque vd et identifier le signe de la statistique (en vue de la pondération)
    pvalue<-1:(nb_DV)
    tstat<-1:(nb_DV)
    w<-1:(nb_DV)
    for (j in 1:(nb_DV)) { i
      res=t.test(x=DVs[,j][fact==1],y=DVs[,j][fact==2])
      pvalue[j]=res$p.value
      tstat[j]=res$statistic
      if (tstat[j]<0){w[j]=-1}
      else {w[j]=1} }
    
    # sélectionner les m DV associées à la plus petite p-valeur
    selectDV <- rep(0,M)
    for(m in 1:M){
      selectDV[m]=which(pvalue==sort(pvalue,partial=m)[m])} 
    
    weighted_ddb<-DVs[,selectDV]*matrix(rep(w[selectDV],each=2*n),2*n,M)
    combination_DV<-apply(weighted_ddb,1,sum)/M

    firstDV_p[i] <- pvalue[1]
    smallest_p[i] <- min(pvalue)
    combi_p[i]<-t.test(x=combination_DV[fact==1],y=combination_DV[fact==2])$p.value

  } 

  freq.sigfirst<-length(firstDV_p[firstDV_p<.05])/nSims
  freq.sigsmallest<-length(smallest_p[smallest_p<.05])/nSims
  freq.combi<-length(combi_p[combi_p<.05])/nSims
  
  par(mfrow=c(3,1))
  # histograms 
  # si on prend systématiquement la première p-valeur
  A=hist(firstDV_p,plot=F,breaks=5)  
  A$counts=A$counts/sum(A$counts)  
  plot(A,freq=T,ylab="Relative Frequency",xlab="p-values",main=paste("p curve","\n","selection of the first DV among",nb_DV,"DVs with r =",correlation,"\n","frequency of p-values under",alpha,"=",round(freq.sigfirst,2)),xlim=c(0,1),col="pink")
  # si on prend systématiquement la plus petite p-valeur
  B=hist(smallest_p,plot=F,breaks=5)   
  B$counts=B$counts/sum(B$counts)
  plot(B,freq=T,ylab="Relative Frequency",xlab="p-values",main=paste("p curve","\n","selection of the DV with the smallest p-value among",nb_DV,"DVs with r =",correlation,"\n","frequency of p-values under",alpha,"=",round(freq.sigsmallest,2)),xlim=c(0,1),col="lightblue")
  # si on prend systématiquement la variable "combi"
  C=hist(combi_p,plot=F,breaks=5)
  C$counts=C$counts/sum(C$counts)
  plot(C,freq=T,ylab="Relative Frequency",xlab="p-values",main=paste("p curve","\n","weighted combination of",M,"DVs among",nb_DV,"DVs with r =",correlation,"\n","frequency of p-values under",alpha,"=",round(freq.combi,2)),xlim=c(0,1),col="lightyellow")
  
  }


phacking_combimoreSig(n=20,corr=0,sds=2,nSims=100000)

phacking_freqSig(n=20,correlation=0,distName=rep("normal",2),sds=2,nSims=100)
phacking_freqSig(n=20,correlation=0.1,distName=rep("normal",2),sds=2)
phacking_freqSig(n=20,correlation=0.2,distName=rep("normal",2),sds=2)
phacking_freqSig(n=20,correlation=0.3,distName=rep("normal",2),sds=2)
phacking_freqSig(n=20,correlation=0.4,distName=rep("normal",2),sds=2)
phacking_freqSig(n=20,correlation=0.5,distName=rep("normal",2),sds=2)
phacking_freqSig(n=20,correlation=0.6,distName=rep("normal",2),sds=2)
phacking_freqSig(n=20,correlation=0.7,distName=rep("normal",2),sds=2)
phacking_freqSig(n=20,correlation=0.8,distName=rep("normal",2),sds=2)
phacking_freqSig(n=20,correlation=0.9,distName=rep("normal",2),sds=2)
# à mesure que la corrélation entre A et B augmente, la combinaison des deux perd en intérêt
# pour p-hacker, toujours plus "productif" d'ajouter une troisième DV que d'envisager la combinaison de deux
# (mais ceci dit, combiner deux VD déjà présentes, ça ne coûte rien)

phackingweighted_freqSig(n=20,correlation=0,distName=rep("normal",2),sds=2)
phackingweighted_freqSig(n=20,correlation=0.1,distName=rep("normal",2),sds=2)
phackingweighted_freqSig(n=20,correlation=0.2,distName=rep("normal",2),sds=2)
phackingweighted_freqSig(n=20,correlation=0.3,distName=rep("normal",2),sds=2)
phackingweighted_freqSig(n=20,correlation=0.4,distName=rep("normal",2),sds=2)
phackingweighted_freqSig(n=20,correlation=0.5,distName=rep("normal",2),sds=2)
phackingweighted_freqSig(n=20,correlation=0.6,distName=rep("normal",2),sds=2)
phackingweighted_freqSig(n=20,correlation=0.7,distName=rep("normal",2),sds=2)
phackingweighted_freqSig(n=20,correlation=0.8,distName=rep("normal",2),sds=2)
phackingweighted_freqSig(n=20,correlation=0.9,distName=rep("normal",2),sds=2)
# à mesure que la corrélation entre A et B augmente, l'intérêt de la combinaison diminue un peu (même si reste nettement mieux qu'en cas de moyennes non pondérées, cf. fonction au dessus)
# pour p-hacker, envisager la combi de plusieurs DV ça cartonne, surtout quand elles ne sont pas corrélées
# (mais ceci dit, d'un point de vue crédibilité scientifique, c'est un peu chelou)
# de 1 ce sont les mm personnse qui répondent aux questions dc presque sûr qu'il y aura au moins un peu de lien
# en plus dur de justifier qu'on a combiné deux variables indépendantes en une seule
 
phacking_combimoreSig(n=20,correlation=0,distName=rep("normal",2),sds=2,nSims=100)
phacking_combimoreSig(n=20,correlation=0.1,distName=rep("normal",2),sds=2)
phacking_combimoreSig(n=20,correlation=0.2,distName=rep("normal",2),sds=2)
phacking_combimoreSig(n=20,correlation=0.3,distName=rep("normal",2),sds=2)
phacking_combimoreSig(n=20,correlation=0.4,distName=rep("normal",2),sds=2)
phacking_combimoreSig(n=20,correlation=0.5,distName=rep("normal",2),sds=2)
phacking_combimoreSig(n=20,correlation=0.6,distName=rep("normal",2),sds=2)
phacking_combimoreSig(n=20,correlation=0.7,distName=rep("normal",2),sds=2)
phacking_combimoreSig(n=20,correlation=0.8,distName=rep("normal",2),sds=2)
phacking_combimoreSig(n=20,correlation=0.9,distName=rep("normal",2),sds=2)
 

