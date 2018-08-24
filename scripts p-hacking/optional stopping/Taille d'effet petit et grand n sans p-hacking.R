for (package in c("psych","parallel")) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}

phacking <- function(n=20, bign=100,m1=0, m2=0, sd_pooled=1, alpha=.05, alternative="two.sided", graph=TRUE, nSims=1000) {
  outmat <- matrix(NA, nrow=nSims, ncol=4)
  
  real_effectsize <- (m1-m2)/sd_pooled
  
  for(i in 1:nSims) {
    g1 <- rnorm(n, m1, sd_pooled)
    g2 <- rnorm(n, m2, sd_pooled)
    dv1 <- c(g1,g2)
    
    g3 <- rnorm(bign, m1, sd_pooled)
    g4 <- rnorm(bign, m2, sd_pooled)
    dv2 <- c(g3,g4)
    
    fact <-c(rep(1,length(g1)),rep(0,length(g2)))
    
    pooled_sd1 <- sqrt(((length(g1)-1)*sd(g1)^2+(length(g2)-1)*sd(g2)^2)/(length(dv1)-2))
    d1 <- (mean(g1)-mean(g2))/pooled_sd1 # cohen's d
    p1 <-t.test(g1, g2, alternative=alternative)$p.value  
    pooled_sd2 <- sqrt(((length(g3)-1)*sd(g3)^2+(length(g4)-1)*sd(g4)^2)/(length(dv2)-2))
    d2 <- (mean(g3)-mean(g4))/pooled_sd2 # cohen's d
    p2 <-t.test(g3, g4, alternative=alternative)$p.value 
    outmat[i,] <- c(d1,p1,d2,p2)
  }
  
  outmat <- data.frame(outmat)
  colnames(outmat) <- c("d1","p1","d2","p2")

    par(mfrow=c(2,1), las=1, font.main=1)
    C=hist(outmat$d1[outmat$p1<alpha],plot=F)
    C$counts=C$counts/sum(C$counts)
    plot(C,freq=T,ylab="Relative Frequency",xlab="cohen's d", main="Effect size with small N",col="lightgreen",xlim=c(-3,1))
    
    D=hist(outmat$d2[outmat$p2<alpha],plot=F)
    D$counts=D$counts/sum(D$counts)
    plot(D,freq=T,ylab="Relative Frequency",xlab="cohen's d", main="Effect size with big N",col="lightgreen",xlim=c(-3,1))
    

  return(outmat)
}

res1 <- phacking(alpha=.1, alternative="two.sided", graph=TRUE, nSims=100000)
res2 <- phacking(alpha=.05, alternative="two.sided", graph=TRUE, nSims=100000)
res3 <- phacking(alpha=.01, alternative="two.sided", graph=TRUE, nSims=100000)

res4 <- phacking(alpha=.1, m1=0, m2=1, alternative="two.sided", graph=TRUE, nSims=100000)
res5 <- phacking(alpha=.05, m1=0, m2=1, alternative="two.sided", graph=TRUE, nSims=100000)
res6 <- phacking(alpha=.01, m1=0, m2=1, alternative="two.sided", graph=TRUE, nSims=100000)

