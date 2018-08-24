############################################################################################
####                                                                                    ####
####                        Generate x DV that are correlated                           ####
####                                                                                    ####
############################################################################################

for (package in c("psych","parallel","smoothmest","fGarch")) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}

# Function to generate samples
get_sample     <- function(distName,n,              # n = sample size; distName = distribution underlying the data
                           lambda=2,                  # arguments for double exponential distribution
                           sds=2,                  # arguments for normal, normal skewed or double exponential distribution
                           bound=3.465,              # arguments for unif distribution
                           df=2)                      # argument for chi square distribution
{
  if (distName=="normal"){ out <- rnorm(n, mean=0, sd=sds)}
  else if (distName=="doublex"){out <- rdoublex(n, mu=0, lambda=lambda/sqrt(2))} # transformation in order that lambda = sd
  else if (distName=="skewpos"){out <- rsnorm(n, mean=0, sd=sds,xi=10)}
  else if (distName=="skewneg"){out <- rsnorm(n, mean=0, sd=sds,xi=-10)}
  else if (distName=="unif"){out <- runif(n, min=-bound, max=bound)}
  else if (distName=="chi2"){out <- rchisq(n, df=df)-2} # so sd=2, mean=0
  return (out)
}

correlDV <- function(n, correlation=.5,nb_DV=4,
                     distName,sds,bound,df) {
  # generate nb_DV variables, stored in a matrix
  X=matrix(0,sum(n),nb_DV)  
  for(v in 1:nb_DV){
    g1 <-get_sample(distName[1],n[1],lambda[1],sds[1],bound[1],df[1])
    g2 <-get_sample(distName[2],n[2],lambda[2],sds[2],bound[2],df[2])
    X[,v] <- c(g1,g2)
  }
  # generate the expected correlation matrix    
  C <- matrix(correlation, nrow = nb_DV, ncol = nb_DV)
  diag(C) <- 1 
  L <- chol(C) # Choleski factorization
  # Return a lower triangular matrix Z such as C=t(Z)%*%Z   
  # induce correlation, without changing X1 (as long as Z is a lower triangular matrix)
  DV <- X %*% L
  return(DV)
}

