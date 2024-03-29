library(maxLik)
###################################################################################################
#This file contains source functions of GBetareg                                                  #
#Given one CpG site, the basic model is: Z ~ beta0 + beta1*phenotype + beta2 * confounder1 + .... #
###################################################################################################
#input:
# Z: DNAm levels: n by L matrix; n is sample size; L is block size
# X: covariates: n by p matrix; n is sample size; p is number of covariates; the first column must be phenotype
#output
# an maxLik object consisting of maximum likelihood, estimate, etc...
get.mle.reg <-function(Z,X){
  X=as.matrix(X)
  Z=as.matrix(Z)
  p=ncol(X)-1
  n=nrow(Z)
  L=ncol(Z)
  log.Z=log(Z)
  log.1mZ=log(1-Z)
  log.denom=log(1+rowSums(Z/(1-Z)))
  starter.beta=mean((colMeans(Z)*(1-colMeans(Z))/apply(Z,MARGIN=2,var)-1)*(1-colMeans(Z)))
  if (p>0){
    starter.alpha=as.matrix(lm((log.Z-log.1mZ) ~ X[,2:(p+1)])[[1]])
    temp=rep(NA,length(p))
    for (i in 2:(p+1)){
      temp[i-1]=mean(starter.alpha[i,])
    }
    starter.alpha=c(starter.alpha[1,1:L],temp)
    starter=c(starter.alpha,starter.beta)
  }else{ 
    starter=c(as.vector(lm((log.Z-log.1mZ) ~ 1)[[1]]),starter.beta)
  }
  
  llk <- function(theta){
    coef.intercept=theta[1:L]
    if (p>0){
      coef.covariate=theta[(L+1):(L+p)]
      coef.beta=rbind(coef.intercept,matrix(coef.covariate,nrow=p,ncol=L))
      beta=theta[p+L+1]
      alpha=exp(X%*%coef.beta)*beta
      if (L!=1){
        nom1=sum(lgamma(rowSums(alpha)+ beta))
        nom2=sum((alpha-1)*log.Z)-sum((alpha+1)*log.1mZ)
        denom1=sum(lgamma(alpha))+n*lgamma(beta)
        denom2=sum((rowSums(alpha)+beta)*log.denom)
        return (nom1+nom2-denom1-denom2)
      }else{
        return (sum(dbeta(Z,alpha,rep(beta,length(alpha)),log=TRUE)))  
      }
    }else{
      coef.beta=matrix(coef.intercept,nrow=1)
      beta=theta[p+L+1]
      alpha.unique=exp(coef.beta)*beta
      alpha=matrix(exp(coef.beta)*beta,nrow=n,ncol=L,byrow=TRUE)
      if (L!=1){
        nom1=n*lgamma(sum(alpha.unique)+beta)
        nom2=sum((alpha-1)*log.Z)-sum((alpha+1)*log.1mZ)
        denom1=sum(lgamma(alpha.unique))*n+n*lgamma(beta)
        denom2=sum((sum(alpha.unique)+beta)*log.denom)
        return (nom1+nom2-denom1-denom2)
      }else{
        return (sum(dbeta(Z,alpha.unique,beta,log=TRUE)))
      }
    }
  }
  A=matrix(rep(0,p+L+1),nrow=1)
  A[p+L+1]=1
  B=as.matrix(0)
  res <- maxLik(llk, start = starter,constraints=list(ineqA=A, ineqB=B),iterlim=10000)
  return (res)
}



#output:joint p value of profile likelihood ratio/Wald method for that block (consiting of L CpG sites)
ht<-function(Z,X.no.int){ #no intercept
  L=ncol(as.matrix(Z))
  X=cbind(1,X.no.int)
  res.constraint=get.mle.reg(Z,X[,-2]) #the first variable is the interesting phenotype
  res.unconstraint=get.mle.reg(Z,X)
  p.lr=(1-pchisq(2*(res.unconstraint$maximum-res.constraint$maximum),1)) #
  p=ncol(as.matrix(X))
  index=L+1
  Fisher=-(res.unconstraint$hessian[index,index])
  p.Wald=1-pchisq(res.unconstraint$estimate[index]^2*Fisher,1)
  return (list(p.Wald=p.Wald, p.lr=p.lr))
  
}


#############################################
###functions 
#############################################
#input: x, y are both n by p matrices 
#output: a vector of length p indicating the column-wise correlation between x and y.

colCors = function(x, y) {
  sqr = function(x) x*x
  if(!is.matrix(x)||!is.matrix(y)||any(dim(x)!=dim(y)))
    stop("Please supply two matrices of equal size.")
  x   = sweep(x, 2, colMeans(x))
  y   = sweep(y, 2, colMeans(y))
  cor = colSums(x*y) /  sqrt(colSums(sqr(x))*colSums(sqr(y)))
  return(cor)
}

#given a vector y, and a randomly generated vector as the same length of y
#return a vector such that cor(y,output)=pho
complement <- function(y, rho, x) {
  if (missing(x)) x <- rnorm(length(y)) # Optional: supply a default if `x` is not given
  y.perp <- residuals(lm(x ~ y))
  return (rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2))
}


#cor(x,rowMeans(y))
#logit function
logit <- function(x) return (log(x/(1-x)))

