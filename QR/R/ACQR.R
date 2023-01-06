dmatrix_ACQ <- function(num,dim){
  return(matrix(rep(num,dim),nrow=1,ncol=dim))
}

###########################CompQreg and AdpComptQreg#############
CompQreg <- function(y,X,n,p,tao1) {
  A1 <-cbind(dmatrix_ACQ(0,p),dmatrix_ACQ(0,p),tao1*dmatrix_ACQ(1,n),(1-tao1)*dmatrix_ACQ(1,n))  #1*(2p+2n)
  B1 <- cbind(dmatrix_ACQ(1,p),dmatrix_ACQ(-1,p),dmatrix_ACQ(0,n),dmatrix_ACQ(0,n))
  B2 <- cbind(X,-X,diag(1,nrow=n,ncol=n),diag(-1,nrow=n,ncol=n))
  B <- rbind(B1,B2)                                                        #(1+n)*(2p+2n)
  H <- rbind(0,matrix(y,ncol=1))                                           #(1+n)*1
  prod <- lp("min", A1, B,"==",H,compute.sens=TRUE)
  beta11 <- prod$solution[1:p]
  beta22 <- prod$solution[(p+1):(2*p)]
  betaresult <- beta11-beta22
  return(betaresult)
}

AdpComptQreg<-function(y,X,n,p,lambda,tao1){
  w <- rep(0,p)
  tildebeta <-CompQreg(y,X,n,p,tao1)

  for(i in 1:p){
    w[i] <- 1/abs(tildebeta[i])
  }
  w[w==Inf]<-999999999 ## Replacing values estimated as Infinite for 999999999
  w <- matrix(w,nrow=1,ncol=p)

  A <- cbind(matrix(lambda*w,nrow=1),matrix(lambda*w,nrow=1),tao1*dmatrix_ACQ(1,n),(1-tao1)*dmatrix_ACQ(1,n))
  B1 <- cbind(dmatrix_ACQ(1,p),dmatrix_ACQ(-1,p),dmatrix_ACQ(0,n),dmatrix_ACQ(0,n))
  B2 <- cbind(X,-X,diag(1,nrow=n,ncol=n),diag(-1,nrow=n,ncol=n))
  B <- rbind(B1,B2)
  H <- rbind(0,matrix(y,ncol=1))
  prod <- lp("min", A, B,"==",H)
  beta11 <- prod$solution[1:p]
  beta22 <- prod$solution[(p+1):(2*p)]
  betaresult <- beta11-beta22
  return(betaresult)
}

###########################find lambda for AdpComptQreg#############
rho_tao_ACQ  <- function(tao1,u) {
  if(u<0) { return(u*(tao1-1))}
  else {return(u*tao1) }
}


HBIC_ACQ <- function(y,X,beta,lambda,tao1){
  n = length(y)
  p = dim(X)[2]
  M_cardi = length(which(beta[1:p]!=0))
  C_n = log(log(n))
  bic = 0
  for(i in 1:n){
    bic <- bic + rho_tao_ACQ(tao1,y[i]-t(X[i,])%*%beta)
  }
  HBIC = log(bic) + M_cardi*C_n*log(p)/n
  return(HBIC)
}


findlambda_ACQ <- function(y,X,beta,lambda,tao1) {
  minBIC=1e10
  for (i in 1:length(lambda)){
    bic=HBIC_ACQ(y,X,beta[,i],lambda[i],tao1)
    if(bic<minBIC) {
      minBIC = bic
      lam=lambda[i]}
  }
  return(lam)
}

###########################AdpComptQreg with HBIC#############
ACQR <- function(y,X,lamseq,ncore,tao1){
  n=length(y)
  p=dim(X)[2]
  betaseq <- matrix(0,nrow=p,ncol=length(lamseq))
  if(ncore>1){
  cl <- makeCluster(ncore)
  registerDoParallel(cl)
  j=0
  betaseq<-foreach(j = 1:length(lamseq),.packages = 'lpSolve',.combine = cbind,.export=c("AdpComptQreg","CompQreg","dmatrix_ACQ","rho_tao_ACQ")) %dopar% {
    AdpComptQreg(y,X,n,p,lamseq[j],tao1) }
  }
  else{
    for(j in 1:length(lamseq)){
      betaseq[,j]<-AdpComptQreg(y,X,n,p,lamseq[j],tao1) }
  }
  lambda=findlambda_ACQ(y,X,betaseq,lamseq,tao1)
  beta <-AdpComptQreg(y,X,n,p,lambda,tao1)
  return(beta)
  stopImplicitCluster()
  stopCluster(cl)
}
