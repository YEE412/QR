dmatrix_ACCQR <- function(num,dim){
  return(matrix(rep(num,dim),nrow=1,ncol=dim))
}

###########################CompQreg and CompadaptiveQreg#############
genmatrix <- function(y,X,n,p,K){
  B<- cbind(dmatrix_ACCQR(1,p),dmatrix_ACCQR(0,K),dmatrix_ACCQR(-1,p),dmatrix_ACCQR(0,K),matrix(0,nrow=1,ncol=(2*K*n)))
  Z <-cbind(X,matrix(0,nrow=n,ncol=K))
  for(k in 1:K){
    Zk=Z
    Zk[,(p+k)]=1
    Bk<-cbind(Zk,-Zk,matrix(0,nrow=n,ncol=(2*(k-1)*n)),diag(1,nrow=n,ncol=n),diag(-1,nrow=n,ncol=n),matrix(0,nrow=n,ncol=(2*(K-k)*n)))
    B <-rbind(B,Bk)
  }    #dimB (1+K*n)*(2(p+k)+2*k*n)
  H<-matrix(c(0,rep(y,K)),ncol=1)   #(1+n*k)*1
  return(list(B=B,H=H))
}

CCompQreg <- function(y,X,n,p,K) {
  tao <- seq(1/(K+1),K/(K+1),by=1/(K+1))
  A1 <- cbind(dmatrix_ACCQR(0,p+K),dmatrix_ACCQR(0,p+K))
  for(k in 1:K){
    A1 <- cbind(A1,dmatrix_ACCQR(1,n)*tao[k],dmatrix_ACCQR(1,n)*(1-tao[k]))
  }
  B <- genmatrix(y,X,n,p,K)$B
  H <- genmatrix(y,X,n,p,K)$H
  prod <- lp("min", A1, B,"==",H)
  beta1 <- prod$solution[1:p]
  beta2 <- prod$solution[(p+K+1):(2*p+K)]
  beta <- beta1-beta2
  return(beta)}

AdpCComptQreg<-function(y,X,n,p,lambda,K) {
  tao <- seq(1/(K+1),K/(K+1),by=1/(K+1))
  w <- rep(0,p)
  tildebeta <-CCompQreg(y,X,n,p,K)

  for(i in 1:p){
    w[i] <- 1/(tildebeta[i]^2)
  }

  w[w==Inf]<-999999999 ## Replacing values estimated as Infinite for 999999999
  w <- matrix(w,nrow=1,ncol=p)

  A <- cbind(matrix(lambda*w,nrow=1),dmatrix_ACCQR(0,K),matrix(lambda*w,nrow=1),dmatrix_ACCQR(0,K))
  for(k in 1:K){
    A <- cbind(A,dmatrix_ACCQR(1,n)*tao[k],dmatrix_ACCQR(1,n)*(1-tao[k]))}
  B <- genmatrix(y,X,n,p,K)$B
  H <- genmatrix(y,X,n,p,K)$H
  prod <- lp("min", A, B,"==",H)
  beta1 <- prod$solution[1:(p+K)]             ##include intercept b to calculate bic
  beta2 <- prod$solution[(p+K+1):(2*p+2*K)]
  beta <- beta1-beta2
  return(beta)
}

rho_tao_ACCQR <- function(tao1,u) {
  if(u<0) { return(u*(tao1-1))}
  else {return(u*tao1) }
}

HBIC_ACCQ <- function(y,X,beta,lambda,K){
  tao <- seq(1/(K+1),K/(K+1),by=1/(K+1))
  n = length(y)
  p = dim(X)[2]
  M_cardi = length(which(beta[1:p]!=0))
  C_n = log(log(n))
  bic = 0
  for(i in 1:n){
    for(k in 1:K) {
      bic <- bic + rho_tao_ACCQR(tao[k],y[i]-t(X[i,])%*%beta[1:p]-beta[p+k])
    }
  }
  HBIC = log(bic/K) + M_cardi*C_n*log(p)/n


  return(HBIC)
}


findlambda_ACCQ <- function(y,X,beta,lambda,K) {
  minBIC=1e10
  for (i in 1:length(lambda)){
    bic=HBIC_ACCQ(y,X,beta[,i],lambda[i],K)
    if(bic<minBIC) {
      minBIC = bic
      lam=lambda[i]}
  }
  return(lam)
}

ACCQR <- function(y,X,lamseq,ncore,K){
  n=length(y)
  p=dim(X)[2]
  betaseq <- matrix(0,nrow=(p+K),ncol=length(lamseq))

  if(ncore>1){
  cl <- makeCluster(ncore)
  registerDoParallel(cl)
  j=0
  betaseq<-foreach(j = 1:length(lamseq),.packages = 'lpSolve',.combine = cbind,.export=c("AdpCComptQreg","CCompQreg","dmatrix_ACCQR","genmatrix","rho_tao_ACCQR")) %dopar% {
    AdpCComptQreg(y,X,n,p,lamseq[j],K) }
  }
  else{
    for(j in 1:length(lamseq)){
      betaseq[,j]<-AdpCComptQreg(y,X,n,p,lamseq[j],K) }
  }
  lambda = findlambda_ACCQ(y,X,betaseq,lamseq,K)
  beta <-AdpCComptQreg(y,X,n,p,lambda,K)[1:p]
  return(beta)
  stopImplicitCluster()
  stopCluster(cl)
}
