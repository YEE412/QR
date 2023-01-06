QR.lasso.cd.r = function(X,y,tau,lambda,beta,maxit,toler)
{
  
  if(missing(toler)){
    toler = 1e-3
    
  }
  
  if(missing(maxit)){
    maxit = 200
  }
  
  if(missing(beta)){
    beta=solve(t(X)%*%X,t(X)%*%y)
  }
  
  if(missing(lambda)){
    lambda=1
  }
  
  beta0=QRCD_R(X,y,beta,toler,maxit,tau)
  betah=QRPCD_R(X,y,beta,beta0,toler,maxit,tau,lambda)
  beta=betah
  b=quantile(y-X%*%beta,tau)
  return(list(beta=beta,b=b))
}