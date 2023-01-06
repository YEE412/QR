QR.lasso.cd.cpp = function(X,y,tau,lambda,beta,maxit,toler)
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
  
  beta=CD(X,y,beta,toler,maxit,tau,lambda)
  b=quantile(y-X%*%beta,tau)
  return(list(beta=beta,b=b))
}