QRMM_R <- function(x,y,beta,toler,maxit,tau){
  n = nrow(x)
  x = cbind(rep(1,n),x)
  p = ncol(x)
  error=10000
  epsilon=0.9999
  iteration=1
  product = matrix(1,p,n)
  xt = t(x)
  
  while(iteration<=maxit & error>toler){
    betaold=beta
    
    r=y-x%*%beta
    v=1-2*tau-r/(abs(r)+epsilon)
    
    W=1/(epsilon+abs(r))
    
    for (i in 1:n)
    {	product[,i]=xt[,i]*W[i]}
    
    delta=solve(product%*%x,xt%*%v);
    beta=beta-delta
    
    
    error=sum(abs(delta))
    iteration = iteration+1
  }
  return (beta)
}