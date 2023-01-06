QRPMM_R <- function(x,y,beta,qrbeta,toler,maxit,tau,lamda){
  n = nrow(x)
  x = cbind(rep(1,n),x)
  p = ncol(x)
  error=10000
  epsilon=0.9999
  iteration=1
  product = matrix(1,p,n)
  epsilonprime = toler*p/2
  xt = t(x)
  qrbeta = append(qrbeta,0,0)
  
  
  while(iteration<=maxit & error>toler){
    betaold=beta
    
    r=y-x%*%beta
    v=1-2*tau-r/(abs(r)+epsilon)
    E=lamda*abs(beta)/beta/(epsilonprime+abs(beta))/qrbeta/qrbeta
    E[1] = 0
    
    W=1/(epsilon+abs(r))
    for (i in 1:n)
    {product[,i]=xt[,i]*W[i]}

    delta=solve(product%*%x-diag(c(E)),xt%*%v-E*beta);
    beta=beta-delta
    
    
    error=sum(abs(delta))
    iteration = iteration+1
  }
  betanew = beta
  return (betanew)
}