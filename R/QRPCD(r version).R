QRPCD_R <-function(x,y,beta,beta0,toler,maxit,tau,lamda) 
{
  n = nrow(x)
  p = ncol(x)
  error=10000
  iteration=1
  
  while (iteration<=maxit & error>toler)
  {
    betaold=beta
    uv=sort(y-x%*%beta)
    
    #replace b0 by the level-tau sample quantile of residuals
    #get (the level-tau sample quantile of residuals) by using linear interpolation
    quantile=(n-1)*tau-floor((n-1)*tau)
    u=quantile*uv[ceiling((n-1)*tau) + 1]+(1-quantile)*uv[floor((n-1)*tau) + 1]
    
    r=y-u-x%*%beta
    signw=(1-sign(r))/2*(1-tau)+(sign(r)+1)*tau/2;
    
    for (j in 1:p)
    {          	
      z=(r+beta[j]*x[,j])/x[,j]
      z = append(z,0)
      order_z = order(z)
      sortz=z[order_z]
      newX=x[,j]*signw
      newX = append(newX,0)
      newX[n+1]=(lamda/beta0[j]/beta0[j])
      w = abs(newX[order_z])
      beta[j] = sortz[(cumsum(w) > sum(w)/2)][1]
    }
    error=sum(abs(beta-betaold))
    iteration = iteration + 1
  }
  
  return (beta)
}