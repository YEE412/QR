#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec CD(arma::mat x, arma::vec y, arma::vec beta,double toler, 
       int maxit, double tau,double lambda) {
  //beta is the initial value of iteration
  //toler is the abs value difference threshold
  //maxit is the toplimit of iteration times
  //tau is the quantile level
  //lambda is the regularization parameter
  
  int n = x.n_rows;
  int p = x.n_cols;
  int iter = 1;
  double error = 10000;
  //res = y - x汕 (after sorted)
  arma::vec res;
  //b0 is the intercept
  double b0;
  //r = y - b0 - x汕
  arma::vec r;
  //牟i =  而 ,when ri > 0; 牟i =  1 - 而 ,when ri < 0
  arma::vec theta;
  //using weighted median of z to update 汕m
  arma::vec z;
  //sortedindex is the sorted order of z
  arma::uvec sortedindex;
  //zsorted is z after sorted
  arma::vec zsorted;
  //w is the weight of z
  arma::vec w;
  //wsorted = w(sortedindex)
  arma::vec wsorted;
  //compute abs value difference by sum(abs(betaold-betanew))
  arma::vec betaold, betanew = beta;
  //betaQR is the estimator (without intercept) obtained from non-regularized quantile regression
  arma::vec betaQR;
  if(lambda != 0){
    betaQR = CD(x,y,beta,toler,maxit,tau,0);
  }
  //istar is the index of weight median of z
  int istar;
  
  //set two stopping criteria
  while(iter <= maxit && error > toler){
    betaold = betanew;
    
    //replace b0 by the level-而 sample quantile of the residuals yi ??? Xi汕
    res = sort(y-x*betanew);
      //get (the level-而 sample quantile of residuals) by using linear interpolation
    b0 = ((n-1)*tau-floor((n-1)*tau))*res(ceil((n-1)*tau))+(1-((n-1)*tau-floor((n-1)*tau)))*res(floor((n-1)*tau));
    
    //r = y - b0 - x汕
    r=y-b0-x*betanew;
    //牟i =  而 ,when ri > 0; 牟i =  1 - 而 ,when ri < 0
    theta = (1-arma::sign(r))/2*(1-tau)+(1+arma::sign(r))*tau/2;
    
    //for fixed m = 1,2,＃＃p, using weighted median of z to update 汕m
    for (int m = 0; m < p; ++m){
      z = (r + x.col(m) * betanew(m)) / x.col(m);
      if (lambda != 0){
        z.insert_rows(n,1);//add 0 to the end of z
      }
      sortedindex = arma::sort_index(z);
      zsorted = z(sortedindex);
      
      //wi = abs(x_im * theta_m)
      w = abs(x.col(m)%theta);
      if (lambda != 0){
        w.insert_rows(n,1);//add lambda/betaQR(m)/betaQR(m) to the end of w
        w(n)=abs(lambda/betaQR(m)/betaQR(m));
      }
      wsorted = w(sortedindex);
      arma::uvec index = arma::find(arma::cumsum(wsorted)>(arma::sum(wsorted)/2),1,"first");
      istar = int(index(0));
      
      betanew(m) = zsorted(istar);
    }
    error=arma::sum(abs(betanew-betaold));
    iter++;
  }
  return betanew;
}
