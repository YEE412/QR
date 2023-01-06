#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec QRPMMCPP (arma::mat xr,arma::vec y,arma::vec beta,
                    double toler,int maxit,double tau,double lamda){
  //beta is the initial value of iteration
  //toler is the abs value difference threshold
  //maxit is the toplimit of iteration times
  //tau is the quantile level
  //lambda is the regularization parameter
  
  arma:: mat x = (xr),product;
  //xt is transpose of x
  arma:: mat xt;
  arma:: vec W,v,r,E;
  //refresh beta by delta.betanew=betaold-delta
  arma:: vec delta;
  arma:: vec betaold ,betanew = beta;
  int n=x.n_rows;
  x.insert_cols( 0, arma::ones(n) );
  int p=x.n_cols;
  // while error > toler,iteration continue
  double error=10000,epsilon=0.9999;
  int iteration=1;
  double epsilonprime=toler*p/2;
  //betaQR is the estimator (without intercept) obtained from non-regularized quantile regression
  arma::vec qrbeta;
  product.ones(p,n);
  if(lamda != 0){
    qrbeta = QRPMMCPP(xr,y,beta,toler,maxit,tau,0);
    qrbeta.shed_row(0);
    qrbeta.insert_rows(0,1);
  }
  xt=x.t();
  //set two stopping criteria
  while (iteration<=maxit&& error>toler)
  {
    betaold = betanew;
    r=y-x*betaold;
    v=1-2*tau-r/(arma::abs(r)+epsilon);
    if(lamda != 0 ){
      E=lamda*arma::abs(betaold)/betaold/(epsilonprime+arma::abs(betaold))/qrbeta/qrbeta;
      E(0)=0;
    }
    
    W=1/(epsilon+arma::abs(r));
    int i;
    //parallel computing
#pragma omp parallel for schedule(static) \
    private(i) shared(product,xt,W)
      for (i=0;i<n;i++)
      {product.col(i)=xt.col(i)*W(i);}
      //delta is different when lambda =0 and lambda!=0
      if(lamda != 0){
        delta=arma::solve(product*x-diagmat(E),xt*v-E%betaold);
      }
      else{
        delta=arma::solve(product*x,xt*v);
      }
      //iterate beta
      betanew=betaold-delta;
      //compute error
      error=arma::sum(arma::abs(delta));
      iteration++;
  }
  return betanew;
}

