#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix g_cpp(NumericVector x){
  NumericVector thresh(x.length(),exp(-20));
    //Function ifelse_cpp("ifelse");
  NumericMatrix M(x.length(),1);
  M(_,0)= ifelse(x>thresh,x,thresh);
  return M;
}
// [[Rcpp::export]]
NumericVector basis_cpp(NumericVector tau,double l, double u){
  NumericVector Q(tau.length());
  
  if(l<0.5){
    Q = Rcpp::qnorm5(tau)-R::qnorm5(u,0,1,TRUE,FALSE);
    Q[tau>u] = 0;
    Q[tau<l] = R::qnorm5(l,0,1,TRUE,FALSE)-R::qnorm5(u,0,1,TRUE,FALSE);
  } 
  if(l>=0.5){
    Q = Rcpp::qnorm5(tau)-R::qnorm5(l,0,1,TRUE,FALSE);
    Q[tau<l] = 0;
    Q[tau>u] = R::qnorm5(u,0,1,TRUE,FALSE)-R::qnorm5(l,0,1,TRUE,FALSE);
  }
  
  return Q;
}

// [[Rcpp::export]]
NumericMatrix makeB_cpp(NumericVector tau, int L){
  NumericMatrix B(tau.length(),L);
  for(int j=0;j<L;j++){
    B(_,j)=basis_cpp(tau,(double) j/L,(double) (j+1)/L);
  }
  return B;
}


// [[Rcpp::export]]
arma::mat g_cpp_arma(arma::vec x){
  arma::vec thresh(x.n_elem);
  thresh.fill(exp(-20));
  //Function ifelse_cpp("ifelse");
  arma::mat M(x.n_elem,1);
  x.elem(find(x<=thresh)).fill(exp(-20));
  M.col(0)= x;
  return M;
}


// [[Rcpp::export]]
arma::vec basis_cpp_arma(NumericVector tau,double l, double u){
  NumericVector Q(tau.length());
  if(l<0.5){
    Q = Rcpp::qnorm5(tau)-R::qnorm5(u,0,1,TRUE,FALSE);
    Q[tau>u] = 0;
    Q[tau<l] = R::qnorm5(l,0,1,TRUE,FALSE)-R::qnorm5(u,0,1,TRUE,FALSE);
  } 
  if(l>=0.5){
    Q = Rcpp::qnorm5(tau)-R::qnorm5(l,0,1,TRUE,FALSE);
    Q[tau<l] = 0;
    Q[tau>u] = R::qnorm5(u,0,1,TRUE,FALSE)-R::qnorm5(l,0,1,TRUE,FALSE);
  }
  arma::vec Qa=Q;
  return Qa;
}



// [[Rcpp::export]]
arma::mat makeB_cpp_arma(NumericVector tau, int L){
  arma::mat B(tau.length(),L);
  
  for(int j=0;j<L;j++){
    B.col(j)=basis_cpp_arma(tau,(double) j/L,(double) (j+1)/L);
  }
  return B;
}


// [[Rcpp::export]]
List makeab_cpp(arma::mat X,arma::mat theta,int L){
  arma::mat T(L,X.n_rows);
  T.zeros();
  arma::mat a = T;
  for(int j=0;j<L;j++){
    T.row(j) = trans(g_cpp_arma(X*theta.col(j)));
  }
  arma::mat b = T;
  NumericVector k(L+1);
  k = arma::linspace(0,1,L+1);
  arma::mat qk=makeB_cpp_arma(k,L)*T;
  arma::vec phik = Rcpp::qnorm5(k);
  for(int j=0;j<L;j++){
    if(k[j]<0.5){
      a.row(j) = qk.row(j+1)-phik(j+1)*T.row(j);
    }
    if(k[j]>=0.5){
      a.row(j) = qk.row(j)-phik(j)*T.row(j);
    }
  }
  return List::create(Named("a") = a,Named("b") = b,Named("qk") = qk);
}

// [[Rcpp::export]]
NumericVector dSplitNorm_cpp(arma::vec y,arma::mat X,arma::vec med,
                         arma::mat theta, bool logts=false){
 
  NumericVector dense(y.n_elem);
  int L = theta.n_cols;
  NumericVector ny(y.n_elem);
  ny = y-X*med;
  y=ny;
  NumericVector nsid(y.n_elem);
  nsid = g_cpp_arma(X*theta.col(0));
  if(L==1){
    dense = Rcpp::dnorm(ny,0,nsid);//dnorm
  }
  if(L>1){
    List ab = makeab_cpp(X,theta,L);
    arma::mat abqk=ab["qk"];
    arma::mat a= ab["a"];
    arma::mat b = ab["b"];
    arma::uvec these;
    unsigned int id;
    for(int j=0;j<L;j++){
      these = arma::find((abqk.row(j).t()<y) && (y<abqk.row(j+1).t()));
      for(unsigned int i=0; i<these.n_elem; i++){
        id =these(i);
        dense(id) = dense(id)+R::dnorm(ny(id),a(j,id),b(j,id),false);
      }
    }
  }
  if(logts){dense = log(dense);}
  
  return dense;
}

// [[Rcpp::export]]
NumericVector qSplitNorm_cpp(NumericVector tau,arma::mat X,arma::vec med,
                             arma::mat theta){
// tau is an n-vector of probabilities
// X is an 1xp matrix of covariates
// median is an p-vector of median coefficnets
// theta is pxL array with X%*%theta[,j]>0 for all j
  int L = theta.n_cols;
  NumericVector mmm(X.n_rows);
  mmm = X*med;
  NumericVector qqq(tau.length());
  if(L==1){
    NumericVector sss(tau.length());
    sss = g_cpp_arma(X*theta.col(0));
     qqq = Rcpp::qnorm5(tau,mmm,sss);
  }
  if(L>1){
    arma::mat T(L,X.n_rows);
    T.zeros();
    for(int j=0;j<L;j++){
      T.row(j) = trans(g_cpp_arma(X*theta.col(j)));
    }
    //arma::mat BBB = makeB_cpp_arma(tau,L);
    NumericVector BBB(tau.length());
    BBB= makeB_cpp_arma(tau,L)*T;
    qqq =rep_len(mmm,tau.length()) + BBB;
  }
  
  return qqq;
}

// [[Rcpp::export]]
NumericVector pSplitNorm_cpp(arma::vec y,arma::mat X,arma::vec med,
                             arma::mat theta){
// y is an n-vector of responses
// X is an nxp matrix of covariates
// median is an p-vector of median coefficnets
// theta is pxL array with X%*%theta[,j]>0 for all j
  NumericVector ppp(y.n_elem);
  ppp = y*0.0;
  int L = theta.n_cols;
  NumericVector ny(y.n_elem);
  ny = y-X*med;
  y=ny;
  NumericVector nsid(y.n_elem);
  nsid = g_cpp_arma(X*theta.col(0));
  if(L==1){
    ppp = Rcpp::pnorm5(ny,0.0,nsid);//pnorm
  }
  
  if(L>1){
    List ab = makeab_cpp(X,theta,L);
    arma::mat qk=ab["qk"];
    arma::mat a= ab["a"];
    arma::mat b = ab["b"];
    arma::uvec these;
    
    unsigned int id;
    for(int j=0;j<L;j++){
      these = arma::find(y>qk.row(j+1).t());
      
      for(unsigned int i=0; i<these.n_elem; i++){
        id =these(i);
        ppp(id) = ppp(id)+R::pnorm(qk(j+1,id),a(j,id),b(j,id),true,false)
          -R::pnorm(qk(j,id),a(j,id),b(j,id),true,false);
      }
      
      these = arma::find((y>=qk.row(j).t())&&(y<qk.row(j+1).t()));
      for(unsigned int i=0; i<these.n_elem; i++){
        id =these(i);
        ppp(id) = ppp(id)+R::pnorm(y(id),a(j,id),b(j,id),true,false)
          -R::pnorm(qk(j,id),a(j,id),b(j,id),true,false);
      }
      
    }
  }
  return ppp;
}

// [[Rcpp::export]]
arma::mat covexp_cpp(arma::mat d,double range,double ratio){
  arma::mat C(size(d));
  C.zeros();
  C = exp(-d/range);
  if(ratio<1.0){
    C = ratio*C + (1-ratio)*eye(size(d));
  }
  return C;
}

// [[Rcpp::export]]
List UPDATEMEDIAN(arma::mat Y,arma::cube X,arma::mat med,arma::vec mn_med,arma::vec tau_med,
                  arma::mat canz,arma::vec curll,arma::vec MHmed,arma::cube theta,
                  arma::vec attmed,arma::vec accmed,arma::mat P,arma::mat Pz,arma::mat z){

  for(unsigned int k=0; k<Y.n_cols;k++){
    double Vz = Pz(k,k);
    arma::rowvec pzk = Pz.row(k);
    pzk.shed_col(k);
    arma::mat zk = z;
    zk.shed_col(k);
    //NumericVector zck1;
    arma::vec zck;
    arma::mat Mz;
    Mz= - pzk*zk.t()/Vz;
    
    arma::vec Mzv = arma::trans(Mz.row(0));
    arma::vec sig(Mz.n_elem);
    sig.fill(1/sqrt(Vz));
    
    double VVV = 0.0;
    arma::rowvec RRR(med.n_cols-1);
    double MMM = 0.0;
    double canll = 0.0;
    arma::mat canmed = med;
    arma::rowvec medk;
    arma::rowvec pk = P.row(k);
    pk.shed_col(k);
    
    arma::vec yk = Y.col(k);
    arma::mat xk = X.col(k);
    arma::mat thetak = theta.col(k);
    arma::vec canmedk(med.n_rows);
    arma::vec canzk(canz.n_rows);
    
    arma::vec tem1(canz.n_rows);
    arma::vec tem2(canz.n_rows);
    
    //NumericVector canzk1;
    double R = 0.0;
    for(unsigned int j=0;j< X.n_slices; j++){
      attmed(j) =attmed(j)+1; 
      VVV = tau_med(j)*P(k,k);
      medk = med.row(j);
      medk.shed_col(k);
      RRR = medk - mn_med(j);
      MMM = tau_med(j)*P(k,k)*mn_med(j) - tau_med(j)*sum(pk*RRR.t());
      
      canmed = med;
      canmed(j,k) = R::rnorm(med(j,k),MHmed(j));
      canmedk = canmed.col(k);
      
      canll = sum(dSplitNorm_cpp(yk,xk,canmedk,thetak,true)); //wrong
      
      canzk =Rcpp::qnorm5(pSplitNorm_cpp(yk,xk,canmedk,thetak));
      canz.col(k) = canzk;
      zck = z.col(k);
      
      tem1 = log(arma::normpdf(canzk,Mzv,sig)) - log(arma::normpdf(zck,Mzv,sig));
      tem2 = -log(arma::normpdf(canzk)) + log(arma::normpdf(zck));
      R = sum(tem1)+sum(tem2)+canll-curll(k)+
        log(arma::normpdf(canmed(j,k),MMM/VVV,1/sqrt(VVV)))
        -log(arma::normpdf(med(j,k),MMM/VVV,1/sqrt(VVV)));
      R = exp(R);
      if(R!=0.0){
        if(R::runif(0.0,1.0)<R){
          med = canmed;
          curll(k)=canll;
          z.col(k) = canzk;
          accmed(j) = accmed(j)+1;
        }
      }
    }
    
  }
  
  return List::create(Named("med") = med,Named("curll") = curll,Named("z") = z,Named("canz") = canz,
                            Named("attmed") = attmed,Named("accmed") = accmed);
  
}


// [[Rcpp::export]]
List UPDATETHETA(arma::mat Y,arma::cube X,arma::mat med,arma::mat mn_theta,arma::mat tau_theta,
                  arma::mat canz,arma::vec curll,arma::mat MHtheta,arma::cube theta,
                  arma::mat atttheta,arma::mat acctheta,arma::mat P,arma::mat Pz,arma::mat z){
  
  for(unsigned int k=0; k<Y.n_cols;k++){
    double Vz = Pz(k,k);
    arma::rowvec pzk = Pz.row(k);
    pzk.shed_col(k);
    arma::mat zk = z;
    zk.shed_col(k);
    //NumericVector zck1;
    arma::vec zck;
    arma::mat Mz;
    Mz= - pzk*zk.t()/Vz;
    
    arma::vec Mzv = arma::trans(Mz.row(0));
    arma::vec sig(Mz.n_elem);
    sig.fill(1/sqrt(Vz));
    
    double VVV = 0.0;
    arma::rowvec RRR(med.n_cols-1);
    double MMM = 0.0;
    double canll = 0.0;
    arma::cube cantheta = theta;
    arma::mat canthetak;
    arma::rowvec pk = P.row(k);
    pk.shed_col(k);
    
    arma::vec yk = Y.col(k);
    arma::mat xk = X.col(k);
    arma::mat thetak = theta.col(k);
    //arma::vec canmedk(med.n_rows);
    arma::vec canzk(canz.n_rows);
    
    arma::vec tem1(canz.n_rows);
    arma::vec tem2(canz.n_rows);
    arma::vec medk;
    //NumericVector canzk1;
    double R = 0.0;
    for(unsigned int l=0; l<theta.n_slices;l++){
      
    for(unsigned int j=0;j< X.n_slices; j++){
      atttheta(j,l) =atttheta(j,l)+1; 
      VVV = tau_theta(j,l)*P(k,k);
      
      thetak = theta.subcube(arma::span(j),arma::span(),arma::span(l));
      thetak.shed_col(k);
      RRR = thetak - mn_theta(j,l);
      MMM = tau_theta(j,l)*P(k,k)*mn_theta(j,l) - tau_theta(j,l)*sum(pk*RRR.t());
      
      cantheta = theta;
      cantheta(j,k,l) = R::rnorm(theta(j,k,l),MHtheta(j,l));
      canthetak = cantheta.col(k);
      medk = med.col(k);
      canll = sum(dSplitNorm_cpp(yk,xk,medk,canthetak,true)); //wrong
      canzk =Rcpp::qnorm5(pSplitNorm_cpp(yk,xk,medk,canthetak));
      canz.col(k) = canzk;

      zck = z.col(k);
      
      tem1 = log(arma::normpdf(canzk,Mzv,sig)) - log(arma::normpdf(zck,Mzv,sig));
      tem2 = -log(arma::normpdf(canzk)) + log(arma::normpdf(zck));
      R = sum(tem1)+sum(tem2)+canll-curll(k)+
        log(arma::normpdf(cantheta(j,k,l),MMM/VVV,1/sqrt(VVV)))
        -log(arma::normpdf(theta(j,k,l),MMM/VVV,1/sqrt(VVV)));
      R = exp(R);
      if(R!=0.0){
        if(R::runif(0.0,1.0)<R){
          theta = cantheta;
          curll(k)=canll;
          z.col(k) = canzk;
          acctheta(j,l) = acctheta(j,l)+1;
        }
      }
    }
    }
    
  }
  
  return List::create(Named("theta") = theta,Named("curll") = curll,Named("z") = z,Named("canz") = canz,
                            Named("atttheta") = atttheta,Named("acctheta") = acctheta);
  
}



//Rcpp::sourceCpp("SpatialQR_cpp.cpp")

///*** R
//timesTwo(42)
//*/
