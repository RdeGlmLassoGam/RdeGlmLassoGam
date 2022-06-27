
// includes from the plugin
#include <RcppArmadillo.h>
#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes


// declarations
extern "C" {
SEXP fit_gam_sp_cpp( SEXP Ry, SEXP RB, SEXP RrS, SEXP Rfamily) ;
}

// definition

SEXP fit_gam_sp_cpp( SEXP Ry, SEXP RB, SEXP RrS, SEXP Rfamily ){
BEGIN_RCPP
// Ry, RB, RrS, Rfamily

// fit.am.sp

int max_count = 50;
arma::colvec y = Rcpp::as<arma::colvec>(Ry);
arma::mat B = Rcpp::as<arma::mat>(RB);
arma::mat rS = Rcpp::as<arma::mat>(RrS);
Rcpp::List family(Rfamily);

std::string link = Rcpp::as<std::string>(family["link"]);
Function linkfun(Rcpp::as<Function>(family["linkfun"]));
Function linkinv(Rcpp::as<Function>(family["linkinv"]));
Function mu_eta(Rcpp::as<Function>(family["mu.eta"]));
Function variance(Rcpp::as<Function>(family["variance"]));

int q, n;
q = B.n_cols;
n = B.n_rows;

arma::mat X1 = arma::join_cols(B,rS);
arma::colvec eta(n);

// set initial values
if (link.compare("log")==0){
    eta = Rcpp::as<arma::colvec>(linkfun(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y+0.001))));
} else if (link.compare("logit")==0){
    eta = Rcpp::as<arma::colvec>(linkfun((Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap((y+0.5)/2)))));
} else if (link.compare("identity")==0){
    eta = Rcpp::as<arma::colvec>(linkfun(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y))));
}
//std::cout<<link.compare("log");

double norm = 0;
double old_norm = 1;
int count = 1;
arma::colvec mu(n);
arma::colvec weight(n+q);
arma::colvec mu_eta_val(n);
arma::colvec variance_val(n);
arma::colvec z(n+q);
arma::colvec m(n+q);
arma::colvec coef(q);

// while loop
while (((std::abs(norm-old_norm)/norm)>1e-4)&&(count<=max_count)){
    mu = Rcpp::as<arma::colvec>(linkinv(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(eta))));
    mu_eta_val = Rcpp::as<arma::colvec>(mu_eta(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(eta))));
    variance_val = Rcpp::as<arma::colvec>(variance(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(mu))));
    weight = join_cols(mu_eta_val/arma::sqrt(variance_val),arma::ones(q));
    z = join_cols((y-mu)/mu_eta_val+eta,arma::zeros(q));
    coef = arma::solve(arma::diagmat(weight)*X1,weight%z);
    m = X1*coef;
    eta = m.subvec(0,n-1);
    old_norm = norm;
    norm = arma::norm((z.subvec(0,n-1)-eta),2)/std::sqrt(n);
    count = count+1;
}
return Rcpp::List::create(Rcpp::Named("coefficients")=coef,Rcpp::Named("fitted.values")=Rcpp::as<arma::colvec>(linkinv(eta)));
END_RCPP
}



