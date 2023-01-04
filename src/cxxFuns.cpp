// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' @title Log-Likelihood of Matrix-Variate Normal Data
//' @description \loadmathjax
//' Log-likelihood of i.i.d data from \mjsdeqn{N(0, \Psi, \Sigma)}.
//'
//' @param x Array.
//' @param Psi,Sig Row and column convariance matrice.
//'
//' @return Log-likelihood value.
//'
//' @noRd
double cxx_logLik(const arma::cube & x, const arma::mat & Psi, const arma::mat & Sig){
    unsigned int r0 = x.n_rows ;
    unsigned int c0 = x.n_cols ;
    unsigned int n0 = x.n_slices ;

    // `arma::inv_sympd(A)` throws warning if A is asymmetric
    arma::mat PsiInv = arma::inv_sympd(Psi) ;
    arma::mat SigInv = arma::inv_sympd(Sig) ;

    double logTr = 0.0 ;
    for (unsigned int i=0; i < n0; i++)
        logTr += arma::accu( (PsiInv*x.slice(i)) % (x.slice(i)*SigInv) ) ;

    // `arma::log_det(val, sign, A)` leads to `|A| = exp(val)*sign`
    double d1, d2, temp ;
    arma::log_det(d1, temp, Psi) ;
    arma::log_det(d2, temp, Sig) ;

    return -0.5*(n0*r0*c0*log(2*M_PI) + n0*c0*d1 + n0*r0*d2 + logTr) ;
}



//' @title MLE of Row and Column Covariance Matrices
//' @description \loadmathjax
//' It supposes the data are sampled independently from
//'   \mjsdeqn{N(0, \Psi, \Sigma)}.
//'
//' @param x Array.
//' @param maxIter Maximal step of iterations.
//' @param tol Tolerance.
//'
//' @return `list(Psi, PsiInv, Sig, SigInv, logLik)`.
//'
//' @noRd
List cxx_mle(const arma::cube & x, unsigned int maxIter=100, double tol=1.0e-8){
    unsigned int r0 = x.n_rows ;
    unsigned int c0 = x.n_cols ;
    unsigned int n0 = x.n_slices ;

    arma::mat Psi = arma::eye(r0, r0), PsiInv = arma::eye(r0, r0) ;
    arma::mat Sig = arma::eye(c0, c0), SigInv = arma::eye(c0, c0) ;

    double logLik = cxx_logLik(x, Psi, Sig) ;
    double logLikOld = logLik ;
    double err = tol + 1.0 ;
    unsigned int i=0, k=0 ;
    
    while (err > tol && k <= maxIter) {
        Psi.zeros() ;
        for (i=0; i < n0; i++)
            Psi += x.slice(i) * SigInv * x.slice(i).t() / (n0*c0) ;
        Psi /= Psi(0, 0) ;
        PsiInv = arma::inv_sympd(Psi) ;

        Sig.zeros() ;
        for (i=0; i < n0; i++)
            Sig += x.slice(i).t() * PsiInv * x.slice(i) / (n0*r0) ;
        SigInv = arma::inv_sympd(Sig) ;

        logLikOld = logLik ;
        logLik = cxx_logLik(x, Psi, Sig) ;
        err = fabs(logLik - logLikOld) ;
        k++ ;
    }

    return Rcpp::List::create(
        ::Named("Psi")=Psi, ::Named("PsiInv") = PsiInv, 
        ::Named("Sig")=Sig, ::Named("SigInv") = SigInv,
        ::Named("logLik")=logLik
    ) ;
}



//' @title Mean Matrices
//' @description
//' If `isTRUE(flag(i, j))`, `M1(i, j)` and `M2(i, j)` are estimated as
//' the pooled mean value.
//'
//' @param x1 See [get_suppSet()].
//' @param x2 See [get_suppSet()].
//' @param flag Logical matrix returned from [get_suppSet()].
//'
//' @return `list(M1, M2)`.
//'
//' @noRd
// [[Rcpp::export]]
List cxx_mean(const arma::cube& x1, const arma::cube& x2, const LogicalMatrix& flag){
    unsigned int r0 = x1.n_rows ;
    unsigned int c0 = x1.n_cols ;
    double n1 = x1.n_slices ;
    double n2 = x2.n_slices ;
    double pi01 = n1/(n1+n2) ;
    double pi02 = n2/(n1+n2) ;

    unsigned int iRow=0, iCol=0 ;
    arma::mat mu1=mean(x1, 2), M1=mu1 ;
    arma::mat mu2=mean(x2, 2), M2=mu2 ;

    for (iRow=0; iRow < r0; iRow++) { for (iCol=0; iCol < c0; iCol++) {
        if (flag(iRow, iCol)) {
            M1(iRow, iCol) = pi01*mu1(iRow, iCol) + pi02*mu2(iRow, iCol) ;
            M2(iRow, iCol) = pi01*mu1(iRow, iCol) + pi02*mu2(iRow, iCol) ;
        }
    }}

    return List::create(::Named("M1")=M1, ::Named("M2")=M2) ;
}



//' @title Center Data
//'
//' @param x1 See [get_suppSet()].
//' @param x2 See [get_suppSet()].
//' @param flag The result returned from [get_suppSet()].
//'
//' @return Array.
//'
//' @noRd
arma::cube cxx_center_data(const arma::cube & x1, const arma::cube & x2, const LogicalMatrix & flag){
    unsigned int r0 = x1.n_rows ;
    unsigned int c0 = x1.n_cols ;
    unsigned int n1 = x1.n_slices ;
    unsigned int n2 = x2.n_slices ;

    unsigned int i=0 ;
    List M = cxx_mean(x1, x2, flag) ;
    arma::mat M1 = M["M1"] ;
    arma::mat M2 = M["M2"] ;
    arma::cube dat = arma::zeros(r0, c0, n1+n2) ;

    for (i=0; i < n1; i++)
        dat.slice(i) = x1.slice(i) - M1 ;

    for (i=0; i < n2; i++)
        dat.slice(n1+i) = x2.slice(i) - M2 ;

    return dat ;
}



//' @title MLE of Row and Column Precision Matrices
//'
//' @param x1  See [get_suppSet()].
//' @param x2  See [get_suppSet()].
//' @param flag Matrix returned from [get_suppSet()].
//'
//' @return See the value of [cxx_mle()].
//'
//' @noRd
// [[Rcpp::export]]
List cxx_prec (arma::cube x1, arma::cube x2, LogicalMatrix flag) {
    return cxx_mle( cxx_center_data(x1, x2, flag) ) ;
}
