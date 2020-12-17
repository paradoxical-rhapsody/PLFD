// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' @title Log-Likelihood of Matrix-Variate Normal Data
//' @description It supposes that the data are independently
//'   sampled from `N(0, Psi, Sigma)`.
//'
//' @param x Array of *r x c x n*.
//' @param Psi Row convariance matrix.
//' @param Sig Column covariance matrix.
//'
//' @noRd
// [[Rcpp::export]]
double cxx_matNormal_lik (const arma::cube& x, const arma::mat& Psi, const arma::mat& Sig) {
	unsigned int r0 = x.n_rows ;
	unsigned int c0 = x.n_cols ;
	unsigned int n0 = x.n_slices ;

	// `arma::inv_sympd(A)` throws warning if A is asymmetric even because of round-off.
	// `arma::inv(A)` is for general square matrix.
	arma::mat PsiInv = arma::inv_sympd(Psi) ;
	arma::mat SigInv = arma::inv_sympd(Sig) ;

	double logTr = 0.0 ;
	for (unsigned int i=0; i < n0; i++) {
		logTr += arma::accu((PsiInv * x.slice(i)) % (x.slice(i) * SigInv)) ;
	}

	double d1, d2, temp ; // log_det of Psi and Sig
	arma::log_det(d1, temp, Psi) ; // arma::log_det(val, sign, A) -> |A| = exp(val)*sign
	arma::log_det(d2, temp, Sig) ;

	return -0.5*(n0*r0*c0*log(2*PI) + n0*c0*d1 + n0*r0*d2 + logTr) ;
}

// MLE of Psi/Sigma/PsiInv/SigmaInv of N_{rxc}(0, Psi, Sigma).
// Return list(logLik, Psi, Sig, PsiInv, SigInv)
List cxx_matNormal_mle (arma::cube x, unsigned int maxIter=100, double tol=1.0e-8) {
	unsigned int r0 = x.n_rows ;
	unsigned int c0 = x.n_cols ;
	unsigned int n0 = x.n_slices ;

	arma::mat Psi = arma::eye(r0, r0), PsiInv = arma::eye(r0, r0) ;
	arma::mat Sig = arma::eye(c0, c0), SigInv = arma::eye(c0, c0) ;

	double logLik = cxx_matNormal_lik(x, Psi, Sig) ;
	double logLikOld = logLik ;
	double err = tol + 1.0 ;
	unsigned int k = 0 ;
	double normConst = 1.0 ;
	while (err > tol && k <= maxIter) {
		Psi.zeros() ;	// reset as full zero for accumulation
		for (unsigned int i=0; i < n0; i++) {
			Psi += x.slice(i) * SigInv * x.slice(i).t() / (n0*c0) ;
		}
		normConst = Psi(0, 0) ;
		Psi /= normConst ;
		PsiInv = arma::inv(Psi) ;

		Sig.zeros() ;
		for (unsigned int i=0; i < n0; i++) {
			Sig += x.slice(i).t() * PsiInv * x.slice(i) / (n0*r0) ;
		}
		SigInv = arma::inv(Sig) ;

		logLikOld = logLik ;
		logLik = cxx_matNormal_lik(x, Psi, Sig) ;
		err = fabs(logLik - logLikOld) ;
		k++ ;
	}

	return Rcpp::List::create(
		::Named("logLik")=logLik,
		::Named("Psi")=Psi, ::Named("PsiInv") = PsiInv, 
		::Named("Sig")=Sig, ::Named("SigInv") = SigInv
	) ;
}
 
// If `(flag)_{ij} = FALSE`, (M1)_{ij} and M2{ij} are within group mean 
// value; otherwize they are estimated as pooled mean value. Return 
// numeric array of `r x c x 2`.
// [[Rcpp::export]]
arma::cube cxx_mean (arma::cube x1, arma::cube x2, LogicalMatrix flag) {
	unsigned int r0 = x1.n_rows ;
	unsigned int c0 = x1.n_cols ;
	unsigned int n1 = x1.n_slices ;
	unsigned int n2 = x2.n_slices ;
	double pi01 = (double) n1/(n1+n2) ;
	double pi02 = (double) n2/(n1+n2) ;
	arma::mat mu1 = mean(x1, 2) ;
	arma::mat mu2 = mean(x2, 2) ;

	arma::cube M = arma::zeros(r0, c0, 2) ;
	M.slice(0) = mu1 ;
	M.slice(1) = mu2 ;
	for (unsigned int iRow=0; iRow < r0; iRow++) {
		for (unsigned int iCol=0; iCol < c0; iCol++) {
			if (flag(iRow, iCol)) {
				M(iRow, iCol, 0) = pi01*mu1(iRow, iCol) + pi02*mu2(iRow, iCol) ;
				M(iRow, iCol, 1) = pi01*mu1(iRow, iCol) + pi02*mu2(iRow, iCol) ;
			}
		}
	}

	return M ;
}

// Centeralize samples by substracting samples by their mean matrice and bind 
// them together. `flag` is logical matrix of `r x c` and is used to estimate
// the mean matrice. Return Numeric array of `r x c x (n1+n2)`.
arma::cube cxx_centralize_samples (arma::cube x1, arma::cube x2, LogicalMatrix flag) {
	unsigned int r0 = x1.n_rows ;
	unsigned int c0 = x1.n_cols ;
	unsigned int n1 = x1.n_slices ;
	unsigned int n2 = x2.n_slices ;
	arma::cube M = cxx_mean(x1, x2, flag) ;

	arma::cube dat = arma::zeros(r0, c0, n1+n2) ;
	for (unsigned int i=0; i < n1; i++) {
		dat.slice(i) = x1.slice(i) - M.slice(0) ;
	}
	for (unsigned int i=0; i < n2; i++) {
		dat.slice(n1+i) = x2.slice(i) - M.slice(1) ;
	}

	return dat ;
}

// Use the centralized samples from `x1` and `x2` to obtain the MLE of the row 
// and the column covariance matrices and their inverse matrices.
// [[Rcpp::export]]
List cxx_prec (arma::cube x1, arma::cube x2, LogicalMatrix flag) {
    return cxx_matNormal_mle(cxx_centralize_samples(x1, x2, flag)) ;
}
