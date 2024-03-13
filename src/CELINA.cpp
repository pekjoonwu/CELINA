// Function to estimate the null model - Celina 
#include <fstream>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <bitset>


using namespace std;
using namespace arma;
using namespace Rcpp;

#define ARMA_DONT_PRINT_ERRORS

// Do eigen decomposition of symetric matrix
// [[Rcpp::export]]
SEXP SysMatEigen(SEXP Min) {
    try {
        arma::fmat M = as<fmat>(Min);
        arma::fvec eigval = zeros<fvec>(M.n_rows);
        arma::fmat eigvec = zeros<fmat>(size(M));
        eig_sym(eigval, eigvec, M, "dc");
        const uvec idx = find(eigval < 1e-8);
        arma::fvec tmp_value = ones<arma::fvec>(idx.n_elem);
        eigval.elem(idx) = tmp_value * 1e-8;
        arma::fmat M1 = eigvec.each_row() % eigval.t();
        M = M1 * eigvec.t();
        // return values
        return List::create(Named("eigval") = eigval, Named("eigvec") = eigvec, Named("relatedness") = M);
        //return List::create(Named("eigval") = eigval, Named("eigvec") = eigvec);
    }
    catch (std::exception &ex)
    {
        forward_exception_to_r(ex);
    }
    catch (...)
    {
        ::Rf_error("C++ exception (unknown reason)...");
    }
    return R_NilValue;
} // end func

//*************************************************************************************//
//                             AVERAGE INFORMATION METHOD                              //
//*************************************************************************************//
// Average information algorithm for estimating variance components 
// for the null model
// Yin Outcome vector
// Xin Cell type identity vector or cell type proportion vector
// Win covariates matrix
// numVarin number of variance parameter
// tauin Initial value for variance component
// fixtauin Variance component to be optimized
// tolin Tolerance
// [[Rcpp::export]]
SEXP LmmAICovarites(SEXP Yin, SEXP Xin, SEXP Win, SEXP numVarin, SEXP tauin, 
           SEXP fixtauin, SEXP tolin) { 
  try {
    // Initialization
    vec Y = as<vec>(Yin); 
    vec X = as<vec>(Xin); 
    mat W = as<mat>(Win);
    int numVar = Rcpp::as<int>(numVarin);
    vec tau = as<vec>(tauin); 
    const uvec fixtau = as<uvec>(fixtauin); // indicator vector for which parameters to be estimated 
    int num_cov2 = sum(fixtau == 0); // number of variance parameter to be estimated
    const double tol = Rcpp::as<double>(tolin);
    uvec ZERO = (tau < tol);
    size_t n_cells = X.n_elem;// number of locations
    cube PHI(n_cells, 1, numVar);
    
    // to get the H inverse
    vec Hinv(n_cells, fill::zeros);
    vec one_vec = ones<vec>(n_cells);
    Hinv = tau(1) * square(X);
    Hinv += tau(2) * one_vec;
    Hinv = 1.0/(Hinv + 1e-5);
    PHI.slice(0) = square(X);
    PHI.slice(1) = one_vec;
    
    // Calculate some identities
    vec HinvY = Hinv % Y;
    mat HinvW = W.each_col() % Hinv;
    mat WtHinvW = W.t() * HinvW;
    mat WtHinvW_inv = inv_sympd(WtHinvW);
    // vec P_diag = Hinv - (HinvW % HinvW) * WtHinvW_inv;
    vec alpha = HinvW.t() * Y;
    alpha = WtHinvW_inv * alpha;
    vec PY = HinvY - HinvW * alpha;
    double tracePA;
    // vec tmp_identites = ;
    
    // AI update
    if (num_cov2 > 0) {
      const uvec idxtau = find(fixtau == 0); 
      mat AImat(num_cov2, num_cov2); // average information matrix
      vec score(num_cov2); // score vector
      for (size_t i = 0; i < num_cov2; i++) {
        vec APY = PHI.slice(idxtau(i) - 1) % PY;
        vec PAPY = Hinv % APY - HinvW * (WtHinvW_inv * (HinvW.t() * APY));
        tracePA = sum(Hinv % PHI.slice(idxtau(i) - 1)) - 
          trace(WtHinvW_inv * (HinvW.t() * (HinvW.each_col() % PHI.slice(idxtau(i) - 1) )));
        score(i) = dot(Y, PAPY) - tracePA;
        
        for (size_t j = 0; j <= i; j++) {
          AImat(i, j) = dot(PY, PHI.slice(idxtau(j) - 1) % PAPY);
          if (j != i) {
            AImat(j, i) = AImat(i, j);
          }
        } //end for j
      } // end for i
      
      mat AImat_inv = pinv(AImat);
      vec Dtau = AImat_inv * score;
      vec tau0 = tau;
      
      //NR update
      tau.elem(idxtau) = tau0.elem(idxtau) + Dtau; 
      tau.elem(find(ZERO % (tau < tol))).zeros();
      double step = 1.0;
      while (any(tau < 0.0)) {
        step *= 0.5;
        tau.elem(idxtau) = tau0.elem(idxtau) + step * Dtau; // NR update
        tau.elem(find(ZERO % (tau < tol))).zeros();
      }
      tau.elem(find(tau < tol)).zeros();
    } // end fi
    // return values
    return List::create(Named("tau") = tau, Named("cov") = WtHinvW_inv,
                        Named("alpha") = alpha, Named("Py") = PY);
  }
  catch (std::exception &ex)
  {
    forward_exception_to_r(ex);
  }
  catch (...)
  {
    ::Rf_error("C++ exception (unknown reason)...");
  }
  return R_NilValue;
} // end func


// Compute the testing quantities with covariates, float format
// yin Working vector
// Pyin The vector P*y
// Win Covariate matrix, including the intercept and the target cell type proportion
// Xin target cell type proportions
// XKXin Kernel matrix to be tested multiplies with two cell type proportion diagonal matrix
// Din Weight for each gene
// tauin variance component estimated values
// Prightin some matrix part in the P matrix
// [[Rcpp::export]]
List TestStatfast(arma::vec yin, arma::vec Pyin, arma::mat Win, arma::vec Xin, 
                  arma::mat& XKXin, arma::vec tauin, 
                  arma::mat& Prightin) {
  try {
    // Initialization
    size_t n_cells = yin.n_elem;
    
    // to get the H inverse
    vec Hinv(n_cells, fill::zeros);
    vec one_vec = ones<vec>(n_cells);
    Hinv = tauin(1) * square(Xin);
    Hinv += tauin(2) * one_vec;
    Hinv = 1.0/(Hinv + 1e-5);

    // Calculate target PK - faster implementation
    mat HinvK = XKXin.each_col() % Hinv;
    mat WtHinvK = Win.t() * HinvK;
    mat PK = HinvK - Prightin * WtHinvK;
    
    // Calculate the mean and variances of the test statistics and the statistics
    double trace_PK = trace(PK);
    double trace_PKPK = accu(PK % PK.t());
    vec PKPY = PK * Pyin;
    double S0 = 0.5 * dot(yin, PKPY);
    double E = trace_PK / 2.0;
    double I_raw = 0.5 * trace_PKPK;

    // return the results
    return List::create(Named("S0") = S0, Named("E") = E, 
                        Named("I_raw") = I_raw);
  } // end try
  catch (std::exception &ex)
  {
    forward_exception_to_r(ex);
  }
  catch (...)
  {
    ::Rf_error("C++ exception (unknown reason)...");
  }
  return R_NilValue;
} // end func

