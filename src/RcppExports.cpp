// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fast_bart
List fast_bart(arma::mat X_con, arma::mat y, arma::mat Z, arma::mat X_mod, arma::mat X_con_test, arma::mat X_mod_test, double alpha, double beta, double alpha_tau, double beta_tau, arma::mat sigma_mu, arma::mat sigma_tau, int v_0, arma::mat sigma_0, int n_iter, int n_tree, int n_tree_tau, int min_nodesize);
RcppExport SEXP _mvbcf_fast_bart(SEXP X_conSEXP, SEXP ySEXP, SEXP ZSEXP, SEXP X_modSEXP, SEXP X_con_testSEXP, SEXP X_mod_testSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP alpha_tauSEXP, SEXP beta_tauSEXP, SEXP sigma_muSEXP, SEXP sigma_tauSEXP, SEXP v_0SEXP, SEXP sigma_0SEXP, SEXP n_iterSEXP, SEXP n_treeSEXP, SEXP n_tree_tauSEXP, SEXP min_nodesizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X_con(X_conSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_mod(X_modSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_con_test(X_con_testSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_mod_test(X_mod_testSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_tau(alpha_tauSEXP);
    Rcpp::traits::input_parameter< double >::type beta_tau(beta_tauSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_mu(sigma_muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_tau(sigma_tauSEXP);
    Rcpp::traits::input_parameter< int >::type v_0(v_0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_0(sigma_0SEXP);
    Rcpp::traits::input_parameter< int >::type n_iter(n_iterSEXP);
    Rcpp::traits::input_parameter< int >::type n_tree(n_treeSEXP);
    Rcpp::traits::input_parameter< int >::type n_tree_tau(n_tree_tauSEXP);
    Rcpp::traits::input_parameter< int >::type min_nodesize(min_nodesizeSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_bart(X_con, y, Z, X_mod, X_con_test, X_mod_test, alpha, beta, alpha_tau, beta_tau, sigma_mu, sigma_tau, v_0, sigma_0, n_iter, n_tree, n_tree_tau, min_nodesize));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _mvbcf_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _mvbcf_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _mvbcf_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _mvbcf_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mvbcf_fast_bart", (DL_FUNC) &_mvbcf_fast_bart, 18},
    {"_mvbcf_rcpparma_hello_world", (DL_FUNC) &_mvbcf_rcpparma_hello_world, 0},
    {"_mvbcf_rcpparma_outerproduct", (DL_FUNC) &_mvbcf_rcpparma_outerproduct, 1},
    {"_mvbcf_rcpparma_innerproduct", (DL_FUNC) &_mvbcf_rcpparma_innerproduct, 1},
    {"_mvbcf_rcpparma_bothproducts", (DL_FUNC) &_mvbcf_rcpparma_bothproducts, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_mvbcf(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
