

##' The main function for fitting multivariate bcf
##'
##' 
##' 
##' @param X_con The control variables used in the mu trees
##' @param y The outcome variable/s 
##' @param Z The treatment indicator
##' @param X_mod The effect moderators used in the tau trees
##' @param X_con_test The control variables used in the mu trees (test data)
##' @param X_mod_test The effect moderators used in the tau trees (test data)
##' @param alpha The alpha parameter in the tree prior for mu trees
##' @param beta The beta parameter in the tree prior for mu trees
##' @param alpha_tau The alpha parameter in the tree prior for tau trees
##' @param beta_tau The beta parameter in the tree prior for tau trees
##' @param sigma_mu The prior for the terminal node parameters of mu trees
##' @param sigma_tau The prior for the terminal node parameter of tau trees
##' @param v_0 The treatment indicator
##' @param sigma_0 The effect moderators used in the tau trees
##' @param n_iter The number of MCMC iterations
##' @param n_tree The number of mu trees
##' @param n_tree_tau The number of tau trees
##' @param min_nodesize Moves resulting in nodes with fewer observations than this are rejected
##' @return A list of model outputs
##' 
##' @export
run_mvbcf<-function(X_con,
                    y,
                    Z,
                    X_mod,
                    X_con_test,
                    X_mod_test,
                    alpha,
                    beta,
                    alpha_tau,
                    beta_tau,
                    sigma_mu,
                    sigma_tau,
                    v_0,
                    sigma_0,
                    n_iter,
                    n_tree,
                    n_tree_tau,
                    min_nodesize)
{
  fast_bart(X_con,
            y,
            Z,
            X_mod,
            X_con_test,
            X_mod_test,
            alpha,
            beta,
            alpha_tau,
            beta_tau,
            sigma_mu,
            sigma_tau,
            v_0,
            sigma_0,
            n_iter,
            n_tree,
            n_tree_tau,
            min_nodesize)
}