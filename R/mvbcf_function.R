


#/#' A helper function used in the mvbcf package
#/#'
#/#' 
#/#' 
#/#' @param mat matrix to check if it is a valid covariance matrix
#/#' @return A boolean true/false
#/#' 
is_valid_cov_matrix <- function(mat) {
  # Check if matrix is square
  if (!is.matrix(mat) || nrow(mat) != ncol(mat)) {
    return(FALSE)
  }
  
  # Check for symmetry
  if (!all(mat == t(mat))) {
    return(FALSE)
  }
  
  # Check for non-null values (no NA or NaN)
  if (any(is.na(mat)) || any(is.nan(mat))) {
    return(FALSE)
  }
  
  # Check for positive semi-definiteness (eigenvalues must be non-negative)
  eigenvalues <- eigen(mat)$values
  if (any(eigenvalues < -1e-10)) {  # Tolerance for numerical precision
    return(FALSE)
  }
  
  return(TRUE)
}




#/#' A helper function used in the mvbcf package
#/#'
#/#' 
#/#' 
#/#' @param mat matrix to check if it is a valid covariance matrix
#/#' @return A boolean true/false
#/#' 
is_valid_scale_matrix <- function(mat) {
  # Check if the input is a matrix and square
  if (!is.matrix(mat) || nrow(mat) != ncol(mat)) {
    return(FALSE)
  }
  
  # Check for symmetry (matrix equals its transpose)
  if (!all(mat == t(mat))) {
    return(FALSE)
  }
  
  # Check for positive definiteness (all eigenvalues must be positive)
  eigenvalues <- eigen(mat)$values
  if (any(eigenvalues <= 0)) {
    return(FALSE)
  }
  
  return(TRUE)
}


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
##' @param v_0 Prior degrees of freedom for inverse-wishart distribution
##' @param sigma_0 The scale matrix in the prior for inverse-wishart distribution
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
                    X_con_test = X_con,
                    X_mod_test = X_mod,
                    alpha = 0.95,
                    beta = 2,
                    alpha_tau = 0.25,
                    beta_tau = 3,
                    sigma_mu = diag((1)^2/n_tree, ncol(y)),
                    sigma_tau = diag((1)^2/n_tree_tau, ncol(y)),
                    v_0 = ncol(y)+2,
                    sigma_0 = diag(1, ncol(y)),
                    n_iter = 1000,
                    n_tree = 50,
                    n_tree_tau = 20,
                    min_nodesize = 1)
{
  
  stopifnot("X_con is not a numeric matrix" =
              is.matrix(X_con) & is.numeric(X_con))
  
  stopifnot("y is not a numeric matrix" =
              is.matrix(y) & is.numeric(y))
  
  stopifnot("Z is not a numeric vector of 1s and 0s" =
              is.numeric(Z) & (sum(Z==0) + sum(Z==1) == length(Z)))
  
  stopifnot("X_mod is not a numeric matrix" =
              is.matrix(X_mod) & is.numeric(X_mod))
  
  stopifnot("X_con_test is not a numeric matrix" =
              is.matrix(X_con_test) & is.numeric(X_con_test))
  
  stopifnot("X_mod_test is not a numeric matrix" =
              is.matrix(X_mod_test) & is.numeric(X_mod_test))
  
  stopifnot("y does not have same number of rows as X_con" =
              dim(X_con)[1] == dim(y)[1])
  
  stopifnot("Z does not have same length as number of rows in X_con" =
              dim(X_con)[1] == length(Z))
  
  stopifnot("X_mod does not have same number of rows as X_con" =
              dim(X_con)[1] == dim(X_mod)[1])
  
  stopifnot("X_con_test does not have same number of rows as X_mod_test" =
              dim(X_con_test)[1] == dim(X_mod_test)[1])
  
  stopifnot("X_con_test does not have same number of cols as X_con" =
              dim(X_con)[2] == dim(X_con_test)[2])
  
  stopifnot("X_mod_test does not have same number of cols as X_mod" =
              dim(X_con)[2] == dim(X_con_test)[2])
  
  stopifnot("alpha must be a single numeric value" =
              length(alpha)==1 & is.numeric(alpha))
  
  stopifnot("beta must be a single numeric value" =
              length(beta)==1 & is.numeric(beta))
  
  stopifnot("alpha_tau must be a single numeric value" =
              length(alpha_tau)==1 & is.numeric(alpha_tau))
  
  stopifnot("beta_tau must be a single numeric value" =
              length(beta_tau)==1 & is.numeric(beta_tau))
  
  stopifnot("sigma_mu must be a valid covariance matrix with number of rows and columns equal to the number of columns in y" =
              is_valid_cov_matrix(sigma_mu) & ncol(sigma_mu) == ncol(y))
  
  stopifnot("sigma_tau must be a valid covariance matrix with number of rows and columns equal to the number of columns in y" =
              is_valid_cov_matrix(sigma_tau) & ncol(sigma_tau) == ncol(y))
  
  stopifnot("v_0 must be a whole number greater than or equal to 1" =
              round(v_0)==v_0 & v_0>=1)
  
  stopifnot("sigma_0 must be a valid scale matrix with number of rows and columns equal to the number of columns in y" =
              is_valid_scale_matrix(sigma_0) & ncol(sigma_0) == ncol(y))
  
  stopifnot("n_iter must be a whole number greater than or equal to 1" =
              round(n_iter)==n_iter & n_iter>=1)
  
  stopifnot("n_tree must be a whole number greater than or equal to 1" =
              round(n_tree)==n_tree & n_tree>=1)
  
  stopifnot("n_tree_tau must be a whole number greater than or equal to 1" =
              round(n_tree_tau)==n_tree_tau & n_tree_tau>=1)
  
  stopifnot("min_nodesize must be a whole number greater than or equal to 1" =
              round(min_nodesize)==min_nodesize & min_nodesize>=1)
  
  Z_mat<-matrix(rep(Z, ncol(y)), ncol=ncol(y), byrow=F) 
  
  fast_bart(X_con,
            y,
            Z_mat,
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