
##' A helper function for getting tau predictions
##'
##' @param fitted_model The fitted mvbcf model to get predictions from
##' @param outcome_number The outcome variable to get predictions for
##' @param train_or_test Indicate if predictions should be for the "train" or "test" data
##' @return A vector of individual treatment effect estimates for the specified outcome
##' 
##' @export
get_tau_preds<-function(fitted_model,
                     outcome_number,
                     train_or_test)
{
  stopifnot("fitted_model must be a fitted mvbcf model output" =
              is.list(fitted_model) & names(fitted_model)==c("predictions", 
                                                             "predictions_tau", 
                                                             "sigmas", 
                                                             "predictions_test", 
                                                             "predictions_tau_test"))
  
  stopifnot("outcome_number is not in the valid range 1 to ncol(y)" =
              outcome_number %in% c(1:dim(fitted_model$sigmas)[1]))
  
  stopifnot("train_or_test must be 'train' or 'test'" =
              train_or_test %in% c("train", "test"))
  
  if(train_or_test=="train")
  {
    return(rowMeans(fitted_model$predictions_tau[,outcome_number,]))
  }
  if(train_or_test=="test")
  {
    return(rowMeans(fitted_model$predictions_tau_test[,outcome_number,]))
  }
}

##' A helper function for getting the posterior of the ATE
##'
##' @param fitted_model The fitted mvbcf model to get the ATE posterior from
##' @param outcome_number The outcome variable to the ATE posterior for
##' @param train_or_test Indicate if the posterior should be for the "train" or "test" data
##' @return The posterior of the ATE for the specified outcome variable
##' 
##' @export
get_ate_post<-function(fitted_model,
                       outcome_number,
                       train_or_test)
{
  stopifnot("fitted_model must be a fitted mvbcf model output" =
              is.list(fitted_model) & names(fitted_model)==c("predictions", 
                                                             "predictions_tau", 
                                                             "sigmas", 
                                                             "predictions_test", 
                                                             "predictions_tau_test"))
  
  stopifnot("outcome_number is not in the valid range 1 to ncol(y)" =
              outcome_number %in% c(1:dim(fitted_model$sigmas)[1]))
  
  stopifnot("train_or_test must be 'train' or 'test'" =
              train_or_test %in% c("train", "test"))
  
  if(train_or_test=="train")
  {
    return(colMeans(fitted_model$predictions_tau[,outcome_number,]))
  }
  if(train_or_test=="test")
  {
    return(colMeans(fitted_model$predictions_tau_test[,outcome_number,]))
  }
}

##' A helper function for getting the posterior of specified sigma values. If
##' the first_outcome and second_outcome are the same, then the residual variance
##' for that outcome is returned. Otherwise the residual covariance for the
##' first_outcome and second_outcome is returned.
##'
##' @param fitted_model The fitted mvbcf model to get the posterior from
##' @param first_outcome The first outcome variable
##' @param second_outcome The second outcome variable
##' @return The posterior of the specified residual variance or covariance
##' 
##' @export
get_sigma_post<-function(fitted_model,
                        first_outcome,
                        second_outcome)
{
  stopifnot("fitted_model must be a fitted mvbcf model output" =
              is.list(fitted_model) & names(fitted_model)==c("predictions", 
                                                             "predictions_tau", 
                                                             "sigmas", 
                                                             "predictions_test", 
                                                             "predictions_tau_test"))
  
  stopifnot("first_outcome is not in the valid range 1 to ncol(y)" =
              first_outcome %in% c(1:dim(fitted_model$sigmas)[1]))
  
  stopifnot("second_outcome is not in the valid range 1 to ncol(y)" =
              second_outcome %in% c(1:dim(fitted_model$sigmas)[1]))
  
  return(fitted_model$sigmas[first_outcome,second_outcome,])
  
}



##' A helper function for getting mu predictions
##'
##' @param fitted_model The fitted mvbcf model to get predictions from
##' @param outcome_number The outcome variable to get predictions for
##' @param train_or_test Indicate if predictions should be for the "train" or "test" data
##' @return A vector of individual mu (prognostic effect) estimates for the specified outcome
##' 
##' @export
get_mu_preds<-function(fitted_model,
                     outcome_number,
                     train_or_test)
{
  stopifnot("fitted_model must be a fitted mvbcf model output" =
              is.list(fitted_model) & names(fitted_model)==c("predictions", 
                                                             "predictions_tau", 
                                                             "sigmas", 
                                                             "predictions_test", 
                                                             "predictions_tau_test"))
  
  stopifnot("outcome_number is not in the valid range 1 to ncol(y)" =
              outcome_number %in% c(1:dim(fitted_model$sigmas)[1]))
  
  stopifnot("train_or_test must be 'train' or 'test'" =
              train_or_test %in% c("train", "test"))
  
  if(train_or_test=="train")
  {
    return(rowMeans(fitted_model$predictions[,outcome_number,]))
  }
  if(train_or_test=="test")
  {
    return(rowMeans(fitted_model$predictions_test[,outcome_number,]))
  }
}


##' A helper function for getting y predictions
##'
##' @param fitted_model The fitted mvbcf model to get predictions from
##' @param outcome_number The outcome variable to get predictions for
##' @param train_or_test Indicate if predictions should be for the "train" or "test" data
##' @param Z A binary treatment status indicator: 1 for treatment 0 for control
##' @return A vector of individual y estimates for the specified outcome
##' 
##' @export
get_y_preds<-function(fitted_model,
                      outcome_number,
                      train_or_test,
                      Z)
{
  stopifnot("fitted_model must be a fitted mvbcf model output" =
              is.list(fitted_model) & names(fitted_model)==c("predictions", 
                                                             "predictions_tau", 
                                                             "sigmas", 
                                                             "predictions_test", 
                                                             "predictions_tau_test"))
  
  stopifnot("outcome_number is not in the valid range 1 to ncol(y)" =
              outcome_number %in% c(1:dim(fitted_model$sigmas)[1]))
  
  stopifnot("train_or_test must be 'train' or 'test'" =
              train_or_test %in% c("train", "test"))
  
  
  stopifnot("Z is not a numeric vector of 1s and 0s" =
              is.numeric(Z) & (sum(Z==0) + sum(Z==1) == length(Z)))
  
  
  if(train_or_test=="train")
  {
    stopifnot("Z must have same length as number of training observations" =
                length(Z)==dim(fitted_model$predictions)[1])
    
    return(rowMeans(fitted_model$predictions[,outcome_number,]+Z*fitted_model$predictions_tau[,outcome_number,]))
  }
  if(train_or_test=="test")
  {
    stopifnot("Z must have same length as number of test observations" =
                length(Z)==dim(fitted_model$predictions_test)[1])
    
    return(rowMeans(fitted_model$predictions_test[,outcome_number,]+Z*fitted_model$predictions_tau_test[,outcome_number,]))
  }
}