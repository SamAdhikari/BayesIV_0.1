#' Gibbs_Reg: Gibbs sampler function to sample slope coefficients in
#' linear regression with normal noise
#'
#' @param Y : dependent variable
#' @param X : covariates
#' @param sigmasq_prior : variance of normal prior
#' @param sigmasq_prior_alpha : variance of the normal prior for \code{alpha}
#' @param sigmasq : variance of \code{Y}
#' @param PP : number of covariates
#' @return A draw from the posterior distribution 
#' @export


Gibbs_Reg = function(Y,X,sigmasq_prior,sigmasq_prior_alpha,sigmasq,PP)
{
  sigma_prior = diag(c(sigmasq_prior_alpha,rep(sigmasq_prior,(PP-1))),PP)
  precision_prior = solve(sigma_prior)
  precision_n_inv = solve(t(X)%*%X + precision_prior+diag(0.01,PP))
#  print(precision_n_inv)
  postMean = (precision_n_inv)%*%t(X)%*%Y
  postVar = diag(sigmasq,PP) * precision_n_inv
  return(rnorm(PP,postMean,sqrt(diag(postVar))))  
}
