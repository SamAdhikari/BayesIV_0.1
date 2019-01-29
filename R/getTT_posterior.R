#' getTT_posterior
#'
#' @param fittedModel : model fits of the class 'BayesIV' 
#' @param X : a numeric matrix of covariates (excluding the column with intercept)
#' @param niter : number indicating the size of the MCMC sample
#' @param burnin: number indicating burnin of the MCMC chain
#' @param thin : number to thin the MCMC chain with
#' @return posterior chain of ATT
#' @export 


getTT_posterior = function(fittedModel,X,niter,burnin=0,thin=1)
{
  Chain = seq(burnin,niter,by=thin)
  ##Treatment effect on the treated
  D1hatChain = sapply(Chain,function(x){
    postvals = X%*%fittedModel$BetaT[,x] +
      fittedModel$Theta[,x]*fittedModel$AlphaD[x] +
      Z*fittedModel$Gamma[x] 
    return(pnorm(postvals)) ## normal CDF
  })
  MeanD1hatChain = apply(D1hatChain,2,mean)
  sumD1hat = mean(D1hatChain)
  Y1hat_Chain = sapply(Chain,function(x){
    X%*%fittedModel$Beta1[,x]+
      fittedModel$Theta[,x]*fittedModel$Alpha1[x]})
  Y0hat_Chain = sapply(Chain,function(x){
    X%*%fittedModel$Beta0[,x]+
      fittedModel$Theta[,x]*fittedModel$Alpha0[x]})
  
  NumeratorChain =apply((Y1hat_Chain-Y0hat_Chain),2,mean) 
  
  NumeratorSum = NumeratorChain*MeanD1hatChain
  TT_treatedChain = NumeratorSum/sumD1hat
  return(TT_treatedChain)
}