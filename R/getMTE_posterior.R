#' getMTE_posterior
#' 
#' @author Sam Adhikari
#' @param fittedModel : model fits of the class 'BayesIV' 
#' @param u: margin
#' @param niter : number indicating the size of the MCMC sample
#' @param burnin: number indicating burnin of the MCMC chain
#' @param thin : number to thin the MCMC chain with
#' @return a list with posterior draws of MTE
#' @export

getweights_MTE = function(fittedModel,u,niter,burnin=0,thin=1){
  Chain = seq(burnin,niter,by=thin)
  numerator = sapply(Chain,function(x){
    mean(dnorm(u-fittedModel$Theta[,x]*fittedModel$AlphaD[x]))})
  ##normal density function
  return(numerator/sum(numerator))
}


getMTE_posterior = function(fittedModel,u,X,ATE_Chain,
                            niter,burnin=0,thin=1)
{
  Chain = seq(burnin,niter,by=thin)
  weights = getweights_MTE(fittedModel=fittedModel,
                           u=u,niter=niter,burnin=burnin,thin=thin)
  MTE_u = weights*ATE_Chain
  return(MTE_u)
}
