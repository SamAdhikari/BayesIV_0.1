#' getATE_posterior_NDPM
#'
#' @author Sam Adhikari
#' @param fittedModel : model fits of the class 'BayesIV' 
#' @param X : a numeric matrix of covariates (excluding the column with intercept)
#' @param niter : number indicating the size of the MCMC sample
#' @param burnin: number indicating burnin of the MCMC chain
#' @param thin : number to thin the MCMC chain with
#' @return a list with posterior draws of Y1hat, Y0hat and ATE
#' @export


getATE_posterior_NDPM = function(fittedModel,X,niter,burnin=0,thin=1)
{
    Chain = seq(burnin,niter,by=thin)
    ##post processing of the chain to make sure the error has expectation 0
    mean1 = apply(fittedModel$muY1,2,mean)
    mean0 = apply(fittedModel$muY0,2,mean)
    
    Y1hat_Chain = sapply(Chain,function(x){
        X%*%fittedModel$Beta1[,x]+fittedModel$Theta[,x]*fittedModel$Alpha1[x] +
            mean1[x]})
    Y0hat_Chain = sapply(Chain,function(x){
        X%*%fittedModel$Beta0[,x]+fittedModel$Theta[,x]*fittedModel$Alpha0[x] +
            mean0[x]})
    ##ATE at each step of the chain
    ATE_Chain_N_DPM = apply(Y1hat_Chain-Y0hat_Chain,2,mean)
    return(list(Y1hat = Y1hat_Chain,Y0hat=Y0hat_Chain,ATE=ATE_Chain_N_DPM))
}



