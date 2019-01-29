#' getTT_posterior_NDPM
#'
#' @author Sam Adhikari
#' @param fittedModel : model fits of the class 'BayesIV' 
#' @param X0c : a numeric matrix of covariates (excluding the column with intercept)
#' @param XTr : a numeric matrix of covariates
#' @param Z: A vector of instrumental variable 
#' @param niter : number indicating the size of the MCMC sample
#' @param burnin: number indicating burnin of the MCMC chain
#' @param thin : number to thin the MCMC chain with
#' @return posterior chain of ATT
#' @export


getTT_posterior_NDPM = function(fittedModel,XOc,XTr,Z,niter,burnin=0,thin=1)
{
    Chain = seq(burnin,niter,by=thin)
    ##Treatment effect on the treated
    D1hatChain = sapply(Chain,function(x){
        postvals = XTr%*%fittedModel$BetaT[,x] +
            fittedModel$Theta[,x]*fittedModel$AlphaD[x] +
            Z*fittedModel$Gamma[x] 
        return(pnorm(postvals)) ## normal CDF
    })
    MeanD1hatChain = apply(D1hatChain,2,mean)
    sumD1hat = mean(D1hatChain)
    
    ##post processing of the chain to make sure the error has expectation 0
    mean1 = apply(fittedModel$muY1,2,mean)
    mean0 = apply(fittedModel$muY0,2,mean)
    
    Y1hat_Chain = sapply(Chain,function(x){
        XOc%*%fittedModel$Beta1[,x]+fittedModel$Theta[,x]*fittedModel$Alpha1[x] +
            mean1[x]})
    Y0hat_Chain = sapply(Chain,function(x){
        XOc%*%fittedModel$Beta0[,x]+fittedModel$Theta[,x]*fittedModel$Alpha0[x] +
            mean0[x]})
    
    NumeratorChain =apply((Y1hat_Chain-Y0hat_Chain),2,mean) 
    
    NumeratorSum = NumeratorChain*MeanD1hatChain
    TT_treatedChain = NumeratorSum/sumD1hat
    return(TT_treatedChain)
}

