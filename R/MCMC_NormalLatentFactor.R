#' mcmcRun_Normal: MCMC sampler for IV Analysis with Normal prior on the latent factors
#'
#' @author : Sam Adhikari
#' @import msm
#' @import DPpackage
#' @import Rcpp
#' @param Yobs: Input numeric vector of observed outcome of length n, can be binary or continuous. 
#' @param Tr : Binary numeric vector of treatment indicator.
#' @param X : numeric matrix of covariates of dimension n-By-p
#' @param Z : numeric vector of instrumental variable of length n
#' @param niter : number of MCMC sampler to be run
#' @return runModel : list of posterior samples for each parameter
#' @examples \dontrun{obj = mcmcRun_Normal(Yobs,Tr,X,Z,niter)}
#' @export mcmcRun_Normal


##MCMC sampler for a simple latent variable treatment effect model
#library(label.switching)
#library(msm) ##for sampling from truncated normal distribution


##Function to run the MCMC sampler
##
##
mcmcRun_Normal = function(Yobs,Tr,X,Z,niter)
{
    n = length(Yobs)
    pp = dim(X)[2]
    #set priors
    sigmasq_prior = 10000
    sigmasq_prior_alpha = 100
    sigmasq = 1
    # acc= rep(0,n)
    # tune = rep(1,n)
    # 
    S = 15
    R = 10
    
    ##Initialization
    Beta1 = rnorm(pp,0,1)
    Alpha1 = rnorm(1,0,1)
    Beta0 = rnorm(pp,0,1)
    Alpha0 = rnorm(1,0,1)
    Theta = rnorm(n,0,1)
    
    AlphaD = rnorm(1,0,1)
    BetaT = rnorm(pp,0,1)
    Gamma = rnorm(1,0,1)
    
    Trstar = Tr
    mean1 = AlphaD*Theta[which(Tr==1)] +
        Gamma*Z[which(Tr==1)] + X[which(Tr==1),]%*%BetaT
    mean0 = AlphaD*Theta[which(Tr==0)] +
        Gamma*Z[which(Tr==0)] + X[which(Tr==0),]%*%BetaT
    
    Trstar[which(Tr==0)] = rtnorm(length(which(Tr==0)),
                                  mean=mean0,
                                  sd=1, upper=0)
    Trstar[which(Tr==1)] = rtnorm(length(which(Tr==1)),
                                  mean=mean1,
                                  sd=1, lower=0)

    runModel = MCMClatentTrEff_Normal(Yobs=Yobs,X=X,Z=Z,Tr=Tr,Trstar=Trstar,
                                      n=n,pp=pp,
                                      sigmasq_prior=sigmasq_prior,
                                      sigmasq_prior_alpha = sigmasq_prior_alpha,
                                      sigmasq=sigmasq,
                                      Theta=Theta,Beta1=Beta1,Beta0=Beta0,
                                      Alpha1=Alpha1,Alpha0=Alpha0,
                                      Gamma=Gamma,BetaT=BetaT,AlphaD=AlphaD,
                                      niter=niter,S=S,R=R)
    return(runModel)
    
}



ThetaUpdateNormal = function(Yobs,X,Z,Tr,Trstar,Theta,sigmasq0,sigmasq1,Beta1,Alpha1,
                             Beta0,Alpha0,AlphaD,Gamma,BetaT,
                             n=n,mu=mu,tau=1)
{
    for(ii in 1:n){
        if(Tr[ii] == 1){
            expY = (sum(X[ii,]%*%Beta1)+Alpha1*Theta[ii])
            expTrstar = AlphaD*Theta[ii]+Gamma*Z[ii]+sum(X[ii,]%*%BetaT)
            postVar = 1/(AlphaD^2+(1/tau)+(Alpha1^2/sigmasq1))
            
            postMean = mu + postVar*(AlphaD*(Trstar[ii]-expTrstar)+(Alpha1/sigmasq1)*(Yobs[ii]-expY))
            Theta[ii] = rnorm(1,mean = postMean,sd = sqrt(postVar))
        }
        
        if(Tr[ii] == 0){   
            expY = (sum(X[ii,]%*%Beta0)+Alpha0*Theta[ii])
            expTrstar = AlphaD*Theta[ii]+Gamma*Z[ii]+sum(X[ii,]%*%BetaT)
            postVar = 1/(AlphaD^2+(1/tau)+(Alpha0^2/sigmasq0))
            
            postMean = mu + postVar*(AlphaD*(Trstar[ii]-expTrstar)+(Alpha0/sigmasq0)*(Yobs[ii]-expY))
            Theta[ii] = rnorm(1,mean = postMean,sd = sqrt(postVar))
        }
    }
    return(Theta)
}
##MCMC sampler
MCMClatentTrEff_Normal = function(Yobs,X,Z,Tr,Trstar,n,pp,
                                  sigmasq_prior,sigmasq_prior_alpha,
                                  sigmasq,
                                  Theta,Beta1,Beta0,
                                  Alpha1,Alpha0,
                                  AlphaD,Gamma,BetaT,
                                  niter,S,R)
{
    
    Y1 = Yobs[which(Tr == 1)]
    Y0 = Yobs[which(Tr == 0)]
    Beta1Final = array(NA,dim=c(pp,niter))
    Alpha1Final = rep(NA,niter)
    Beta0Final = array(NA,dim=c(pp,niter))
    Alpha0Final = rep(NA,niter)
    ThetaFinal =  array(NA,dim=c(n,niter))
    AlphaDFinal = rep(NA,niter)
    GammaFinal = rep(NA,niter)
    BetaTFinal = array(NA,dim=c(pp,niter))
    
    sigmasq0 = sigmasq1 = sigmasq
    for(iter in 1:niter){
        XOutcome = as.matrix(data.frame(Theta,X))
        XOutcome1 = XOutcome[which(Tr == 1),]
        XOutcome0 = XOutcome[which(Tr == 0),]
        TrOutcome = as.matrix(data.frame(Theta,Z,X))
        
        posteriorBetasOutcome1 = Gibbs_Reg(Y1,XOutcome1,sigmasq_prior,
                                           sigmasq_prior_alpha,sigmasq1,pp+1)
        
        Alpha1 = posteriorBetasOutcome1[1]
        Beta1 = posteriorBetasOutcome1[-1]
        
        posteriorBetasOutcome0 = Gibbs_Reg(Y0,XOutcome0,sigmasq_prior,
                                           sigmasq_prior_alpha,sigmasq0,pp+1)
        Alpha0 = posteriorBetasOutcome0[1]
        Beta0 = posteriorBetasOutcome0[-1]
        
        posteriorBetasTr = Gibbs_Reg(Trstar,TrOutcome,sigmasq_prior=0.5,
                                     sigmasq_prior_alpha,1,pp+2)
        AlphaD = posteriorBetasTr[1]
        Gamma = posteriorBetasTr[2]
        BetaT = posteriorBetasTr[-c(1,2)]
        
        mean1 = AlphaD*Theta[which(Tr==1)] +
            Gamma*Z[which(Tr==1)] +X[which(Tr==1),]%*% BetaT
        mean0 = AlphaD*Theta[which(Tr==0)] +
            Gamma*Z[which(Tr==0)] + X[which(Tr==0),]%*%BetaT
        
        Trstar[which(Tr==0)] = rtnorm(length(which(Tr==0)),
                                      mean=mean0,
                                      sd=1,  upper=0)
        Trstar[which(Tr==1)] = rtnorm(length(which(Tr==1)),
                                      mean=mean1,
                                      sd=1, lower=0)
        
        
        Rhat1 = R + length(Y1)/2
        Rhat0 = R + length(Y0)/2
        
        Shat1 = S + 1/2*( sum((Y1 - (XOutcome1[,-1]%*%Beta1+ XOutcome1[,1]*Alpha1))^2))
        
        Shat0 = S + 1/2 *(sum((Y0 - (XOutcome0[,-1]%*%Beta0+ XOutcome0[,1]*Alpha0))^2))
        
        sigmasq1 = 1/(rgamma(1,Rhat1,Shat1))
        sigmasq0 = 1/(rgamma(1,Rhat0,Shat0))
        
        ThetaUpdt = ThetaUpdateNormal(Yobs=Yobs,X=X,Z=Z,Tr=Tr,Trstar=Trstar,
                                      Theta=Theta,sigmasq1=sigmasq1,
                                      sigmasq0=sigmasq0,
                                      Beta1=Beta1,Alpha1=Alpha1,
                                      Beta0=Beta0,Alpha0=Alpha0,AlphaD=AlphaD,
                                      Gamma=Gamma,BetaT=BetaT,
                                      n=n,mu=0,tau = 1)
        Theta = ThetaUpdt
        #random sign switch for identification
        switch_sign = rbinom(n=1,1,0.5) - 1
        Alpha1 = switch_sign*Alpha1
        Alpha0 = switch_sign*Alpha0
        AlphaD = switch_sign*AlphaD
        Theta = switch_sign*Theta
        
        Beta1Final[,iter] = Beta1
        Alpha1Final[iter] = Alpha1
        Beta0Final[,iter] = Beta0
        Alpha0Final[iter] = Alpha0
        ThetaFinal[,iter] = Theta
        AlphaDFinal[iter] = AlphaD
        GammaFinal[iter]= Gamma
        BetaTFinal[,iter] = BetaT
    }
    # acc = acc/niter
    draws = list(Beta1=Beta1Final,Alpha1=Alpha1Final,Beta0=Beta0Final,
                 Alpha0=Alpha0Final,Theta=ThetaFinal,
                 AlphaD=AlphaDFinal,Gamma=GammaFinal,BetaTFinal=BetaTFinal)    
    return(draws)
}

