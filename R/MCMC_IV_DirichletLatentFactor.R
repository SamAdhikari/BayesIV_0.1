#'mcmcRunDP: MCMC sampler for IV Analysis with Dirichlet process mixture prior on the latent factors
#'
#' @author Sam Adhikari
#' @import msm
#' @import DPpackage
#' @import Rcpp
#' @param Yobs: Input numeric vector of observed outcome of length n, can be binary or continuous. 
#' @param Tr : Binary numeric vector of treatment indicator.
#' @param X : numeric matrix of covariates of dimension n-By-p
#' @param Z : numeric vector of instrumental variable of length n
#' @param niter : number of MCMC sampler to be run
#' @return runModel : list of posterior samples for each parameter, acceptance rate and tuning parameter
#' @examples \dontrun{obj = mcmcRunDP(Yobs,Tr,X,Z,niter)}
#' @export mcmcRunDP


## library(msm) ##for sampling from truncated normal distribution
## library(DPpackage) ##we will use this package to sample from the prior distribution of the
##                     ##latent factors


##MCMC sampler for IV Analysis with Dirichlet process mixture prior on the 
##      latent factors
#  prior specification	
#  a list giving the prior information.
# a0 and b0: hyperparameters for prior distribution of
#   the precision parameter of the Dirichlet process prior, 
#alpha: the value of the precision parameter 
#    (it must be specified if a0 is missing, see details below),
#nu2 and psiinv2: the hyperparameters of the 
#     inverted Wishart prior distribution for the scale matrix,
#Psi1: the inverted Wishart part of the baseline distribution,
#tau1 and tau2: the hyperparameters for the gamma prior distribution 
#      of the scale parameter k0 of the normal part of the baseline distribution, 
#m2 and s2: the mean and the covariance of the normal prior 
#      for the mean, m1 of the normal component of the baseline distribution,
#nu1 and psiinv1 (it must be specified if nu2 is missing)
#      giving the hyperparameters of the inverted Wishart part 
#      of the baseline distribution and,
#m1 giving the mean of the normal part of the baseline 
#       distribution (it must be specified if m2 is missing, 
#       see details below) and,
#k0 giving the scale parameter of the normal part of the
#         baseline distribution (it must be specified if tau1
#          is missing, see details below).

# we constrain the variance of the normal distribution
# by using constraining parameters on the inverse wishart prior
# distribution by specifyin 'psiinv' and 'nu1'
#prior1 = list(alpha=1,m1=rep(0,1),psiinv1=diag(1,1),nu1=1,
#     tau1=1,tau2=100)

##Function to run the MCMC sampler
mcmcRunDP = function(Yobs,Tr,X,Z,niter)
{
    n = length(Yobs)
    pp = dim(X)[2]
    if(is.null(pp))pp = 1
    #set priors
    sigmasq_prior = 1000
    sigmasq_prior_alpha = 100
    sigmasq = 1
    S = 15
    R = 10
    ##Initialization
    Beta1 = rnorm(pp,0,1)
    Alpha1 = rnorm(1,0,1)
    Beta0 = rnorm(pp,0,1)
    Alpha0 = rnorm(1,0,1)
    Theta = rnorm(n,0,0.1)
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
    runModel = MCMClatentTrEff(Yobs=Yobs,X=X,Z=Z,Tr=Tr,Trstar=Trstar,
                               n=n,pp=pp,
                               sigmasq_prior=sigmasq_prior,
                               sigmasq_prior_alpha=sigmasq_prior_alpha,
                               sigmasq=sigmasq,
                               Theta=Theta,Beta1=Beta1,Beta0=Beta0,
                               Alpha1=Alpha1,Alpha0=Alpha0,
                               Gamma=Gamma,BetaT=BetaT,AlphaD=AlphaD,
                               niter=niter,S=S,R=R)
    return(runModel)
}



ThetaUpdate = function(Yobs,X,Z,Tr,Trstar,Theta,sigmasq0,sigmasq1,Beta1,Alpha1,
                       Beta0,Alpha0,AlphaD,Gamma,BetaT,
                       n=n,mu=mu,tau=tau)
{
    for(ii in 1:n){
        if(Tr[ii] == 1){
            expY = (sum(X[ii,]%*%Beta1)+Alpha1*Theta[ii])
            expTrstar = AlphaD*Theta[ii]+Gamma*Z[ii]+sum(X[ii,]%*%BetaT)
            postVar = 1/(AlphaD^2+(1/tau[ii])+(Alpha1^2/sigmasq1))
            postMean = mu[ii] + postVar*(AlphaD*(Trstar[ii]-expTrstar)+
                                             (Alpha1/sigmasq1)*(Yobs[ii]-expY))
            Theta[ii] = rnorm(1,mean = postMean,sd = sqrt(postVar))
        }
        if(Tr[ii] == 0){
            expY = (sum(X[ii,]%*%Beta0)+Alpha0*Theta[ii])
            expTrstar = AlphaD*Theta[ii]+Gamma*Z[ii]+sum(X[ii,]%*%BetaT)
            postVar = 1/(AlphaD^2+(1/tau[ii])+(Alpha0^2/sigmasq0))
            postMean = mu[ii] + postVar*(AlphaD*(Trstar[ii]-expTrstar)+
                                             (Alpha0/sigmasq0)*(Yobs[ii]-expY))
            Theta[ii] = rnorm(1,mean = postMean,sd = sqrt(postVar))
        }
    }
    return(Theta)
}

##################
##MCMC sampler
##
MCMClatentTrEff = function(Yobs,X,Z,Tr,Trstar,n,pp,
                           sigmasq_prior,sigmasq_prior_alpha,
                           sigmasq,
                           Theta,Beta1,Beta0,
                           Alpha1,Alpha0,
                           AlphaD,Gamma,BetaT,
                           niter,S,R){
    Y1 = Yobs[which(Tr == 1)]
    Y0 = Yobs[which(Tr == 0)]
    ##place holders to save the posterior draws
    Beta1Final = array(NA,dim=c(pp,niter))
    Alpha1Final = rep(NA,niter)
    Beta0Final = array(NA,dim=c(pp,niter))
    Alpha0Final = rep(NA,niter)
    ThetaFinal =  array(NA,dim=c(n,niter))
    AlphaDFinal = rep(NA,niter)
    GammaFinal = rep(NA,niter)
    BetaTFinal = array(NA,dim=c(pp,niter))
    
    sigmasq0 = sigmasq1 = sigmasq
    ##initialization of the DPM
    state = NULL
    status = TRUE
    for(iter in 1:niter){
        XOutcome = as.matrix(data.frame(Theta,X))
        TrOutcome = as.matrix(data.frame(Theta,Z,X))
        # TrOutcome = as.matrix(data.frame(Z,X))
        ##select the covariates and latent factor for each treatment group
        ##it needs to be updated after updating theta
        XOutcome1 = XOutcome[which(Tr == 1),]
        XOutcome0 = XOutcome[which(Tr == 0),]
        #updating coefficients for the treated
        posteriorBetasOutcome1 = Gibbs_Reg(Y=Y1,X=XOutcome1,
                                           sigmasq_prior=sigmasq_prior,
                                           sigmasq_prior_alpha=sigmasq_prior_alpha,
                                           sigmasq=sigmasq1,#n=length(which(Tr==1)),
                                           PP= (pp+1))
        Alpha1 = posteriorBetasOutcome1[1]
        Beta1 = posteriorBetasOutcome1[-1]
        #updating coefficients for the control
        posteriorBetasOutcome0 = Gibbs_Reg(Y=Y0,X=XOutcome0,
                                           sigmasq_prior=sigmasq_prior,
                                           sigmasq_prior_alpha=sigmasq_prior_alpha,
                                           sigmasq=sigmasq0,#n=length(which(Tr==0)),
                                           PP= (pp+1))
        Alpha0 = posteriorBetasOutcome0[1]
        Beta0 = posteriorBetasOutcome0[-1]
        #updating coefficients for the selection process
        posteriorBetasTr = Gibbs_Reg(Trstar,TrOutcome,
                                     0.5,
                                     sigmasq_prior_alpha,
                                     1,#n,
                                     pp+2)
        AlphaD = posteriorBetasTr[1]
        
        #   AlphaD = AlphaD
        Gamma = posteriorBetasTr[2]
        BetaT = posteriorBetasTr[-c(1,2)]
        ## updating the intermediate variable in the probit
        ##  model for selection
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
        ##updating the variance of the errors in outcomes 
        Rhat1 = R + length(Y1)/2
        Rhat0 = R + length(Y0)/2
        Shat1 = S + 1/2*( sum((Y1 -
                                   (XOutcome1[,-1]%*%Beta1+ XOutcome1[,1]*Alpha1))^2))
        Shat0 = S + 1/2 *(sum((Y0 - 
                                   (XOutcome0[,-1]%*%Beta0+ XOutcome0[,1]*Alpha0))^2))
        sigmasq1 = 1/(rgamma(1,Rhat1,Shat1))
        sigmasq0 = 1/(rgamma(1,Rhat0,Shat0))
        
        #sampling the parameters of the dirichlet mixture prior and membership vector
        #       ##Using package 'DPpackage'
        # Initial state
        # MCMC parameters
        nburn <- 0
        nsave <- 1
        nskip <- 0
        ndisplay <- 0
        mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,
                     ndisplay=ndisplay)
        
        prior1 = list(a0=100,b0=10,m1=rep(0,1),psiinv1=diag(5,1),nu1=1,
                      tau1=1,tau2=100)
                      
        fit1.1 = DPdensity(y=Theta,prior=prior1,mcmc=mcmc,
                           state=state,status=status)
        state = fit1.1$state
        status=FALSE #continuation of previous chain
        mu = DPrandom(fit1.1)[[1]][,1] #group specific mean
        tau = DPrandom(fit1.1)[[1]][,2] #group specific variance
        
        ThetaUpdt = ThetaUpdate(Yobs=Yobs,X=X,Z=Z,Tr=Tr,Trstar=Trstar,
                                Theta=Theta,sigmasq1=sigmasq1,
                                sigmasq0=sigmasq0,
                                Beta1=Beta1,Alpha1=Alpha1,
                                Beta0=Beta0,Alpha0=Alpha0,AlphaD=AlphaD,
                                Gamma=Gamma,BetaT=BetaT,
                                n=n,mu=mu,tau=tau)
        Theta = ThetaUpdt
        
        #random sign switch for identification
        switch_sign = rbinom(n=1,1,0.5) 
        if(switch_sign==0)switch_sign = -1
        
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
    #    acc = acc/niter
    draws = list(Beta1=Beta1Final,Alpha1=Alpha1Final,Beta0=Beta0Final,
                 Alpha0=Alpha0Final,Theta=ThetaFinal,
                 AlphaD=AlphaDFinal,Gamma=GammaFinal,
                 BetaT=BetaTFinal)    
    return(draws)
    rm(list=ls())
    gc(verbose = FALSE)
}
