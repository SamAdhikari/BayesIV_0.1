#' mcmcRun_Normal_DPM: MCMC sampler for IV Analysis with Dirichlet process mixture prior on the errors and Normal prior on the latent factors
#'
#' @author Sam Adhikari
#' @import msm
#' @import DPpackage
#' @import Rcpp
#' @param Yobs : Input numeric vector of observed outcome of length n, can be binary or continuous. 
#' @param Tr : Binary numeric vector of treatment indicator.
#' @param X : numeric matrix of covariates of dimension n-By-p
#' @param Z : numeric vector of instrumental variable of length n
#' @param niter : number of MCMC sampler to be run
#' @param priors: list to specify the parameters for the prior distribution. If NULL, default specification is used. 
#' @param initialVals: list of initial values for each parameter. If NULL, default setting for initialization is used.
#' @return runModel : list of posterior samples for each parameter
#' @examples { 
#'	obj = mcmcRun_Normal_DPM(Yobs=Yobs,Tr=D,X=X[,-1],Z=Z,niter=10,priors=NULL,initialVals=NULL)
#'	}
#' @export mcmcRun_Normal_DPM



########
########
##Wrapper function to run the MCMC sampler
########
mcmcRun_Normal_DPM = function(Yobs,Tr,X,Z,niter,priors=NULL,
                              initialVals=NULL)
{
    n = length(Yobs)
    pp = dim(X)[2]
    if(is.null(dim(X))){pp = 1}
    
    if(is.null(priors)){
        #set priors
        sigmasq_prior = 100
        sigmasq_prior_alpha = 100
        muTheta = 0
        tauTheta = 1
        #Default DPM priors
        prior0 = list(a0=1,b0=1,m1=rep(0,1),psiinv1=diag(5,1),nu1=1,
                      tau1=1,tau2=10)
        prior1 = list(a0=1,b0=1,m1=rep(0,1),psiinv1=diag(5,1),nu1=1,
                      tau1=1,tau2=10)
    }else{
        sigmasq_prior = priors$sigmasq_prior
        sigmasq_prior_alpha = priors$sigmasq_prior_alpha
        muTheta = priors$muTheta
        tauTheta = priors$tauTheta
        prior0 = priors$prior0
        prior1 = priors$prior1
    }
    
    if(is.null(initialVals)){
        ##Initialization
        Beta1 = rnorm(pp,0,1)
        Alpha1 = rnorm(1,0,1)
        Beta0 = rnorm(pp,0,1)
        Alpha0 = rnorm(1,0,1)
        Theta = rnorm(n,0,1)
        AlphaD = rnorm(1,0,1)
        BetaT = rnorm(pp+1,0,1)
        Gamma = rnorm(1,0,1)
        state0 = state1 = NULL
        status0 = status1 = TRUE
        
    }else{
        Beta1 = initialVals$Beta1
        Alpha1 = initialVals$Alpha1
        Beta0 = initialVals$Beta0
        Alpha0 = initialVals$Alpha0
        Theta = initialVals$Theta
        AlphaD = initialVals$AlphaD
        BetaT = initialVals$BetaT
        Gamma = initialVals$Gamma
        state1 = initialVals$state1
        state0 = initialVals$state0
        status1 = FALSE
        status0 = FALSE
    }
    
    XTr = as.matrix(data.frame(Int = rep(1,length(Tr)),X))
    Trstar = Tr
    if(pp==0){
        mean1 = AlphaD*Theta[which(Tr==1)] +
            Gamma*Z[which(Tr==1)] +XTr[which(Tr==1),]* BetaT
        
        mean0 = AlphaD*Theta[which(Tr==0)] +
            Gamma*Z[which(Tr==0)] + XTr[which(Tr==0),]* BetaT
    }
    if(pp > 0){
        mean1 = AlphaD*Theta[which(Tr==1)] +
            Gamma*Z[which(Tr==1)] +XTr[which(Tr==1),]%*% BetaT
        mean0 = AlphaD*Theta[which(Tr==0)] +
            Gamma*Z[which(Tr==0)] + XTr[which(Tr==0),]%*%BetaT
    }
    Trstar[which(Tr==0)] = rtnorm(length(which(Tr==0)),
                                  mean=mean0,
                                  sd=1, upper=0)
    Trstar[which(Tr==1)] = rtnorm(length(which(Tr==1)),
                                  mean=mean1,
                                  sd=1, lower=0)
    #Run MCMC
    runModel = MCMClatentTrEff_Normal_DPM(
        Yobs=Yobs,X=X,Z=Z,Tr=Tr,Trstar=Trstar,
        n=n,pp=pp,
        sigmasq_prior=sigmasq_prior,
        sigmasq_prior_alpha = sigmasq_prior_alpha,
        prior0=prior0,prior1=prior1,
        muTheta=muTheta,tauTheta=tauTheta,
        Theta=Theta,Beta1=Beta1,Beta0=Beta0,
        Alpha1=Alpha1,Alpha0=Alpha0,
        Gamma=Gamma,BetaT=BetaT,AlphaD=AlphaD,
        niter=niter,
        state1=state1,state0=state0,
        status1=status1,status0=status0)
    return(runModel)
}

UpdateCoefficients_Normal_DPM = function(YY,XX,Betas,mu0,sigmasq0,muY,sigmasqY){
    for(kk in 1:length(Betas)){
        betanotK = Betas[-kk]
        XXnotK = XX[,-kk]
        XK = XX[,kk]
        muPart = muY + XXnotK%*%betanotK
        posterior_var = 1/(sum((XK*XK)/sigmasqY) + 1/sigmasq0) 
        posterior_mean = (sum((YY-muPart)*XK /sigmasqY) + mu0/sigmasq0) *
            posterior_var
        Betas[kk] = rnorm(1,posterior_mean,sd=sqrt(posterior_var))
    }
    #  print(Betas)
    return(Betas)
}

##MCMC sampler
MCMClatentTrEff_Normal_DPM = function(Yobs,X,Z,Tr,Trstar,n,pp,
                                      sigmasq_prior,sigmasq_prior_alpha,
                                      prior0,prior1,
                                      muTheta,tauTheta,
                                      Theta,Beta1,Beta0,
                                      Alpha1,Alpha0,
                                      AlphaD,Gamma,BetaT,
                                      niter,
                                      state1,state0,
                                      status1,status0)
{
    Y1 = Yobs[which(Tr == 1)]
    Y0 = Yobs[which(Tr == 0)]
    XOutcome = as.matrix(data.frame(Theta,X))
    XOutcome1 = XOutcome[which(Tr == 1),]
    XOutcome0 = XOutcome[which(Tr == 0),]
    XTr = as.matrix(data.frame(rep(1,length(Tr)),X))
    TrOutcome = as.matrix(data.frame(Theta,Z,XTr))
    Beta1Final = array(NA,dim=c(pp,niter))
    Alpha1Final = rep(NA,niter)
    Beta0Final = array(NA,dim=c(pp,niter))
    Alpha0Final = rep(NA,niter)
    ThetaFinal =  array(NA,dim=c(n,niter))
    AlphaDFinal = rep(NA,niter)
    GammaFinal = rep(NA,niter)
    BetaTFinal = array(NA,dim=c(pp+1,niter))
    Y1hatFinal = rep(NA,niter)
    Y0hatFinal = rep(NA,niter)
    sigmasqY0Final = muY0Final = 
        array(NA,dim=c(length(which(Tr==0)),niter))
    sigmasqY1Final = muY1Final = 
        array(NA,dim=c(length(which(Tr==1)),niter))
    loglikeAll = rep(NA,niter)
    # MCMC parameters
    nburn = 0
    nsave = 1
    nskip = 0
    ndisplay = 0
    mcmc = list(nburn=nburn,nsave=nsave,nskip=nskip,
                ndisplay=ndisplay)
    for(iter in 1:niter){
        ##update coefficients for treated group
        param1 = c(Alpha1,Beta1)
        #NP residual regression for treated outcome
        Youtcome1 = Y1-(XOutcome1%*%param1)
        ##Code to specify grid on density estimations
        ##Default options run into error when sample size is large
        # Youtcome1[which(is.na(Youtcome1))] = mean(Youtcome1,na.rm = TRUE)
        miny <- min(Youtcome1,na.rm = TRUE)
        maxy <- max(Youtcome1,na.rm = TRUE)
        vary <- var(Youtcome1,na.rm=TRUE)
        if(is.na(vary)|is.nan(vary)|!is.finite(vary)) {vary = 10E300} #place holder for variance
        if(is.na(miny)|is.nan(miny)|!is.finite(miny)) {miny = -10E300}
        if(is.na(maxy)|is.nan(maxy)|!is.finite(maxy)) {maxy = 10E300}
        from =miny-0.25*sqrt(vary)-50000
        to = maxy+0.25*sqrt(vary)+50000
          
        grid <- seq(from=from[1],to=to[1],length.out=5000)
        
        fit1 = DPdensity(y=Youtcome1,grid=grid,
                         prior=prior1,mcmc=mcmc,
                         state=state1,status=status1,na.action = na.omit)
        
        state1 = fit1$state
        status1=FALSE #continuation of previous chain
        ccount1 = seq(0,(length(fit1$save.state$randsave)-4)/2,by=1)
        muY1 = fit1$save.state$randsave[1,2*ccount1 + 1]
        sigmasqY1 = fit1$save.state$randsave[1,2*ccount1 + 2]
        
        UpdateBetasOutcome1 = UpdateCoefficients_Normal_DPM(YY=Y1,
                                                            XX=XOutcome1,
                                                            Betas =param1,
                                                            mu0=0,sigmasq0= sigmasq_prior,
                                                            muY=muY1,sigmasqY=sigmasqY1)
        posteriorBetasOutcome1 = UpdateBetasOutcome1
        Alpha1 = posteriorBetasOutcome1[1]
        Beta1 = posteriorBetasOutcome1[-1]
        param1 = c(Alpha1,Beta1)
        loglikelihoodOld1 = sum(sapply(1:length(Youtcome1),function(x){dnorm(Youtcome1[x],
                                                 (XOutcome1[x,]*param1+
                                                      muY1[x]),sqrt(sigmasqY1[x]),log=TRUE)}))
        ##update coefficients for control group
        ##NP residual regression 
        param0 = c(Alpha0,Beta0)
        Youtcome0 = Y0-(XOutcome0%*%param0)
        miny0 <- min(Youtcome0,na.rm = TRUE)
        maxy0 <- max(Youtcome0,na.rm = TRUE)
        vary0 <- var(Youtcome0,na.rm=TRUE)
        ###My trick to make sure we don't run into errors when sample size is big
        if(is.na(vary0)|is.nan(vary0)|!is.finite(vary0)) {vary0 = 10E300} #place holder for variance
        if(is.na(miny0)|is.nan(miny0)|!is.finite(miny0)) {miny0 = -10E300}
        if(is.na(maxy0)|is.nan(maxy0)|!is.finite(maxy0)) {maxy0 = 10E300}
        
        from0 =miny0 - 0.25*sqrt(vary0)-50000
        to0 = maxy0 + 0.25*sqrt(vary0)+50000
        
        grid0 <- seq(from=from0[1],to=to0[1],length.out=5000)
        
        fit0 = DPdensity(y=Youtcome0,grid=grid0,
                         prior=prior0,mcmc=mcmc,
                         state=state0,status=status0,na.action = na.omit)
        state0 = fit0$state
        status0=FALSE #continuation of previous chain
        ccount0 = seq(0,(length(fit0$save.state$randsave)-4)/2,by=1)
        muY0 = fit0$save.state$randsave[1,2*ccount0 + 1]
        sigmasqY0 = fit0$save.state$randsave[1,2*ccount0 + 2]
        muY0 = DPrandom(fit0)[[1]][,1] #group specific mean
        sigmasqY0 = DPrandom(fit0)[[1]][,2] #group specific variance
        UpdateBetasOutcome0 = UpdateCoefficients_Normal_DPM(
            YY=Y0,XX=XOutcome0,
            Betas=param0,
            mu0=0,sigmasq0=sigmasq_prior,
            muY=muY0,sigmasqY=sigmasqY0)
        
        posteriorBetasOutcome0 = UpdateBetasOutcome0
        Alpha0 = posteriorBetasOutcome0[1]
        Beta0 = posteriorBetasOutcome0[-1]
        param0 = c(Alpha0,Beta0)
        loglikelihoodOld0 = sum(sapply(1:length(Youtcome0),function(x){dnorm(Youtcome0[x],
                                                 (XOutcome0[x,]*param0+
                                                      muY0[x]),sqrt(sigmasqY0[x]),log=TRUE)}))
        # ####
        ####Updating parameters for treatment equation
        ####by probit regression
        posteriorBetasTr = Gibbs_Reg(Trstar,TrOutcome,
                                     sigmasq_prior=0.5,
                                     sigmasq_prior_alpha,1,pp+3)
        AlphaD = posteriorBetasTr[1]
        Gamma = posteriorBetasTr[2]
        BetaT = posteriorBetasTr[-c(1,2)]
        if(pp==0){
            mean1 = AlphaD*Theta[which(Tr==1)] +
                Gamma*Z[which(Tr==1)] +XTr[which(Tr==1),]* BetaT 
            mean0 = AlphaD*Theta[which(Tr==0)] +
                Gamma*Z[which(Tr==0)] + XTr[which(Tr==0),]* BetaT 
        }
        if(pp > 0){
            mean1 = AlphaD*Theta[which(Tr==1)] +
                Gamma*Z[which(Tr==1)] +XTr[which(Tr==1),]%*% BetaT
            mean0 = AlphaD*Theta[which(Tr==0)] +
                Gamma*Z[which(Tr==0)] + XTr[which(Tr==0),]%*%BetaT
        }
        Trstar[which(Tr==0)] = rtnorm(length(which(Tr==0)),
                                      mean=mean0,
                                      sd=1,  upper=0)
        Trstar[which(Tr==1)] = rtnorm(length(which(Tr==1)),
                                      mean=mean1,
                                      sd=1, lower=0)
        #Normal CDF
        
        loglikelihoodOldTR = sum(log(pnorm(mean1))) + sum(log(1-pnorm(mean0)))
        
        loglikeAll[iter] = loglikelihoodOld1 + loglikelihoodOld0 + loglikelihoodOldTR        
        
        ##Updating the latent factor Theta
        XTr = as.matrix(data.frame(Int = rep(1,length(Tr)),X))
        muY = rep(0,n)
        muY[which(Tr==1)] = muY1
        muY[which(Tr==0)] = muY0
        sigmaSq = rep(0,n)
        sigmaSq[which(Tr==1)] = sigmasqY1
        sigmaSq[which(Tr==0)] = sigmasqY0
        
        #  ThetaUpdt = ThetaUpdateNormal_DPM(Yobs=Yobs,X=X,Z=Z,Tr=Tr,
        #                                Trstar=Trstar,
        #                               Theta=Theta,sigmasqY1=sigmasqY1,
        #                              sigmasqY0=sigmasqY0,
        #                               muY1=muY1,muY0=muY0,
        #                               Beta1=Beta1,Alpha1=Alpha1,
        #                               Beta0=Beta0,Alpha0=Alpha0,
        #                               AlphaD=AlphaD,
        #                               Gamma=Gamma,BetaT=BetaT,
        #                               n=n,mu=muTheta,tau=tauTheta)
        ThetaUpdt =  ThetaUpdateNormal_DPM(Yobs=Yobs,XX=X,XTr=XTr, 
                                           Z=Z, Tr=Tr, Trstar=Trstar, Theta=Theta, 
                                           sigmaSq=sigmaSq, muY=muY, Beta1=Beta1, 
                                           Alpha1=Alpha1, Beta0=Beta0, Alpha0=Alpha0, 
                                           AlphaD=AlphaD, Gamma=Gamma, BetaT=BetaT, 
                                           n=n, pp=pp, mu=muTheta, tau=tauTheta) 
        
        Theta = ThetaUpdt
        #print(Theta)
        #random sign switch for identification
        switch_sign = rbinom(n=1,1,0.5) - 1
        Alpha1 = switch_sign*Alpha1
        Alpha0 = switch_sign*Alpha0
        AlphaD = switch_sign*AlphaD
        Theta = switch_sign*Theta
        
        ##Updating the covariates matrix with new Theta
        XOutcome = as.matrix(data.frame(Theta,X))
        XOutcome1 = XOutcome[which(Tr == 1),]
        XOutcome0 = XOutcome[which(Tr == 0),]
        TrOutcome = as.matrix(data.frame(Theta,Z,XTr))
        
        
        
        Beta1Final[,iter] = Beta1
        Alpha1Final[iter] = Alpha1
        Beta0Final[,iter] = Beta0
        Alpha0Final[iter] = Alpha0
        ThetaFinal[,iter] = Theta
        AlphaDFinal[iter] = AlphaD
        GammaFinal[iter]= Gamma
        BetaTFinal[,iter] = BetaT
        # Y1hatFinal[iter] = Y1hat
        # Y0hatFinal[iter] = Y0hat
        muY1Final[,iter] = muY1
        muY0Final[,iter] = muY0
        sigmasqY0Final[,iter] = sigmasqY0
        sigmasqY1Final[,iter] = sigmasqY1
    }
    # acc = acc/niter
    draws = list(Beta1=Beta1Final,Alpha1=Alpha1Final,Beta0=Beta0Final,
                 Alpha0=Alpha0Final,Theta=ThetaFinal,
                 AlphaD=AlphaDFinal,Gamma=GammaFinal,
                 BetaT=BetaTFinal,
                 muY1=muY1Final,muY0=muY0Final,
                 sigmasqY0 = sigmasqY0Final,sigmasqY1=sigmasqY1Final,
                 state1=state1,state0=state0,
                 status1=status1,status0=status0,
                 loglikeAll=loglikeAll)    
    return(draws)
}


