#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector ThetaUpdateNormal_DPM(NumericVector Yobs, 
                   NumericVector XX, 
                   NumericVector XTr,
                   NumericVector Z, 
                   NumericVector Tr,
                   NumericVector Trstar,
                   NumericVector Theta, 
                   NumericVector sigmaSq,
                   NumericVector muY,
                   NumericVector Beta1, 
                   double Alpha1, 
                   NumericVector Beta0, 
                   double Alpha0, 
                   double AlphaD, 
                   double Gamma, 
                   NumericVector BetaT,
                   int n, int pp, 
                   double mu, 
                   NumericVector tau)
{
    int ii, kk, ll;
    double sumXB, sumXBT, expY, expTrstar, postVar, postMean;
    NumericVector ThetaNew(n);
    double a = 0;
    // NumericVector muY(n);
    // muY[which(Tr == 1)] = muY1;
    // muY[which(Tr==0)] = muY0;
    // NumericVector si
    for(ii = 0 ; ii <= (n-1) ; ii++){
        if(Tr[ii] == 1){
            sumXB = 0.0;
            for(kk = 0; kk <= (pp-1); kk++){
                sumXB += XX[ii+(n-1)*kk]*Beta1[kk];
            }
            sumXBT = 0.0;
            for(ll = 0; ll<=pp; ll++){
                sumXBT += XTr[ii+(n-1)*kk] * BetaT[kk];
            }
            expY = sumXB + Alpha1*Theta[ii] + muY[ii];
            expTrstar = sumXBT + Gamma*Z[ii] + AlphaD*Theta[ii] ;
            postVar = 1/((AlphaD*AlphaD)+ 1/tau[ii] + 
                (Alpha1*Alpha1)/sigmaSq[ii]);
            postMean = mu + postVar*(AlphaD*(Trstar[ii]-expTrstar)+ 
                (Alpha1/sigmaSq[ii])*(Yobs[ii]-expY));
            a =  as<double>(rnorm(1,postMean,sqrt(postVar)));
            ThetaNew[ii] = a;
        }
        if(Tr(ii) == 0){
            sumXB = 0.0;
            for( kk = 0; kk <= (pp-1); kk++){
                sumXB += XX[ii+(n-1)*kk]*Beta0[kk];
            }
            sumXBT = 0.0;
            for( ll = 0; ll<=pp; ll++){
                sumXBT += XTr[ii+(n-1)*kk] * BetaT[kk];
            }
            expY = sumXB + Alpha0*Theta[ii] + muY[ii];
            expTrstar = sumXBT + Gamma*Z[ii] + AlphaD*Theta[ii] ;
            postVar = 1/((AlphaD*AlphaD)+ 1/tau[ii] + 
                (Alpha0*Alpha0)/sigmaSq[ii]);
            postMean = mu+ postVar*(AlphaD*(Trstar[ii]-expTrstar)+ 
                    (Alpha0/sigmaSq[ii])*(Yobs[ii]-expY));
            a =  as<double>(rnorm(1,postMean,sqrt(postVar)));
            ThetaNew[ii] = a;
        }           
    }
    return Theta;
}




// ##Gibbs sampler function to sample slope coefficients in
// ## linear regression with normal noise
// Gibbs_Reg = function(Y,X,sigmasq_prior,sigmasq_prior_alpha,sigmasq,PP)
// {
//     sigma_prior = diag(c(sigmasq_prior_alpha,rep(sigmasq_prior,(PP-1))),PP)
//     precision_prior = solve(sigma_prior,tol=1E-100)
//     precision_n_inv = solve(t(X)%*%X + precision_prior,tol=1E-100)
//     postMean = (precision_n_inv)%*%t(X)%*%Y
//     postVar = diag(sigmasq,PP) * precision_n_inv
//     return(rnorm(PP,postMean,sqrt(diag(postVar))))
// }
