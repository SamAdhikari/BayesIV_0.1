#' plotATEs: Function to plot posterior chain of ATEs
#' 
#' @param xx
#' @param pdfname
#' @param int.quantiles
#' @param effect_true
#' @param effect_obs
#' @param labels.axis
#' @export

plotATEs = function(xx,pdfname=NULL,int.quantiles,effect_true,effect_obs,
                    labels.axis=NULL)
{
  ##also add observed mean difference in these plots
  if(!is.null(pdfname)){
    pdf(file=paste(pdfname,'.pdf',sep=''),
        height=5,width=10)}
  plot(xx, int.quantiles[,3], pch=".",
       ylim=c(min(int.quantiles,effect_true ), 
              max(int.quantiles,effect_true )),
       col="white", xlab="",
       ylab="",
       main="Treatment Effects",xaxt='n',cex.main=2,cex.axis=1.5)
  axis(side=1,at = xx-0.08,labels=labels.axis,cex.axis=2)
  rect((xx-0.05), int.quantiles[,1], (xx+0.05),
       int.quantiles[,3], col='black', border=NA)
  segments((xx-0.15), int.quantiles[,2], (xx+0.15),
           int.quantiles[,2], lwd=2, col='blue')
  segments((xx-0.15), effect_true, (xx+0.15),
           effect_true, lwd=2, col='red')
  segments((xx-0.15),effect_obs, (xx+0.15),
           effect_obs, lwd=2, col='green')
  legend('topleft',legend=c('True','Observed','Estimated'),
         col=c('red','green','blue'),lty=1,lwd=3,cex=2,bty='n') 
  if(!is.null(pdfname))dev.off()
}