np.gmeta <-
function(Thetahat, se, alpha=c(0.025,0.975), n, m, band_pwr = 0.5, resample=200, B=40, len=10)
{
 out<-.Fortran("confdistint",K=as.integer(length(Thetahat)),
               Thetahat = as.double(Thetahat),se = as.double(se),
               n=as.double(n),resample=as.integer(resample),B=as.integer(B),len=as.integer(len),
               m=as.integer(m),band_pwr=as.double(band_pwr),
               shrink=as.double(0),smoothlist=as.double(rep(0,len)),
               eligids=as.integer(rep(0,len^2)),distvs=as.double(rep(0,len^2)),
               finalxi1=as.double(rep(0,resample)),finalxi2=as.double(rep(0,resample)),finalxi3=as.double(rep(0,resample))
               )
 percentiles = matrix(0,3,length(alpha))
 colnames(percentiles) = as.character(alpha); rownames(percentiles) = c("min.unif","max.smooth","mean.smooth")
 percentiles[1,] = quantile(out$finalxi1,alpha)
 percentiles[2,] = quantile(out$finalxi2,alpha)
 percentiles[3,] = quantile(out$finalxi3,alpha)
 aa=out$eligids;aa=aa[aa!=0]
 return( list(percentiles = percentiles, shrink = out$shrink, smoothlist= out$smoothlist,
             distance = matrix(out$distvs,10,10), elig.ind = arrayInd(aa,c(len,len)))             
       )
}

