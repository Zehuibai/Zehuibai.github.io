ssc.logrank.Hmisc <- function(tref,   
                              n,     
                              mc,
                              mi,
                              accrual,
                              tmin,   
                              noncomp.c=0,
                              noncomp.i=0,
                              alpha=.05,  
                              nc, ni,
                              GeorgeDesu = TRUE,
                              pr=TRUE)
{
  if(mc>1)
    stop("mc (control group mortality ) should be a fraction")
  
  ## Find mortality in intervention group
  # mi <- (1-r/100)*mc
  if(mi>1)
    stop("mi (intervention group mortality) should be a fraction")
  
  if(missing(nc) | missing(ni)) {
    nc <- n/2; ni <- n/2
  } else n <- nc+ni
  
  if(pr) {
    cat("\nAccrual duration:",accrual,"years  Minimum follow-up:",tmin,"years\n")
    cat("\nTotal sample size:",n,"\n")
    cat("\nAlpha=",alpha,"\n")
    d <- c("Control","Intervention")
    m <- c(mc,mi)
    names(m) <- d
    cat("\n",tref,"-year Mortalities (Events Rate)\n",sep=""); print(m)
  }
  
  ## Find exponential hazards for all groups
  lamc <- -logb(1-mc)/tref
  lami <- -logb(1-mi)/tref
  
  if(pr) {
    lam <- c(lamc,lami)
    names(lam) <- d
    cat("\nHazard Rates\n");
    print(lam)
  }
  
  ## Find probability that a subject will have her event observed during
  ## the study, for all groups
  tmax <- tmin+accrual
  pc <- if(accrual==0)
    1-exp(-lamc*tmin)
  else
    1-1/accrual/lamc*(exp(-tmin*lamc)-exp(-tmax*lamc))
  
  pi <- if(accrual==0)
    1-exp(-lami*tmin)
  else
    1-1/accrual/lami*(exp(-tmin*lami)-exp(-tmax*lami))
  
  if(pr) {
    p <- c(pc,pi)
    names(p) <- d
    cat("\nProbabilities of an Event During Study\n")
    print(p)
  }
  
  ## Find expected number of events, all groups
  mc <- pc*nc
  mi <- pi*ni
  
  if(pr) {
    m <- c(mc,mi)
    names(m) <- d
    cat("\nExpected Number of Events\n")
    print(round(m,1))
  }
  
  ## Find expected value of observed log hazard ratio
  delta <- logb(lami/lamc)
  if(pr)
    cat("\nHazard ratio:",format(exp(delta)),"\n")
  
  if(noncomp.c+noncomp.i>0) {
    if(pr)
      cat("\nDrop-in rate (controls):",noncomp.c,
          "%\nNon-adherence rate (intervention):",noncomp.i,"%\n",sep="")
    
    delta <- delta * (1 - (noncomp.c+noncomp.i)/100)
    if(pr)
      cat("\nEffective hazard ratio with non-compliance:",
          format(exp(delta)),"\n")
  }
  
  ## Find its variance
  ## Schoenfeld approximates the variance of the log hazard ratio by 4/m, where m is the total number of events
  ## George-Desu method uses the slightly better 1/m1 + 1/m2
  if(GeorgeDesu){
    v <- 1/mc + 1/mi
  } else{
    ## Get same as /sasmacro/samsizc.sas if use 4/(mc+mi)
    v <- 4/(mc + mi)
  }
  
  sd <- sqrt(v)
  if(pr)
    cat("Standard deviation of log hazard ratio:",format(sd),"\n")
  if(pr)
    if(GeorgeDesu){
      cat("\nApproximation method of variance of the log hazard ratio based on Peterson B, George SL: Controlled Clinical Trials 14:511–522; 1993.","\n", "\n")
    } else{
      cat("\nApproximation method of variance of the log hazard ratio based on Schoenfeld D: Biometrics 39:499–503; 1983.","\n","\n")
    }
  
  z <- -qnorm(alpha/2)
  
  if(pr){
    cat("Power:",1 - (pnorm(z - abs(delta)/sd) - pnorm(-z - abs(delta)/sd)))
  }else{
    c(Power = 1 - (pnorm(z - abs(delta)/sd) - pnorm(-z - abs(delta)/sd)))
  }
  
}

# ssc.logrank.Hmisc(tref=5,n=950,mc=0.18,mi=0.10,
#                   accrual=1.5,tmin=5,
#                   noncomp.c=10,noncomp.i=15)
# 
# 
# ssc.logrank.Hmisc(tref=5,n=950,mc=0.18,mi=0.10,
#                   accrual=1.5,tmin=5,
#                   noncomp.c=10,noncomp.i=15,
#                   pr=F)
# 
# ssc.logrank.Hmisc(tref=5,nc=450, ni=500,
#                   mc=0.18,mi=0.10,
#                   accrual=1.5,tmin=5,
#                   noncomp.c=10,noncomp.i=15,GeorgeDesu = F)
