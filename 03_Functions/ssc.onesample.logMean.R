## ssc.onesample.logMean

ssc.onesample.logMean <- function(HR, sig.level = 0.05, power = 0.8, pr=TRUE) {
  
  z.alpha = qnorm(sig.level, lower.tail=F)
  z.beta = qnorm(1-power, lower.tail=F)
  num = (z.alpha + z.beta)^2
  denom = (log(HR))^2
  if(pr){
    cat("\nHazard ratio:",format(HR),"\n")
    cat("Alpha (one-sided):",sig.level,"\n")
    cat("Power:",power*100,"%\n")
    cat("\nLog-mean based approach","\n")
    cat("Expected number of events:", ceiling(num/denom))
  }else{
    ceiling(num/denom)
  }
    
}
## Log-mean based approach
## Expected number of events
# ssc.onesample.logMean(HR = 1.5, sig.level = 0.05, power = 0.8)
# ssc.onesample.logMean(HR = 1.5, sig.level = 0.05, power = 0.8, pr = F)

prob.event = function(lambda, accrual, followup){
  term1 = exp(-lambda*followup) - exp(-lambda*(accrual + followup))
  probDeath = 1 - (term1/(accrual*lambda))
  probDeath
}

ssc.onesample.logMean2 <- function(HR, 
                                   sig.level = 0.05, 
                                   power = 0.8, 
                                   lambda, 
                                   accrual, 
                                   followup) {
  
  z.alpha = qnorm(sig.level, lower.tail=F)
  z.beta = qnorm(1-power, lower.tail=F)
  num = (z.alpha + z.beta)^2
  denom = (log(HR))^2
  
  ## Probability of death
  pd = prob.event(lambda=lambda, accrual=accrual, followup=followup)
  cat("Expected number of events:", ceiling(num/denom))
  cat("\n")
  cat("Probability of event:", round(pd,3))
  cat("\n")
  cat("Expected number of patients:", ceiling( ceiling(num/denom)/pd) )
}

# ssc.onesample.logMean2(HR = 1.5, sig.level = 0.05, power = 0.8, lambda=0.10, accrual=2, followup=3)
