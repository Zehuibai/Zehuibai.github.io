expLikeRatio = function(d, sig.level, power){
  num = qchisq(sig.level, df=(2*d), lower.tail=F)
  denom = qchisq(power, df=(2*d), lower.tail=F)
  HR = num/denom
  HR
}

ssc.onesample.LR <- function(HR, sig.level = 0.05, power = 0.8) {
  LR = function(d, sig.level, power, HR){
    expLikeRatio(d, sig.level, power) - HR
  }
  # Find the root for the function LR(d)
  result = uniroot(f = LR, lower = 1, upper = 1000,
                   sig.level = sig.level, power = power, HR = HR)
  cat("Expected number of events:", ceiling(result$root))
}


# ssc.onesample.LR(HR = 1.5, sig.level = 0.05, power = 0.8)