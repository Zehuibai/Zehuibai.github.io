## based on pwrss.t.2means

# use welch.df = TRUE for unequal variances and unequal n in independent samples t test
# allows r between pre and post to be different from 0, pwr.t.test(...,paired = TRUE) assumes paired.r = 0
# allows diferent sd for pre and post
# default specification assumes standardized means (or mean difference)
# margin is defined as the minimum mu1 - mu2 that is practically relevant
# suggest margin = t_critical * SE by default?
ssc.twosamplemean.ttest <- function (mu1, mu2 = 0, margin = 0,
                            sd1 = ifelse(paired, sqrt(1/(2*(1-paired.r))), 1), sd2 = sd1,
                            kappa = 1, paired = FALSE, paired.r = 0.50,
                            alpha = 0.05, welch.df = FALSE,
                            alternative = c("not equal", "greater", "less",
                                            "equivalent", "non-inferior", "superior"),
                            n2 = NULL, power = NULL, verbose = TRUE)
{
  
  welch_df <- function(sd1, sd2, n1, n2) {
    (sd1^2 / n1 + sd2^2 / n2)^2 /
      (sd1^4 / (n1^2 * (n1 - 1)) + sd2^4 / (n2^2 * (n2 - 1)))
  }
  
  if (length(alternative) > 1)
    alternative <- alternative[1]
  if (is.null(n2) & is.null(power))
    stop("`n2` and `power` cannot be `NULL` at the same time",
         call. = FALSE)
  if (!is.null(n2) & !is.null(power))
    stop("one of the `n2` or `power` should be `NULL`",
         call. = FALSE)
  if (paired & isTRUE(welch.df))
    warning("Welch test does not apply to paired samples",
            call. = FALSE)
  if (!is.null(n2) & is.null(power))
    requested <- "power"
  if (is.null(n2) & !is.null(power))
    requested <- "n2"
  
  if (alternative == "not equal") {
    
    if (margin != 0)
      warning("`margin` argument is ignored")
    
    # sample size
    if (is.null(n2)) {
      HA_H0 <- mu1 - mu2
      beta <- 1 - power
      i <- 0
      tol <- 0.01
      conv <- FALSE
      n2.0 <- 30
      while(i <= 100 & conv == FALSE){
        n1.0 <- n2.0 * kappa
        
        if(paired) {
          df <- n2.0 - 1
        } else {
          ifelse(welch.df,
                 df <- welch_df(sd1, sd2, n1.0, n2.0),
                 df <- n1.0 + n2.0 - 2)
        }
        
        if(df<= 0 | is.infinite(df)){break}
        M <- qt(alpha / 2, df = df, ncp = 0, lower.tail = FALSE) +
          qt(beta, df = df, ncp = 0, lower.tail = FALSE)
        
        if(paired) {
          n2 = M^2 * (sd1^2 + sd2^2 - 2 * sd1 * sd2 * paired.r) / (HA_H0)^2
        } else {
          n2 = M^2 * (sd1^2 / kappa + sd2^2) / (HA_H0)^2
        }
        
        
        if(abs(n2 - n2.0) < tol){conv <- TRUE}
        n2.0 <- (n2 + n2.0)/2
        i <- i + 1
      }
      ifelse(df > 0, n2 <- n2.0, n2 <- NA)
      n1 <- n2 * kappa
    }
    
    # power
    if (is.null(power)) {
      
      HA_H0 <- mu1 - mu2
      if(paired) {
        df <- n2 - 1
        lambda <- (HA_H0) / sqrt((sd1^2 + sd2^2 - 2 * sd1 * sd2 * paired.r) / n2)
      } else {
        n1 <- n2 * kappa
        ifelse(welch.df,
               df <- welch_df(sd1, sd2, n1, n2),
               df <- n1 + n2 - 2)
        lambda = (HA_H0) / sqrt(sd1^2 / n1 + sd2^2 / n2)
      }
      
      power <- 1 - pt(qt(alpha / 2, df = df, ncp = 0, lower.tail = FALSE), df = df, ncp = lambda) +
        pt(-qt(alpha / 2, df = df, ncp = 0, lower.tail = FALSE), df = df, ncp = lambda)
    }
    
  }
  else if (alternative == "greater" | alternative ==
           "less") {
    
    if (margin != 0)
      warning("`margin` argument is ignored", call. = FALSE)
    if (alternative == "greater" & (mu1 < mu2))
      stop("alternative = 'greater' but mu1 < mu2",
           call. = FALSE)
    if (alternative == "less" & (mu1 > mu2))
      stop("alternative = 'less' but mu1 > mu2",
           call. = FALSE)
    
    # sample size
    if (is.null(n2)) {
      HA_H0 <- mu1 - mu2
      beta <- 1 - power
      i <- 0
      tol <- 0.01
      conv <- FALSE
      n2.0 <- 30
      while(i <= 100 & conv == FALSE){
        
        if(paired) {
          df <- n2.0 - 1
        } else {
          n1.0 <- n2.0 * kappa
          ifelse(welch.df,
                 df <- welch_df(sd1, sd2, n1.0, n2.0),
                 df <- n1.0 + n2.0 - 2)
        }
        
        if(df <= 0 | is.infinite(df)){break}
        
        M <- qt(alpha, df = df, ncp = 0, lower.tail = FALSE) +
          qt(beta, df = df, ncp = 0, lower.tail = FALSE)
        
        if(paired) {
          n2 = M^2 * (sd1^2 + sd2^2 - 2 * sd1 * sd2 * paired.r) / (HA_H0)^2
        } else {
          n2 = M^2 * (sd1^2 / kappa + sd2^2) / (HA_H0)^2
        }
        
        if(abs(n2 - n2.0) < tol){conv <- TRUE}
        n2.0 <- (n2 + n2.0)/2
        i <- i + 1
      }
      ifelse(df > 0, n2 <- n2.0, n2 <- NA)
      n1 <- n2 * kappa
    }
    
    # power
    if (is.null(power)) {
      HA_H0 <- mu1 - mu2
      if(paired) {
        df <- n2 - 1
        lambda <- (HA_H0) / sqrt((sd1^2 + sd2^2 - 2 * sd1 * sd2 * paired.r) / n2)
      } else {
        n1 <- n2 * kappa
        ifelse(welch.df,
               df <- welch_df(sd1, sd2, n1, n2),
               df <- n1 + n2 - 2)
        lambda = (HA_H0) / sqrt(sd1^2 / n1 + sd2^2 / n2)
      }
      
      
      if(alternative == "less") lambda <- abs(lambda)
      power <- 1 - pt(qt(alpha, df = df, ncp = 0, lower.tail = FALSE), df = df, ncp = lambda)
    }
    
  }
  else if (alternative == "non-inferior" | alternative ==
           "superior") {
    
    # sample size
    if (is.null(n2)) {
      beta <- 1 - power
      HA_H0 <- mu1 - mu2 - margin
      i <- 0
      tol <- 0.01
      conv <- FALSE
      n2.0 <- 30
      while(i <= 100 & conv == FALSE){
        n1.0 <- n2.0 * kappa
        
        if(paired) {
          df <- n2.0 - 1
        } else {
          ifelse(welch.df,
                 df <- welch_df(sd1, sd2, n1.0, n2.0),
                 df <- n1.0 + n2.0 - 2)
        }
        
        if(df <= 0 | is.infinite(df)){break}
        M <- qt(alpha, df = df, ncp = 0, lower.tail = FALSE) +
          qt(beta, df = df, ncp = 0, lower.tail = FALSE)
        
        if(paired) {
          n2 = M^2 * (sd1^2 + sd2^2 - 2 * sd1 * sd2 * paired.r) / (HA_H0)^2
        } else {
          n2 = M^2 * (sd1^2 / kappa + sd2^2) / (HA_H0)^2
        }
        
        if(abs(n2 - n2.0) < tol){conv <- TRUE}
        n2.0 <- (n2 + n2.0)/2
        i <- i + 1
      }
      ifelse(df > 0, n2 <- n2.0, n2 <- NA)
      n1 <- n2 * kappa
    }
    
    # power
    if (is.null(power)) {
      HA_H0 <- mu1 - mu2 - margin
      n1 <- n2 * kappa
      if(paired) {
        df <- n2 - 1
        lambda <- (HA_H0) / sqrt((sd1^2 + sd2^2 - 2 * sd1 * sd2 * paired.r) / n2)
      } else {
        ifelse(welch.df,
               df <- welch_df(sd1, sd2, n1, n2),
               df <- n1 + n2 - 2)
        lambda = (HA_H0) / sqrt(sd1^2 / n1 + sd2^2 / n2)
      }
      
      if(alternative == "non-inferior") lambda <- abs(lambda)
      power <- 1 - pt(qt(alpha, df = df, ncp = 0, lower.tail = FALSE), df = df, ncp = lambda)
    }
    
  }
  else if (alternative == "equivalent") {
    
    # sample size
    if (is.null(n2)) {
      beta <- 1 - power
      HA_H0 <- abs(mu1 - mu2) - margin
      i <- 0
      tol <- 0.01
      conv <- FALSE
      n2.0 <- 30
      while(i <= 100 & conv == FALSE){
        n1.0 <- n2.0 * kappa
        
        if(paired) {
          df <- n2.0 - 1
        } else {
          ifelse(welch.df,
                 df <- welch_df(sd1, sd2, n1.0, n2.0),
                 df <- n1.0 + n2.0 - 2)
        }
        
        if(df <= 0 | is.infinite(df)){break}
        M <- qt(alpha, df = df, ncp = 0, lower.tail = FALSE) +
          qt(beta / 2, df = df, ncp = 0, lower.tail = FALSE)
        
        if(paired) {
          n2 = M^2 * (sd1^2 + sd2^2 - 2 * sd1 * sd2 * paired.r) / (HA_H0)^2
        } else {
          n2 = M^2 * (sd1^2 / kappa + sd2^2) / (HA_H0)^2
        }
        
        if(abs(n2 - n2.0) < tol){conv <- TRUE}
        n2.0 <- (n2 + n2.0)/2
        i <- i + 1
      }
      ifelse(df > 0, n2 <- n2.0, n2 <- NA)
      n1 <- n2 * kappa
    }
    
    # power
    if (is.null(power)) {
      HA_H0 <- abs(mu1 - mu2) - margin
      if(paired) {
        df <- n2 - 1
        lambda <- (HA_H0) / sqrt((sd1^2 + sd2^2 - 2 * sd1 * sd2 * paired.r) / n2)
      } else {
        n1 <- n2 * kappa
        ifelse(welch.df,
               df <- welch_df(sd1, sd2, n1, n2),
               df <- n1 + n2 - 2)
        lambda = (HA_H0) / sqrt(sd1^2 / n1 + sd2^2 / n2)
      }
      
      power <- 2 * (1 - pt(qt(alpha, df = df, ncp = 0, lower.tail = FALSE), df = df, ncp = lambda)) - 1
      if(power < 0) stop("design is not feasible", call. = FALSE)
    }
    
  }
  else {
    stop("`alternative` should be in c('not equal', 'greater', 'less', 'equivalent', 'non-inferior', 'superior')",
         call. = FALSE)
  }
  
  
  if(paired) {
    n2 <- ceiling(n2)
    ncp <- (HA_H0) / sqrt((sd1^2 + sd2^2 - 2 * sd1 * sd2 * paired.r) / n2)
    df <- n2 - 1
  } else {
    n1 <- ceiling(n1)
    n2 <- ceiling(n2)
    ncp <- (HA_H0) / sqrt(sd1^2 / n1 + sd2^2 / n2)
    ifelse(welch.df,
           df <- welch_df(sd1, sd2, n1, n2),
           df <- n1 + n2 - 2)
  }
  ifelse(paired, n <- n2, n <- c(n1 = n1, n2 = n2))
  
  hypothesis <- alternative
  
  if(verbose) {
    cat(ifelse(paired,
               " Difference between Two means \n (Paired Samples t Test) \n",
               " Difference between Two means \n (Independent Samples t Test) \n"),
        switch(hypothesis,
               `not equal` = "H0: mu1 = mu2 \n HA: mu1 != mu2 \n",
               `greater` = "H0: mu1 = mu2 \n HA: mu1 > mu2 \n",
               `less` = "H0: mu1 = mu2 \n HA: mu1 < mu2 \n",
               `non-inferior` = "H0: mu1 - mu2 <= margin \n HA: mu1 - mu2 > margin \n",
               `superior` = "H0: mu1 - mu2 <= margin \n HA: mu1 - mu2 > margin \n",
               `equivalent` = "H0: |mu1 - mu2| >= margin \n HA: |mu1 - mu2| < margin \n"),
        "------------------------------ \n",
        " Statistical power =", round(power, 3), "\n",
        if(paired) {
          c(" n =", ceiling(n))
        } else {
          c(" n1 =", ceiling(n1), "\n  n2 =", ceiling(n2))
        }, "\n",
        "------------------------------ \n",
        "Alternative =", dQuote(alternative),"\n",
        "Degrees of freedom =", round(df, 2), "\n",
        "Non-centrality parameter =", round(ncp, 3), "\n",
        "Type I error rate =", round(alpha, 3), "\n",
        "Type II error rate =", round(1 - power, 3), "\n")
  }
  
  invisible(structure(list(parms = list(mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2, kappa = kappa, welch.df = welch.df,
                                        paired = paired, paired.r = paired.r,
                                        alpha = alpha, margin = margin, alternative = alternative, verbose = verbose),
                           test = "t",
                           df = df,
                           ncp =  ncp,
                           power = power,
                           n = n),
                      class = c("pwrss", "t", "2means")))
}