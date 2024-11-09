#' @title Sample Size Calculation for the Comparison of Means in Phase III Clinical Trials
#'
#' @description This function aims to calculate sample size for the comparison of means
#' in Phase III clinical trials.
#'
#' @param design The design of the clinical trials.
#'   \cr 1L
#'     \cr Testing for equality
#'   \cr 2L
#'     \cr Superiority trial
#'   \cr 3L
#'     \cr Non-inferiority trial
#'   \cr 4L
#'     \cr Equivalence trial.
#' @param ratio The ratio between the number of subjects in the treatment arm and that in the control arm
#' @param alpha Type I error rate
#' @param power Statistical power of the test (1-type II error rate)
#' @param sd The standard deviation of observed outcomes in both arms
#' @param theta The true mean difference between two arms
#' @param delta The prespecified superiority, non-inferiority or equivalence margin
#'
#' @return samplesize
#'
#' @usage ssc_meancomp(design, ratio, alpha, power, sd, theta, delta)
#'
#' @importFrom stats qnorm
#'
#' @examples
#' ##The comparison of means, a non-inferiority trial and the non-inferiority margin is -0.05
#' ssc_meancomp(design = 3L, ratio = 1, alpha = 0.05, power = 0.8, sd = 0.1, theta = 0, delta = -0.05)
#'
#' @references
#' Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 2nd Ed. Chapman & Hall/CRC Biostatistics Series.
#'
#' Yin, G. 2012. Clinical Trial Design: Bayesian and Frequentist Adaptive Methods. John Wiley & Sons.
#'
#' @export

##Sample size calculation for mean comparison
ssc_meancomp<-function(design = c(1L,2L,3L,4L), ratio = 1, alpha = 0.05, power = 0.8, sd = NULL,theta = NULL, delta = NULL){
  if (!(design %in% 1L:4L))
    stop("Unrecognized study design, 1: Test for equility, 2: Superiority trial,
         3: Non-inferiority trial, 4: Equivalence trial")
  ##Numerator
  if (design==1L) {
    numerator <- (abs(qnorm(alpha/2)) + abs(qnorm(1- power)))^2 * sd^2 * (1 + 1/ratio)
  }
  ##Note that we unify the superiority and non-inferiority trials in this version (refer to Chow et al.)
  ##In the last version, there was no superiority margin and the non-inferiority margin was positive (refer to Yin).
  if (design %in% 2L:3L) {
    numerator <- (abs(qnorm(alpha)) + abs(qnorm(1 - power)))^2 * sd^2 * (1 + 1/ratio)
  }
  if (design == 4L) {
    numerator <- (abs(qnorm(alpha)) + abs(qnorm((1 - power)/2)))^2 * sd^2 * (1 + 1/ratio)
  }

  ##Denominator
  if (design == 1L) {
    denom<-theta^2
  }
  if (design%in% 2L:3L) {
    denom<-(theta - delta)^2
  }
  if (design==4L) {
    denom<-(delta-abs(theta))^2
  }

  ##Sample size equation, n2 means number of subjects in the control arm
  n2<-ceiling(numerator/denom)
  n1<-ratio*n2

  ##Calculate the sample size
  samplesize<- data.frame(Treatment = n1, Control = n2)
  return(samplesize)
}
