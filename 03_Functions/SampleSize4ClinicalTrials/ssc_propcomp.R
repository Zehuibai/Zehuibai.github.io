#' @title Sample Size Calculation for the Comparison of Proportions in Phase III Clinical Trials
#'
#' @description This function aims to calculate sample size for the comparison of proportions
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
#' @param ratio The ratio between the number of subjects in the treatment arm and that in the control arm.
#' @param alpha Type I error rate
#' @param power Statistical power of the test (1-type II error rate)
#' @param p1 The true mean response rate of the treatment arm
#' @param p2 The true mean response rate of the control arm
#' @param delta The prespecified superiority, non-inferiority or equivalence margin
#'
#' @return samplesize
#'
#' @usage ssc_propcomp(design, ratio, alpha, power, p1, p2, delta)
#'
#' @importFrom stats qnorm
#'
#' @examples
#' ##The comparison of proportions, an equivalence trial and the equivalence margin is 0.2
#' ssc_propcomp(design = 4L, ratio = 1, alpha = 0.05, power = 0.8, p1 = 0.75, p2 = 0.80, delta = 0.2)
#' @references
#' Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 2nd Ed. Chapman & Hall/CRC Biostatistics Series.
#'
#' Yin, G. 2012. Clinical Trial Design: Bayesian and Frequentist Adaptive Methods. John Wiley & Sons.
#'
#' @export

##Sample size calculation for the comparison of proportions
ssc_propcomp<-function(design = c(1L,2L,3L,4L), ratio = 1, alpha = 0.05, power = 0.8, p1 = NULL, p2 = NULL, delta = NULL){
  if (!(design %in% 1L:4L))
    stop("Unrecognized study design, 1: Test for equility, 2: Superiority trial,
         3: Non-inferiority trial, 4: Equivalence trial")
  ##Numerator
  if (design==1L) {
    numerator <- (abs(qnorm(alpha/2)) + abs(qnorm(1- power)))^2 * (p1*(1 - p1)/ratio + p2*(1 - p2))
  }
  if (design%in%2L:3L) {
    numerator <- (abs(qnorm(alpha)) + abs(qnorm(1- power)))^2 * (p1*(1 - p1)/ratio + p2*(1 - p2))
  }
  if (design==4L) {
    numerator <- (abs(qnorm(alpha)) + abs(qnorm((1- power)/2)))^2 * (p1*(1 - p1)/ratio + p2*(1 - p2))
  }


  ##Denominator
  if (design == 1L) {
    denom <- (p1 - p2)^2
  }
  if (design %in% 2L:3L)
    denom <- ((p1 - p2) - delta)^2
  if (design == 4L)
    denom <- (delta - abs(p1 - p2))^2

  ##n4 means number of sujects in the control arm
  n4<-ceiling(numerator/denom)
  n3<-ratio*n4

  ##Calculate the sample size
  samplesize<- data.frame(Treatment = n3, Control = n4)
  return(samplesize)
}
