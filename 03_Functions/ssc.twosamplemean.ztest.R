## Sample size calculation for mean comparison
## Based on SampleSize4ClinicalTrials Package

# design
  # 1L:  Testing for equality
  # 2L:  Superiority trial
  # 3L:  Non-inferiority trial
  # 4L:  Equivalence trial
# ratio: The ratio between the number of subjects in the treatment arm and that in the control arm
# alpha: Type I error rate
# power: Statistical power of the test (1-type II error rate)
# sd:    The standard deviation of observed outcomes in both arms
# theta: The true mean difference between two arms
# delta: The prespecified superiority, non-inferiority or equivalence margin
 

ssc.twosamplemean.ztest<-function(design = c(1L,2L,3L,4L), 
                       ratio = 1, 
                       alpha = 0.05, 
                       power = 0.8, 
                       sd = NULL,
                       theta = NULL, 
                       delta = NULL,
                       pr=TRUE){
  
  if (!(design %in% 1L:4L))
    stop("Unrecognized study design, 1: Test for equility, 2: Superiority trial,
         3: Non-inferiority trial, 4: Equivalence trial")
  
  ##Numerator
  if (design==1L) {
    numerator <- (abs(qnorm(alpha/2)) + abs(qnorm(1- power)))^2 * sd^2 * (1 + 1/ratio)
    Test_Design <- "Test for equilit"
  }
  ##Note that we unify the superiority and non-inferiority trials in this version (refer to Chow et al.)
  ##In the last version, there was no superiority margin and the non-inferiority margin was positive (refer to Yin).
  if (design %in% 2L:3L) {
    numerator <- (abs(qnorm(alpha)) + abs(qnorm(1 - power)))^2 * sd^2 * (1 + 1/ratio)
    if(design==2L){Test_Design <- "Test for superiority"}
    if(design==3L){Test_Design <- "Test for non-inferiority"}
  }
  if (design == 4L) {
    numerator <- (abs(qnorm(alpha)) + abs(qnorm((1 - power)/2)))^2 * sd^2 * (1 + 1/ratio)
    Test_Design <- "Test for equivalence"
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
  n2<-numerator/denom
  n1<-ratio*n2

  ##Calculate the sample size
  samplesize<- data.frame(Treatment = n1, Control = n2)
  SS_name <- c("Treatment Group","Control Group")
  names(samplesize) <- SS_name
  
  if(pr){
    cat("\nSample Size Calculation for Clinical Trials:", Test_Design ,"\n")
    cat("Alpha",alpha,"\n")
    cat("Power",power*100,"%\n")
    cat("Ratio between subjects in the Treatment Group and in the Control Group:",ratio,"\n")
    cat("Standard deviation:",sd,"\n")
    cat("Mean difference between two arms:",theta,"\n")
    if (design==2L) {cat("Superiority Margin",delta,"\n")}
    if (design==3L) {cat("Non-inferiority margin",delta,"\n")}
    if (design==4L) {cat("Equivalence margin",delta,"\n")}
    cat("\n")
    print(samplesize) 
  }else{
    return(samplesize)
  }
}

 
# ssc.twosamplemean.ztest(design = 1L, ratio = 1, alpha = 0.05, power = 0.8, sd = 0.1, theta = 0.2)
# ssc.twosamplemean.ztest(design = 2L, ratio = 1, alpha = 0.05, power = 0.8, sd = 0.1, theta = 0, delta =  0.05)
# ssc.twosamplemean.ztest(design = 3L, ratio = 1, alpha = 0.05, power = 0.8, sd = 0.1, theta = 0, delta = -0.05)
# ssc.twosamplemean.ztest(design = 4L, ratio = 1, alpha = 0.05, power = 0.8, sd = 0.1, theta = 0, delta =  0.05)
# 
# ssc.twosamplemean.ztest(design = 3L, ratio = 2, alpha = 0.05, power = 0.8, sd = 0.1, theta = 0, delta = -0.05, pr=F)

