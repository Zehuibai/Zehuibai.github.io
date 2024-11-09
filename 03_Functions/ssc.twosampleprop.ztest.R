## Sample size calculation for proportions comparison
## Based on SampleSize4ClinicalTrials Package

# design
  # 1L: Testing for equality
  # 2L: Superiority trial
  # 3L: Non-inferiority trial
  # 4L: Equivalence trial
# ratio: The ratio between the number of subjects in the treatment arm and that in the control arm.
# alpha: Type I error rate
# power: Statistical power of the test (1-type II error rate)
# p1:    The true mean response rate of the treatment arm
# p2:    The true mean response rate of the control arm
# delta: The prespecified superiority, non-inferiority or equivalence margin

 
ssc.propcomp<-function(design = c(1L,2L,3L,4L), 
                       ratio = 1, 
                       alpha = 0.05, 
                       power = 0.8, 
                       p1 = NULL, 
                       p2 = NULL, 
                       delta = NULL,
                       pr=TRUE){
  
  if (!(design %in% 1L:4L))
    stop("Unrecognized study design, 1: Test for equility, 2: Superiority trial,
         3: Non-inferiority trial, 4: Equivalence trial")
  ##Numerator
  if (design==1L) {
    numerator <- (abs(qnorm(alpha/2)) + abs(qnorm(1- power)))^2 * (p1*(1 - p1)/ratio + p2*(1 - p2))
    Test_Design <- "Test for equilit"
  }
  if (design%in%2L:3L) {
    numerator <- (abs(qnorm(alpha)) + abs(qnorm(1- power)))^2 * (p1*(1 - p1)/ratio + p2*(1 - p2))
    if(design==2L){Test_Design <- "Test for superiority"}
    if(design==3L){Test_Design <- "Test for non-inferiority"}
  }
  if (design==4L) {
    numerator <- (abs(qnorm(alpha)) + abs(qnorm((1- power)/2)))^2 * (p1*(1 - p1)/ratio + p2*(1 - p2))
    Test_Design <- "Test for equivalence"
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
  SS_name <- c("Treatment Group","Control Group")
  names(samplesize) <- SS_name
  
  if(pr){
    cat("\nSample Size Calculation for Clinical Trials:", Test_Design ,"\n")
    cat("Alpha",alpha,"\n")
    cat("Power",power*100,"%\n")
    cat("Ratio between subjects in the Treatment Group and in the Control Group:",ratio,"\n")
    cat("Mean response rate of the Treatment Group:",p1,"\n")
    cat("Mean response rate of the Control Group:",p2,"\n")
    if (design==2L) {cat("Superiority Margin",delta,"\n")}
    if (design==3L) {cat("Non-inferiority margin",delta,"\n")}
    if (design==4L) {cat("Equivalence margin",delta,"\n")}
    cat("\n")
    print(samplesize) 
  }else{
    return(samplesize)
  }
}

# ssc.propcomp(design = 1L, ratio = 2, alpha = 0.05, power = 0.8, p1 = 0.75, p2 = 0.80)
# ssc.propcomp(design = 2L, ratio = 2, alpha = 0.05, power = 0.8, p1 = 0.75, p2 = 0.80, delta = 0.1)
# ssc.propcomp(design = 3L, ratio = 2, alpha = 0.05, power = 0.8, p1 = 0.75, p2 = 0.80, delta = -0.1)
# ssc.propcomp(design = 4L, ratio = 1, alpha = 0.05, power = 0.8, p1 = 0.75, p2 = 0.80, delta = 0.2)
# 
# ssc.propcomp(design = 4L, ratio = 1, alpha = 0.05, power = 0.8, p1 = 0.75, p2 = 0.80, delta = 0.2, pr=F)
