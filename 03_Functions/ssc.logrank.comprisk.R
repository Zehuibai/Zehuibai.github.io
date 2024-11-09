 
## Logrank Tests Accounting for Competing Risks
## https://www.ncss.com/wp-content/themes/ncss/pdf/Procedures/PASS/Logrank_Tests_Accounting_for_Competing_Risks.pdf

## S_ev1: Survival probability for the event of interest in group 1
## S_ev2: Survival probability for the event of interest in group 2
## S_cr: Competing risks 0.05
## P1: Proportion of subjects in group 1, the control group
## P2: Proportion of subjects in group 2, the treatment group
## R: Accrual Time 4 years
## T: Total time of trial 7 years

ssc.logrank.comprisk <- function(S_ev1, 
                                 S_ev2, 
                                 S_cr, 
                                 P1, 
                                 P2, 
                                 T, 
                                 R,
                                 alpha=0.05, 
                                 power=0.8,
                                 pr=TRUE){
  
  ## Quantile from standard normal distribution
  z_power <- qnorm(power, mean = 0, sd = 1)
  z_alpha <- qnorm(1-alpha, mean = 0, sd = 1)
  
  ## Hazard Ratio
  H_ev1 <- -log(S_ev1)/R
  H_ev2 <- -log(S_ev2)/R
  H_cr <- -log(S_cr)/R
  HR <- H_ev2/H_ev1
  
 
  ## Probability of observing event of interest in the group
  P_ev1 <- (H_ev1/(H_ev1+H_cr))*(1-((exp(-(T-R)*(H_ev1+H_cr))-exp(-T*(H_ev1+H_cr)))/(R*(H_ev1+H_cr))))
  P_ev2 <- (H_ev2/(H_ev2+H_cr))*(1-((exp(-(T-R)*(H_ev2+H_cr))-exp(-T*(H_ev2+H_cr)))/(R*(H_ev2+H_cr))))
  P_ev <- P1*P_ev1+P2*P_ev2
  
  ## Overall sample size required
  N <- (1/(P1*P2*P_ev))*(((z_alpha+z_power)/log(HR))^2)

  ## total number of events for the risk factor of interest
  E <- N*P_ev
  
  if(pr){
    cat("\nSample Size Calculation using Logrank Tests Accounting for Competing Risks\n")
    cat("Alpha",alpha,"\n")
    cat("Power",power*100,"%\n")
    cat("\n")
    cat("Accrual time of survival rate observed:",R,"years\n")
    cat("Total time of tria:",T,"years\n")
    cat("Follow-Up Time:",T-R,"years\n")
    cat("\n")
    cat("Survival probability for the event of interest in group 1:",S_ev1,"\n")
    cat("Survival probability for the event of interest in group 2:",S_ev2,"\n")
    cat("Hazard Ration:",HR,"\n")
    cat("\n")
    cat("Competing risks probability:",S_cr,"\n")
    cat("\n")
    cat("Proportion of subjects in group 1:",P1,"\n")
    cat("Proportion of subjects in group 2:",P2,"\n")
    cat("\n")
    cat("The probability of observing the event of interest in a subject during the study for the group 1:", P_ev1,"\n")
    cat("The probability of observing the event of interest in a subject during the study for the group 2:", P_ev2,"\n")
    cat("\n")
    cat("The number of events required for the group 1:",ceiling(N*P1*P_ev1),"\n")
    cat("The number of events required for the group 2:",ceiling(N*P2*P_ev2),"\n")
    cat("The total number of events required for the study:",ceiling(N*P1*P_ev1) + ceiling(N*P2*P_ev2),"\n")
    cat("\n")
    cat("The sample sizes for the group 1:",ceiling(N*P1),"\n")
    cat("The sample sizes for the group 2:",ceiling(N*P2),"\n")
    cat("The total sample size of both groups combined:",ceiling(N*P1)+ ceiling(N*P2),"\n")
  }else{
    sample <- c(ceiling(N*P1*P_ev1), 
                ceiling(N*P2*P_ev2),
                ceiling(N*P1*P_ev1) + ceiling(N*P2*P_ev2),
                ceiling(N*P1), 
                ceiling(N*P2),
                ceiling(N*P1)+ ceiling(N*P2))
    names(sample) <- c("N1_Event",
                       "N2_Event",
                       "Total Event", 
                       "N1_Patient",
                       "N2_Patient",
                       "Total Patient")
    return(sample)
  }
}
 
## Two sided 0.05
# ssc.logrank.comprisk(S_ev1 = 0.50, 
#                      S_ev2=0.7071068,
#                      S_cr=0.4, 
#                      P1=0.5, 
#                      P2=0.5, 
#                      T=5, 
#                      R=3,
#                      alpha=0.025,
#                      power=0.6162274)
# 
# ssc.logrank.comprisk(S_ev1 = 0.50, 
#                      S_ev2=0.7071068,
#                      S_cr=0.4, 
#                      P1=0.5, 
#                      P2=0.5, 
#                      T=5, 
#                      R=3,
#                      alpha=0.025,
#                      power=0.6162274,
#                      pr=F)
