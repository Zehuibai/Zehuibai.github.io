### Setup file to be sourced at beginning of each new chapter

### Load libraries

# install.packages(c("tidyverse", "mmrm", "emmeans", "nlme", "gtsummary", "rbmi", "gt"))
# ## get the file path
# rstudioapi::getSourceEditorContext()$path
# ## get the wd path
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# getwd()


library(tidyverse)
library(mmrm)
library(emmeans)
library(gtsummary)
library(rbmi)
library(gt)
library("Hmisc")
library("naniar")
library("sjPlot")
library("kableExtra")

### Load data

### This is the small example dataset, created from the low dropout
all <- haven::read_sas("data/all2.sas7bdat")
colnames(all) <- tolower(colnames(all))

drp_grp <- all %>%
  dplyr::filter(!is.na(chgdrop)) %>%
  dplyr::group_by(subject) %>%
  dplyr::summarise(
    dropout_grp = dplyr::case_when(
      max(time) == 1 ~ "Week 2 Dropout",
      max(time) == 2 ~ "Week 4 Dropout",
      max(time) == 3 ~ "Completer"
    ),
    .groups = "drop"
  )

all2 <- all %>%
  dplyr::full_join(., drp_grp, by = "subject") %>%
  dplyr::mutate(
    aval = change + basval,
    avisit = dplyr::recode(as.character(time), "1" = "Week 2", "2" = "Week 4", "3" = "Week 8"),
    avisit = factor(avisit, levels = c("Week 2", "Week 4", "Week 8")),
    week = ifelse(time==3,8,time*2),
    subject = factor(subject),
    group = factor(trt, levels = 1:2, labels = c("Arm 1","Arm 2"))
  )

rm(drp_grp)

### Larger dataset with low dropout rate
### To be deleted
# low <- haven::read_sas("data/low2.sas7bdat") %>%
#   dplyr::mutate(
#     aval = change + basval
#   )
#
# colnames(low) <- tolower(colnames(low))
### To be deleted

high <- haven::read_sas("data/high2.sas7bdat")
colnames(high) <- tolower(colnames(high))

high2 <- high %>% dplyr::group_by(patient) %>%
dplyr::mutate( drop=max(week), .groups = "drop")

high2 <- high2 %>%
dplyr::mutate(
  aval = change + basval,
  group = factor(trt, levels = 1:2, labels = c("Arm 1","Arm 2")),
  avisit = dplyr::recode(as.character(week), "1" = "Week 1", "2" = "Week 2", "4" = "Week 4", "6" = "Week 6", "8" = "Week 8"),
  avisit = factor(avisit, levels = c("Week 1", "Week 2", "Week 4", "Week 6", "Week 8")),
  dropout_grp = dplyr::recode(as.character(drop), "1" = "Week 1 Drop", "2" = "Week 2 Drop", "4" = "Week 4 Drop", "6" = "Week 6 Drop", "8" = "Completer"),
  dropout_grp = factor(dropout_grp, levels=c("Week 1 Drop","Week 2 Drop","Week 4 Drop","Week 6 Drop","Completer")),
  subject = as.factor(patient)
)
