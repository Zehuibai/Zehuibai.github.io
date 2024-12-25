## -------------------------------
##
## Script name: 
##
## Purpose of script:
##
## Author: Zehui Bai
##
## Date Created: 202X-XX-XX
##
## -------------------------------
##
## Notes:
##   
##
## -------------------------------


# <!-- ---------------------------------------------------------------------- -->
# <!--                        1. Basic system settings                        -->
# <!-- ---------------------------------------------------------------------- -->

## get the R Version
paste(R.Version()[c("major", "minor")], collapse = ".")

## get the file path
rstudioapi::getSourceEditorContext()$path

## get the wd path
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()



# <!-- ---------------------------------------------------------------------- -->
# <!--                    2. load the required packages                       -->
# <!-- ---------------------------------------------------------------------- --> 

packages<-c("tidyverse", "pharmaverse",
            "haven", "readxl", "writexl",
            "kableExtra", "clinUtils")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)


## Load all r functions 
# source_code_dir <- "C:/Users/baiz/Downloads/Data-Analyst-with-R/00 R Function/ZB Function/"  
# file_path_vec <- list.files(source_code_dir, full.names = T)
# for(f_path in file_path_vec){source(f_path)}


# <!-- ---------------------------------------------------------------------- -->
# <!--                    4. Check Folder if Needed                           -->
# <!-- ---------------------------------------------------------------------- -->

## Define the global variable
folder_check <- "No"

if (folder_check == "yes") {
  # List of directories to be checked and potentially created
  directories <- c("Output", "Check", "Log")
  # Loop through each directory in the list
  for (dir in directories) {
    # Check if the directory exists
    if (!dir.exists(dir)) {
      # Create the directory if it does not exist
      dir.create(dir)
      print(paste("Directory created:", dir))
    } else {
      print(paste("Directory already exists:", dir))
    }
  }
}


# <!-- ---------------------------------------------------------------------- -->
# <!--                       4. Load the datasets                             -->
# <!-- ---------------------------------------------------------------------- -->

# <!-- ---------------------------- -->
# <!-- --4.1 Import csv data ------ -->
# <!-- ---------------------------- -->

# pfad <- "~/Desktop/SASUniversityEdition/myfolders/Daten"
# mydata1 <- read.csv(file.path(pfad, "yourcsv_data.csv"), 
#                     sep=";", 
#                     header=TRUE)   

## Import all csv data from folder
# list_csv_files <- list.files(path = "./csvfolder/")
# do.call(rbind, lapply(list_csv_files, function(x) read.csv(x, stringsAsFactors = FALSE)))

# <!-- ---------------------------- -->
# <!-- --4.2 Import xlsx data ----- -->
# <!-- ---------------------------- -->

# library(readxl)
# mydata2 <- read_excel("C:/Users/zbai/Documents/GitHub/R-Projects/SAS/Yimeng/results-text.xlsx")

# <!-- ---------------------------- -->
# <!-- --4.3 Import sas7dbat data - -->
# <!-- ---------------------------- -->

# library(sas7bdat)
# mydata3 <- read.sas7bdat("~/Desktop/SASUniversityEdition/myfolders/Daten/uis.sas7bdat")

## Import all sas7dbat data from SASfolder
# ZB.import.sas.folder("./SASfolder/")




# <!-- ---------------------------------------------------------------------- -->
# <!--                           Start Programming                            -->
# <!-- ---------------------------------------------------------------------- -->

 

# Load data
d <- read_excel('./goniometer.xlsx') |>
  mutate(time=as.numeric(time))

# Bland-Altman plot - Inter-goniometer agreement
df<-d |>
  pivot_wider(names_from = `rater`, values_from = `y`) |>
  mutate(
    time = paste("T", time),
    diff = manual - electro,
    avg = (manual + electro) / 2)
df2 <- df |>
  group_by(time) |>
  summarise(
    g_avg = mean(avg),
    g_diff = mean(diff),
    sd_diff = sd(diff),
    n = n()
  )
df2  
df <- df |>
  left_join(df2)

ggplot(data=df, aes(x = avg, y = diff, col = ifelse(diff>(g_diff+1.96*sd_diff) | diff<(g_diff-1.96*sd_diff), "1", "0"))) +
  geom_point() +
  geom_rug(col="grey") +
  # orange and blue colors do not show legend  
  scale_color_manual(values = c("black", "red")) +
  geom_hline(aes(yintercept = g_diff)) +
  geom_hline(aes(yintercept = g_diff+1.96*sd_diff), linetype = "dashed") +
  geom_hline(aes(yintercept = g_diff-1.96*sd_diff), linetype = "dashed") +
  #geom_smooth(method = "lm", se = F) +
  labs(
    title = "Bland-Altman plot - Inter-goniometer agreement (electro vs. manual)",
    x = "Average",
    y = "Difference",
    caption = "Solid line represents average of differences. \nDashed lines represent limits of agreement (±1.96xSD). \nData outside 95% limits of agreement are labeled"
  ) +
  facet_wrap(time~.) +
  #plot data outside g_diff+1.96*sd_diff
  geom_text(data=df, aes(label = ifelse(diff>(g_diff+1.96*sd_diff) | diff<(g_diff-1.96*sd_diff), id, NA), hjust = -0.5)) +
  scale_y_continuous(limits = c(-10, 10)) +
  theme_bw() +
  theme(legend.position = "none")
# save the plot
ggsave("./inter-goniometer-agreement.png", width = 10, height = 5)

# Bland-Altman plot - Intra-goniometer agreement
df<-NULL
for (i in 1:3){
  c <- combn(1:3,2)[,i]
  d1 <- d |>
    filter(time %in% c) |>
    pivot_wider(names_from = `time`, values_from = `y`) 
  
  names(d1) <- c("id", "goniometer", "A", "B") 
  
  d2 <- d1 |>
    mutate(
      diff = A - B,
      avg = (A + B) / 2,
      comparison=paste("T",c[1],"vs", "T", c[2]))
  
  df <- rbind(df, d2)
}

df2 <- df |>
  group_by(goniometer, comparison) |>
  summarise(
    g_avg = mean(avg),
    g_diff = mean(diff),
    sd_diff = sd(diff),
    n = n()
  )

df <- df |>
  left_join(df2)

ggplot(data=df, aes(x = avg, y = diff, col = ifelse(diff>(g_diff+1.96*sd_diff) | diff<(g_diff-1.96*sd_diff), "1", "0"))) +
  geom_point() +
  geom_rug(col="grey") +
  # orange and blue colors do not show legend  
  scale_color_manual(values = c("black", "red")) +
  geom_hline(aes(yintercept = g_diff)) +
  geom_hline(aes(yintercept = g_diff+1.96*sd_diff), linetype = "dashed") +
  geom_hline(aes(yintercept = g_diff-1.96*sd_diff), linetype = "dashed") +
  #geom_smooth(method = "lm", se = F) +
  labs(
    title = "Bland-Altman plot - Intra-goniometer agreement",
    x = "Average",
    y = "Difference",
    caption = "Solid line represents average of differences. \nDashed lines represent limits of agreement (±1.96xSD). \nData outside 95% limits of agreement are labeled"
  ) +
  # add text to geom_hline
  facet_wrap(goniometer~comparison) +
  geom_text(aes(label = ifelse(diff>(g_diff+1.96*sd_diff) | diff<(g_diff-1.96*sd_diff), id, NA), hjust = -0.5)) +
  theme_bw() +
  theme(legend.position = "none")

# save the plot
ggsave("./intra-goniometer-agreement.png", width = 9, height = 6)

# <!-- ---------------------------------------------------------------------- -->
# <!--                             End Programming                            -->
# <!-- ---------------------------------------------------------------------- -->


