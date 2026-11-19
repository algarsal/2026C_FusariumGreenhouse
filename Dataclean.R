library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(multcompView)
library(emmeans)
library(broom)

# Load data
file_path <- "C:/Users/User/Downloads/Foa4_23_virulence_data_.xlsx"
F1 <- read_excel(file_path, sheet = 1)

