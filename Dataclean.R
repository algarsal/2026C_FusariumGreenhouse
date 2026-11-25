
F1 <- read.csv("data/foa4_gh_data (1).csv", row.names = NULL)

F1 <- F1[-c(7, 12)]



save(F1, file = "data/fusarium.Rdata")


#For every new script 
load("data/fusarium.Rdata")
