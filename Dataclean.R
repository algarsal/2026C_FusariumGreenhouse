
F1 <- read.csv("data/foa4_gh_data (2).csv", row.names = NULL)

F1 <- F1[-c(7, 11)]



save(F1, file = "data/fusarium.Rdata")


#For every new scriot 
load("data/fusarium.Rdata")
