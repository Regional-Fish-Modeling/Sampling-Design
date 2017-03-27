
library(data.table)
#Settings to vary in simulations
nSites <- c(5, 10, 50, 100)
nYears <- c(5, 10, 20)
annualDec <- c(-0.01, -0.025, -0.05)
stage<-c("adult","yoy")
covariates<-c(FALSE,TRUE)

#Matrix of simulations for task assignment
simsToDo <- expand.grid(nSites, nYears, annualDec,stage,covariates)
names(simsToDo) <- c("nSites","nYears","trend","stage","covariates")

#Replicate for the number of iterations (100)
simsToDo<-simsToDo[rep(1:nrow(simsToDo),each=110),]
simsToDo$whichSim<-rep(1:(nrow(simsToDo)/10),each=10)
simsToDo<-data.table(simsToDo)
saveRDS(simsToDo,file="simsToDo.rds")

# demo run
#simsToDo<-simsToDo[c(1,2,11,12)]
#saveRDS(simsToDo,file="simsToDo.rds")
