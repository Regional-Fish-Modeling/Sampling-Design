library(dplyr)
library(data.table)

res<-readRDS("yetiSims/trendSimResults1.rds") %>% data.table()
setkey(res,nYears,nSites,simNum)
res[,trend:=rep(log(c(-0.01,-0.025,-0.05)+1),each=9)]
res[,converged:=all(rHat<1.1),by=.(nSites,nYears,trend,simNum)]

badOnes<-res[parameter=="mu",.(didNotConverge=sum(!converged[!is.na(converged)])/sum(!is.na(converged)),
                               didNotInitialize=sum(is.na(trueValue))/length(trueValue)),
             by=.(nSites,nYears,trend)]

res<-res[!is.na(trueValue)&converged==T]

res[,diffFromTrue:=q2.5>trueValue|q97.5<trueValue]
propSig<-res[parameter=="trend",.(propSig95=sum(q97.5<0)/length(q97.5),
                                  propSig90=sum(q95<0)/length(q95),
                                  propSig80=sum(q90<0)/length(q95),
                                  n=.N),by=.(nYears,nSites,trend)]
propSig[,didNotConverge:=badOnes$didNotConverge]
propSig[,numToRun:=ceiling((100-n)/(1-didNotConverge))]
toRun<-propSig[,.(nYears,nSites,trend,numToRun)]

toRun<-rbind(toRun,data.table(nYears=rep(c(5,10,20),3),
                              nSites=rep(25,9),
                              trend=rep(log(c(-0.01,-0.025,-0.05)+1),each=3),
                              numToRun=rep(120,9)))

toRun<-toRun[rep(1:nrow(toRun),numToRun)]
toRun[,whichSim:=round((1:nrow(toRun))-5.1,-1)/10+1]
saveRDS(toRun,"yetiSims/simsToDo.rds")
