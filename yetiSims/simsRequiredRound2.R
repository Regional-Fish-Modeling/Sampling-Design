res1<-readRDS("yetiSims/trendSimResults1.rds") %>% 
  data.table() %>%
  setkey(nYears,nSites,simNum) %>%
  .[,trend:=rep(log(c(-0.01,-0.025,-0.05)+1),each=9)]
res2<- readRDS("yetiSims/trendSimResults.rds") %>%
  data.table() %>%
  .[,trend:=trueValue[which(parameter=="trend")],by=simNum]
res<-rbind(res1,res2)
rm(list=c("res1","res2"))


setkey(res,nYears,nSites,simNum)
res[,converged:=all(rHat<1.1),by=.(nSites,nYears,trend,simNum)]

badOnes<-res[parameter=="mu"&!is.na(converged),
             .(didNotConverge=sum(!converged)/length(converged)),
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
toRun<-toRun[numToRun>0]
toRun[,numToRun:=numToRun*2]
toRun<-toRun[rep(1:nrow(toRun),numToRun)]
toRun[,whichSim:=round((1:nrow(toRun))-5.1,-1)/10+1]

saveRDS(toRun,"yetiSims/simsToDo.rds")
