library(data.table)
library(dplyr)
library(plotHacks)

res<-readRDS("yetiSims/trendSimResults.rds") %>%
  data.table() %>%
  setkey(nYears,nSites,simNum) %>%
  .[,trend:=trueValue[which(parameter=="trend")],by=simNum]
          

setkey(res,nYears,nSites,simNum)
res[,converged:=all(rHat<1.1),by=.(nSites,nYears,stage,trend,simNum)]

badOnes<-res[parameter=="mu",.(didNotCoverge=sum(!converged[!is.na(converged)])/sum(!is.na(converged)),
                              didNotInitialize=sum(is.na(trueValue))/length(trueValue)),
             by=.(nSites,nYears,trend)]

res<-res[!is.na(trueValue)&converged==T]

res[,diffFromTrue:=q2.5>trueValue|q97.5<trueValue]
propSig<-res[parameter=="trend",.(propSig95=sum(q97.5<0)/length(q97.5),
                                  propSig90=sum(q95<0)/length(q95),
                                  propSig80=sum(q90<0)/length(q95),
                                  n=.N),by=.(nYears,nSites,stage,trend)]

cols=c(rgb(1,0,0,0.6),rgb(0,1,0,0.6),rgb(0,0,1,0.6))
tiff.par("figures/simsSigTrendsAdult.tif",mfrow=c(3,1),height=6,width=3.5)
plot(propSig95~trend,pch=16,col=cols[match(nYears,c(5,10,20))],cex=log(nSites/2),data=propSig[stage=="adult"],
     xlab="true trend",ylab="proportion not overlapping 0 (95% CI)",ylim=c(0,1))
legend(-0.015,0.9,c("5","10","20"),col=cols,pch=19,pt.cex=2,title="Years",bty='n')
legend(-0.02,0.9,c("5","10","50","100"),pt.cex=log(c(25,50,100,150)/10),pch=16,
       col='black',title="Sites",bty='n')
plot(propSig90~trend,pch=19,col=cols[match(nYears,c(5,10,20))],cex=log(nSites/2),data=propSig[stage=="adult"],
     xlab="true trend",ylab="proportion not overlapping 0 (90% CI)",ylim=c(0,1))
plot(propSig80~trend,pch=19,col=cols[match(nYears,c(5,10,20))],cex=log(nSites/2),data=propSig[stage=="adult"],
     xlab="true trend",ylab="proportion not overlapping 0 (80% CI)",ylim=c(0,1))
dev.off()

tiff.par("figures/simsSigTrendsYoy.tif",mfrow=c(3,1),height=6,width=3.5)
plot(propSig95~trend,pch=16,col=cols[match(nYears,c(5,10,20))],cex=log(nSites/2),data=propSig[stage=="yoy"],
     xlab="true trend",ylab="proportion not overlapping 0 (95% CI)",ylim=c(0,1))
legend(-0.015,0.9,c("5","10","20"),col=cols,pch=19,pt.cex=2,title="Years",bty='n')
legend(-0.02,0.9,c("5","10","50","100"),pt.cex=log(c(25,50,100,150)/10),pch=16,
       col='black',title="Sites",bty='n')
plot(propSig90~trend,pch=19,col=cols[match(nYears,c(5,10,20))],cex=log(nSites/2),data=propSig[stage=="yoy"],
     xlab="true trend",ylab="proportion not overlapping 0 (90% CI)",ylim=c(0,1))
plot(propSig80~trend,pch=19,col=cols[match(nYears,c(5,10,20))],cex=log(nSites/2),data=propSig[stage=="yoy"],
     xlab="true trend",ylab="proportion not overlapping 0 (80% CI)",ylim=c(0,1))
dev.off()


