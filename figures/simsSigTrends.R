res<-readRDS("trendSimResults.rds") %>% data.table()
setkey(res,nYears,nSites,simNum)
res[,trend:=rep(log(c(-0.01,-0.025,-0.05)+1),each=9)]
res[,converged:=all(rHat<1.1),by=.(nSites,nYears,trend,simNum)]

badOnes<-res[parameter=="mu",.(didNotCoverge=sum(!converged[!is.na(converged)])/sum(!is.na(converged)),
                              didNotInitialize=sum(is.na(trueValue))/length(trueValue)),
             by=.(nSites,nYears,trend)]

res<-res[!is.na(trueValue)&converged==T]

res[,diffFromTrue:=q2.5>trueValue|q97.5<trueValue]
propSig<-res[parameter=="trend",.(propSig95=sum(q97.5<0)/length(q97.5),
                                  propSig90=sum(q95<0)/length(q95),
                                  propSig80=sum(q90<0)/length(q95),
                                  n=.N),by=.(nYears,nSites,trend)]

cols=c(rgb(1,0,0,0.7),rgb(0,1,0,0.7),rgb(0,0,1,0.7))
tiff.par("figures/simsSigTrends.tif",mfrow=c(3,1),height=6,width=3.5)
plot(propSig95~trend,pch=19,col=cols[match(nYears,c(5,10,20))],cex=log(nSites/10),data=propSig,
     xlab="true trend",ylab="proportion overlapping 0 (95% CI)")
legend(-0.015,0.9,c("5","10","20"),col=cols,pch=19,pt.cex=2,title="Years",bty='n')
legend(-0.02,0.9,c("50","100","150"),pt.cex=log(c(50,100,150)/10),pch=19,
       col='black',title="Sites",bty='n')
plot(propSig90~trend,pch=19,col=cols[match(nYears,c(5,10,20))],cex=log(nSites/10),data=propSig,
     xlab="true trend",ylab="proportion overlapping 0 (90% CI)")
plot(propSig80~trend,pch=19,col=cols[match(nYears,c(5,10,20))],cex=log(nSites/10),data=propSig,
     xlab="true trend",ylab="proportion overlapping 0 (80% CI)")
dev.off()

