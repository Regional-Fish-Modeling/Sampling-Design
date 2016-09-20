# for(b in 1:6){
#   assign(paste0("beta",b),out2$mean$b[b])
# }
beta1<-0.1638
beta2<- -0.3272
beta3<- -0.3010
beta4<- 0.0951
beta5<- -0.2475
beta6<- 0.1385

# mu<-out2$mean$mu
mu<- 2.8500

sdNoCov<-NA
sdCov<-NA

for(i in 1:1000){
fallPrcp<-rnorm(nYears)
fallTmean<-rnorm(nYears)
winterPrcp<-rnorm(nYears)
winterTmean<-rnorm(nYears)
springPrcp<-rnorm(nYears)
springTmean<-rnorm(nYears)

ranYear<-rnorm(nYears,0,sd.year)

y<-mu+
  beta1*fallPrcp+
  beta2*fallTmean+
  beta3*winterPrcp+
  beta4*winterTmean+
  beta5*springPrcp+
  beta6*springTmean+
  ranYear

a1<-lm(y~1)
a2<-lm(y~fallPrcp+fallTmean+winterPrcp+winterTmean+springPrcp+springTmean)

y1<-rnorm(1000,mu,summary(a1)$sigma)
y2<-coef(a2)[1]+
  coef(a2)[2]*rnorm(1000)+
  coef(a2)[3]*rnorm(1000)+
  coef(a2)[4]*rnorm(1000)+
  coef(a2)[5]*rnorm(1000)+
  coef(a2)[6]*rnorm(1000)+
  coef(a2)[7]*rnorm(1000)+
  rnorm(1000,mu,summary(a2)$sigma)
  
sdNoCov[i]<-sd(y1)
sdCov[i]<-sd(y2)
}

plot(density(sdNoCov))
points(density(sdCov),type='l',col="blue")

mean(sdNoCov)
mean(sdCov)
