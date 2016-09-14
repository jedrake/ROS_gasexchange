#---------------------------------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------- 
# Fits beta models for g1 and non-stomatal limitation, using predawn leaf water potential as the predictor.
#---------------------------------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------- 


#- get the leaf water potential data, average by date to combine with g1 and NSL
lwp <- return.gx.vwc.lwp()
lwp.m <- summaryBy(LWP+LWP.md~Species+Treat+gxDate,data=lwp,FUN=mean,keep.names=T,na.rm=T)

#- get the g1 fits
g1values <- returng1()
g1values2 <- merge(g1values,lwp.m,by.x=c("Species","Treat","Date"),by.y=c("Species","Treat","gxDate"))


#- get the NSL fits
NSL <- returnVcmaxa()
NSLpars <- summaryBy(NSL+Photo+Cond+Ci+TDR~Species+gxDate+Treat,data=NSL,FUN=mean,keep.names=T)
NSLpars2 <- merge(NSLpars,lwp.m,by.x=c("Species","Treat","gxDate"),by.y=c("Species","Treat","gxDate"))
NSLpars2 <- subset(NSLpars2,gxDate != as.Date("2012-10-31"))

#---------------------------------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------- 
# fit g1 first

dat <- g1values2
type="g1"
startlist = list(Xlow = 1, Xhigh=7, q = 2)
#- split into list of species
dat.l <- split(dat,dat$Species)

#- preallocate empty lists to catch output
fit.sp <- newdat <-  list()

#- loop over each element of the list, fit data
for (i in 1:length(dat.l)){
  dat.temp <- dat.l[[i]]
  dat.temp$Species <- factor(dat.temp$Species)
  dat.temp$LWPpos <- dat.temp$LWP+11
  
  #- fit the model
  if (type=="g1") dat.temp$Yval <- dat.temp$g1/max(dat.temp$g1)
  if (type=="NSL") dat.temp$Yval <- dat.temp$NSL
  fit.sp[[i]] <- nls(Yval ~ ((LWPpos-Xlow)/(Xhigh-Xlow))^q,start=startlist,data=dat.temp,algorithm="port",
                     lower=c(1,6,1.5),upper=c(5,11,6))
  
  
  # get predicted values and 95% confidence intervals by bootstrapping
  if (type=="g1") newdat[[i]] <- expand.grid(Species=levels(dat.temp$Species), LWPpos=seq(from=coef(fit.sp[[i]])["Xlow"]+0.01,to=max(dat.temp$LWPpos),length.out=99),lower=NA,upper=NA)
  if (type=="NSL") newdat[[i]] <- expand.grid(Species=levels(dat.temp$Species), LWPpos=seq(from=coef(fit.sp[[i]])["Xlow"]+0.01,to=coef(fit.sp[[i]])["Xhigh"],length.out=99),lower=NA,upper=NA)
  newdat[[i]]$wpred <- predict(fit.sp[[i]],newdat[[i]],level=0,se.fit=T)
  
  
  rm(b)
  b <- bootCase(fit.sp[[i]],B=300)
  for(j in 1:nrow(newdat[[i]])){
    rm(b02)
    b02 <- ((newdat[[i]]$LWPpos[j]-b[,"Xlow"])/(b[,"Xhigh"]-b[,"Xlow"]))^b[,"q"]
    
    newdat[[i]]$lower[j] <- unname(quantile(b02,na.rm=T,probs=c(0.025,0.975)))[1]
    newdat[[i]]$upper[j] <- unname(quantile(b02,na.rm=T,probs=c(0.025,0.975)))[2]
    
  }
  
}
newdatg1 <- newdat
fit.spg1.LWP <- fit.sp
rm(newdat)
rm(fit.sp)
g1_TDR_beta2 <- list(fit.spg1.LWP,newdatg1)
#---------------------------------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------- 




#---------------------------------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------- 
# fit NSL next 

dat <- NSLpars2
type="NSL"
startlist = list(Xlow = 1, Xhigh=7, q = 3)
#- split into list of species
dat.l <- split(dat,dat$Species)

#- preallocate empty lists to catch output
fit.sp <- newdat <-  list()

#- loop over each element of the list, fit data
for (i in 1:length(dat.l)){
  dat.temp <- dat.l[[i]]
  dat.temp$Species <- factor(dat.temp$Species)
  dat.temp$LWPpos <- dat.temp$LWP+11
  
  #- fit the model
  if (type=="g1") dat.temp$Yval <- dat.temp$g1/max(dat.temp$g1)
  if (type=="NSL") dat.temp$Yval <- dat.temp$NSL
  if (i <= 3) startlist = list(Xlow = 1, Xhigh=7, q = 3)
  
  if(i <= 3)fit.sp[[i]] <- nls(Yval ~ ((LWPpos-Xlow)/(Xhigh-Xlow))^q,start=startlist,data=dat.temp,algorithm="port",trace=F,
                     lower=c(1,7,1.1),upper=c(2,11,6))
  if(i==4) startlist = list(Xlow = 8.5, Xhigh=11, q = 2)
  
  if(i ==4)fit.sp[[i]] <- nls(Yval ~ ((LWPpos-Xlow)/(Xhigh-Xlow))^q,start=startlist,data=dat.temp,algorithm="port",trace=F,
                              lower=c(8.5,10,1.1),upper=c(9.5,11,6))
  
  
  # get predicted values and 95% confidence intervals by bootstrapping
  if (type=="g1") newdat[[i]] <- expand.grid(Species=levels(dat.temp$Species), LWPpos=seq(from=coef(fit.sp[[i]])["Xlow"]+0.01,to=max(dat.temp$LWPpos),length.out=99),lower=NA,upper=NA)
  if (type=="NSL") newdat[[i]] <- expand.grid(Species=levels(dat.temp$Species), LWPpos=seq(from=coef(fit.sp[[i]])["Xlow"],to=coef(fit.sp[[i]])["Xhigh"],length.out=99),lower=NA,upper=NA)
  newdat[[i]]$wpred <- predict(fit.sp[[i]],newdat[[i]],level=0,se.fit=T)
  
  
  rm(b)
  b <- bootCase(fit.sp[[i]],B=300)
  for(j in 1:nrow(newdat[[i]])){
    rm(b02)
    b02 <- ((newdat[[i]]$LWPpos[j]-b[,"Xlow"])/(b[,"Xhigh"]-b[,"Xlow"]))^b[,"q"]
    
    newdat[[i]]$lower[j] <- unname(quantile(b02,na.rm=T,probs=c(0.025,0.975)))[1]
    newdat[[i]]$upper[j] <- unname(quantile(b02,na.rm=T,probs=c(0.025,0.975)))[2]
    
  }
  
}

newdatNSL <- newdat
rm(newdat)
fit.spNSL.LWP <- fit.sp
rm(fit.sp)
NSL_TDR_beta2 <- list(fit.spNSL.LWP,newdatNSL)
#---------------------------------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------- 

