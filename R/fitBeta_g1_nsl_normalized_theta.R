#-----
#- get the g1 fits
g1values <- returng1()

#- normalize theta between field capacity (0.14) and the wilting point (0.05).
#  Note that MOST of the data are at a volumetric water content that exceeds field capacity!
#  Also, negative x-values cause the Beta model fitting to fail (negative number raised to ~0.5 power)
g1values$TDR_norm <- (g1values$TDR-0.05)/(0.14-0.05)
g1values$TDR_norm2 <- g1values$TDR_norm+2

#- get the NSL fits
NSL <- returnVcmaxa()
NSLpars <- summaryBy(NSL+Photo+Cond+Ci+TDR~Species+gxDate+Treat,data=NSL,FUN=mean,keep.names=T)
NSLpars$TDR_norm <- (NSLpars$TDR-0.05)/(0.14-0.05)

#---------------------------------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------- 
# fit g1 first. Note that none of the fitting work yet. Ugh.

dat <- g1values
type="g1"
startlist = list(Xlow = 3.5, Xhigh=4, q = 0.2)
#- split into list of species
dat.l <- split(dat,dat$Species)

#- preallocate empty lists to catch output
fit.sp <- newdat <-  list()

#- loop over each element of the list, fit data
for (i in 1:length(dat.l)){
  dat.temp <- dat.l[[i]]
  dat.temp$Species <- factor(dat.temp$Species)
  
  #- fit the model
  if (type=="g1") dat.temp$Yval <- dat.temp$g1/max(dat.temp$g1)
  if (type=="NSL") dat.temp$Yval <- dat.temp$NSL
  fit.sp[[i]] <- nls(Yval ~ ((TDR_norm2-Xlow)/(Xhigh-Xlow))^q,start=startlist,data=dat.temp,algorithm="port",
                     lower=c(2,3,0.01),upper=c(3,6,6))
  
  
  # get predicted values and 95% confidence intervals by bootstrapping
  if (type=="g1") newdat[[i]] <- expand.grid(Species=levels(dat.temp$Species), TDR=seq(from=coef(fit.sp[[i]])["Xlow"]+0.01,to=max(dat.temp$TDR),length.out=99),lower=NA,upper=NA)
  if (type=="NSL") newdat[[i]] <- expand.grid(Species=levels(dat.temp$Species), TDR=seq(from=coef(fit.sp[[i]])["Xlow"]+0.01,to=coef(fit.sp[[i]])["Xhigh"],length.out=99),lower=NA,upper=NA)
  newdat[[i]]$wpred <- predict(fit.sp[[i]],newdat[[i]],level=0,se.fit=T)
  
  
  rm(b)
  b <- bootCase(fit.sp[[i]],B=300)
  for(j in 1:nrow(newdat[[i]])){
    rm(b02)
    b02 <- ((newdat[[i]]$TDR[j]-b[,"Xlow"])/(b[,"Xhigh"]-b[,"Xlow"]))^b[,"q"]
    
    newdat[[i]]$lower[j] <- unname(quantile(b02,probs=c(0.025,0.975)))[1]
    newdat[[i]]$upper[j] <- unname(quantile(b02,probs=c(0.025,0.975)))[2]
    
  }
  
}
newdatg1 <- newdat
fit.spg1 <- fit.sp
rm(newdat)
rm(fit.sp)
g1_TDR_beta <- list(fit.spg1,newdatg1)
#---------------------------------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------- 




#---------------------------------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------- 
# fit NSL next 

dat <- NSLpars
type="NSL"
startlist= list(Xlow = 0.0, Xhigh=0.1, q = 0.1)
#- split into list of species
dat.l <- split(dat,dat$Species)

#- preallocate empty lists to catch output
fit.sp <- newdat <-  list()

#- loop over each element of the list, fit data
for (i in 1:length(dat.l)){
  dat.temp <- subset(dat.l[[i]],TDR<0.3)
  if (i==4) dat.temp <- subset(dat.l[[i]],TDR<0.2)
  dat.temp$Species <- factor(dat.temp$Species)
  
  #- fit the model
  if (type=="g1") dat.temp$Yval <- dat.temp$g1/max(dat.temp$g1)
  if (type=="NSL") dat.temp$Yval <- dat.temp$NSL
  fit.sp[[i]] <- nls(Yval ~ ((TDR-Xlow)/(Xhigh-Xlow))^q,start=startlist,data=dat.temp,algorithm="port",
                     lower=c(0,0.01,0.01),upper=c(0.007,0.6,3))
  
  
  # get predicted values and 95% confidence intervals by bootstrapping
  if (type=="g1") newdat[[i]] <- expand.grid(Species=levels(dat.temp$Species), TDR=seq(from=coef(fit.sp[[i]])["Xlow"]+0.01,to=max(dat.temp$TDR),length.out=99),lower=NA,upper=NA)
  if (type=="NSL") newdat[[i]] <- expand.grid(Species=levels(dat.temp$Species), TDR=seq(from=coef(fit.sp[[i]])["Xlow"]+0.01,to=coef(fit.sp[[i]])["Xhigh"],length.out=99),lower=NA,upper=NA)
  newdat[[i]]$wpred <- predict(fit.sp[[i]],newdat[[i]],level=0,se.fit=T)
  
  
  rm(b)
  b <- bootCase(fit.sp[[i]],B=300)
  for(j in 1:nrow(newdat[[i]])){
    rm(b02)
    b02 <- ((newdat[[i]]$TDR[j]-b[,"Xlow"])/(b[,"Xhigh"]-b[,"Xlow"]))^b[,"q"]
    
    newdat[[i]]$lower[j] <- unname(quantile(b02,na.rm=T,probs=c(0.025,0.975)))[1]
    newdat[[i]]$upper[j] <- unname(quantile(b02,na.rm=T,probs=c(0.025,0.975)))[2]
    
  }
  
}

newdatNSL <- newdat
rm(newdat)
fit.spNSL <- fit.sp
rm(fit.sp)
NSL_TDR_beta <- list(fit.spNSL,newdatNSL)
#---------------------------------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------- 

