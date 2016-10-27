#---------------------------------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------- 
# Fits beta models for g1 and non-stomatal limitation, using maximum transpirable water as the predictor.
#---------------------------------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------- 


#- this needs a bunch of work. The new fit_beta model doesn't effectively return "sharp" NSL curves, and 
#  some of the other fits look suspect too.


#- get the g1 fits
g1values <- returng1()



#- get the NSL fits
NSL <- returnVcmaxa()
NSLpars <- summaryBy(NSL+Photo+Cond+Ci+TDR~Species+gxDate+Treat,data=NSL,FUN=mean,keep.names=T)

#---------------------------------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------- 
# fit g1 first

dat <- g1values
dat$TSW <- (dat$TDR)/(max(dat$TDR)-min(dat$TDR)) # normalize TDR data to estimate the transpirable soil water
type="g1"
startlist = list(Xlow = 0.0, Xhigh=0.5, q = 0.4)
#- split into list of species
dat.l <- split(dat,dat$Species)

#- preallocate empty lists to catch output
fit.sp <- newdat <-  list()

#--------------------------------------------------
#-- NLS function that returns "1.0" above Xhigh
fit_beta <- function(Xval,Xlow,Xhigh,q){

  #- predict the Yval given the Xval
  Yval <- ((Xval-Xlow)/(Xhigh-Xlow))^q
  
  #- overwrite the values beyond Xhigh
  Yval[which(Xval>Xhigh)] <- 1
  
  return(Yval)
}
#--------------------------------------------------


#- loop over each element of the list, fit data
for (i in 1:length(dat.l)){
  dat.temp <- dat.l[[i]]
  dat.temp$Species <- factor(dat.temp$Species)
  dat.temp$Xval <- dat.temp$TSW
  #- fit the model
  if (type=="g1") dat.temp$Yval <- dat.temp$g1/max(dat.temp$g1)
  if (type=="NSL") dat.temp$Yval <- dat.temp$NSL
  fit.sp[[i]] <- nls(Yval ~ ((Xval-Xlow)/(Xhigh-Xlow))^q,start=startlist,data=dat.temp,algorithm="port",
                     lower=c(0,0,0.01),upper=c(0.007,1,3))
  #fit.sp[[i]] <- nls(Yval ~ fit_beta(Xval,Xlow,Xhigh,q),start=startlist,data=dat.temp,algorithm="port",
  #                   lower=c(0,0,0.01),upper=c(0.005,1,3))
  
  
  # get predicted values and 95% confidence intervals by bootstrapping
  if (type=="g1") newdat[[i]] <- expand.grid(Species=levels(dat.temp$Species), Xval=seq(from=0.01,to=1,length.out=99),lower=NA,upper=NA)
  if (type=="NSL") newdat[[i]] <- expand.grid(Species=levels(dat.temp$Species), TDR=seq(from=coef(fit.sp[[i]])["Xlow"]+0.01,to=coef(fit.sp[[i]])["Xhigh"],length.out=99),lower=NA,upper=NA)
  newdat[[i]]$wpred <- predict(fit.sp[[i]],newdat[[i]],level=0,se.fit=T)
  
  
  rm(b)
  b <- bootCase(fit.sp[[i]],B=300)
  for(j in 1:nrow(newdat[[i]])){
    rm(b02)
    b02 <- ((newdat[[i]]$Xval[j]-b[,"Xlow"])/(b[,"Xhigh"]-b[,"Xlow"]))^b[,"q"]
    #b02 <- fit_beta(Xval=newdat[[i]]$Xval[j],Xlow=b[,"Xlow"],Xhigh=b[,"Xhigh"],q=b[,"q"])
    newdat[[i]]$lower[j] <- unname(quantile(b02,probs=c(0.025,0.975)))[1]
    newdat[[i]]$upper[j] <- unname(quantile(b02,probs=c(0.025,0.975)))[2]
    
  }
  
}
newdatg1.TSW <- newdat
fit.spg1.TSW <- fit.sp
rm(newdat)
rm(fit.sp)
g1_TDR_beta.TSW <- list(fit.spg1.TSW,newdatg1.TSW)

#---------------------------------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------- 




#---------------------------------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------- 
# fit NSL next 

dat <- NSLpars
dat$TSW <- (dat$TDR)/(max(dat$TDR)-min(dat$TDR)) # normalize TDR data to estimate the transpirable soil water

type="NSL"
startlist= list(Xlow = 0.001, Xhigh=0.4, q = 1)
#- split into list of species
dat.l <- split(dat,dat$Species)

#- preallocate empty lists to catch output
fit.sp <- newdat <-  list()

#- loop over each element of the list, fit data
for (i in 1:length(dat.l)){
  dat.temp <- subset(dat.l[[i]],TSW<0.8)
  #dat.temp <- dat.l[[i]]
  dat.temp$Species <- factor(dat.temp$Species)
  dat.temp$Xval <- dat.temp$TSW
  
  #- fit the model
  if (type=="g1") dat.temp$Yval <- dat.temp$g1/max(dat.temp$g1)
  if (type=="NSL") dat.temp$Yval <- dat.temp$NSL
  fit.sp[[i]] <- nls(Yval ~ ((Xval-Xlow)/(Xhigh-Xlow))^q,start=startlist,data=dat.temp,algorithm="port",trace=F,
                    lower=c(0,0.01,0.01),upper=c(0.01,0.9,3))
  #fit.sp[[i]] <- nls(Yval ~ fit_beta(Xval,Xlow,Xhigh,q),start=startlist,data=dat.temp,algorithm="port",
  #                   lower=c(0,0,0.01),upper=c(0.006,1,3))
  
  
  # get predicted values and 95% confidence intervals by bootstrapping
  if (type=="g1") newdat[[i]] <- expand.grid(Species=levels(dat.temp$Species), TDR=seq(from=coef(fit.sp[[i]])["Xlow"]+0.01,to=max(dat.temp$TDR),length.out=99),lower=NA,upper=NA)
  if (type=="NSL") newdat[[i]] <- expand.grid(Species=levels(dat.temp$Species), Xval=seq(from=0.01,to=1,length.out=99),lower=NA,upper=NA)
  newdat[[i]]$wpred <- predict(fit.sp[[i]],newdat[[i]],level=0,se.fit=T)
  
  
  rm(b)
  b <- bootCase(fit.sp[[i]],B=300)
  for(j in 1:nrow(newdat[[i]])){
    rm(b02)
    #b02 <- ((newdat[[i]]$TDR[j]-b[,"Xlow"])/(b[,"Xhigh"]-b[,"Xlow"]))^b[,"q"]
    b02 <- fit_beta(Xval=newdat[[i]]$Xval[j],Xlow=b[,"Xlow"],Xhigh=b[,"Xhigh"],q=b[,"q"])
    
    newdat[[i]]$lower[j] <- unname(quantile(b02,na.rm=T,probs=c(0.025,0.975)))[1]
    newdat[[i]]$upper[j] <- unname(quantile(b02,na.rm=T,probs=c(0.025,0.975)))[2]
    
  }
  
}

newdatNSL.TSW <- newdat
rm(newdat)
fit.spNSL.TSW <- fit.sp
rm(fit.sp)
NSL_TDR_beta.TSW <- list(fit.spNSL.TSW,newdatNSL.TSW)
#---------------------------------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------- 

