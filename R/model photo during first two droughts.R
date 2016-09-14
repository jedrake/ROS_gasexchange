#-----------------------------------------------------------------------------------------
#--- adaptation of "fitBeta_g1_nsl.R" script.
#    The goal here is to fit the model to teh final drought data, and then predict
#     what we should have observed for the first two droughts. This may be a way to examine
#     drought acclimation.
#-----------------------------------------------------------------------------------------



#- get the g1 fits
g1values <- returng1()

#- get the NSL fits
NSL <- returnVcmaxa()
NSLpars <- summaryBy(NSL+Photo+Cond+Ci+TDR~Species+gxDate+Treat,data=NSL,FUN=mean,keep.names=T)

#---------------------------------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------- 
# fit g1 first, last drought data ONLY

dat <- subset(g1values,Date > as.Date("2013-02-01"))
type="g1"
startlist = list(Xlow = 0.0, Xhigh=0.5, q = 0.6)
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
  if (i ==4) startlist = list(Xlow = 0, Xhigh=0.5, q = 2)
  
  fit.sp[[i]] <- nls(Yval ~ ((TDR-Xlow)/(Xhigh-Xlow))^q,start=startlist,data=dat.temp,algorithm="port",trace=F,
                     lower=c(0,0.01,0.01),upper=c(0.007,0.6,3),nls.control(maxiter=100))
  
  
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

dat <- subset(NSLpars,gxDate>as.Date("2013-02-01"))
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





#---------------------------------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------- 
#- fit the data from the prior few droughts



#---------------------------------------------------------------------------------------------------------- 
#---------------------------------------------------------------------------------------------------------- 
# get the data
ros2 <- return.gx.vwc()
dat <- subset(ros2,gxDate<as.Date("2013-02-01"))


#- subset to just the variables I want, to make things a little easier
dat.all <- dat[,c("Species","Treat","Pot","Date","gxDate","Photo","Cond","VpdL","Tleaf","CO2S","PARi","TDR")]
dat.all <- dat.all[complete.cases(dat.all),]
#------------------------------------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------------------------------------
#- parameter estimates for the B and Tuzets models
params <- data.frame(Species=factor(c("cacu","eusi","eute","pira")),
                     #- Beta functions, soil water content
                     Xl_theta_s = c(.07,.07,0,.0004),Xh_theta_s = c(.307,.398,.6,.321),q_theta_s = c(0.47,.56,.31,.93),
                     Xl_theta_ns = c(.007,.007,.007,.007),xh_theta_ns = c(.197,.145,.133,.119),q_theta_ns=c(.909,.78,.678,.632))
dat.all3 <- merge(dat.all,params,by="Species")

#- calculate the beta terms
dat.all3$Beta_theta_s <- with(dat.all3, ((TDR/100-Xl_theta_s)/(Xh_theta_s-Xl_theta_s))^q_theta_s)
dat.all3$Beta_theta_ns <- with(dat.all3, ((TDR/100-Xl_theta_ns)/(xh_theta_ns-Xl_theta_ns))^q_theta_ns)

#- convert Beta values >1 to 1
dat.all3$Beta_theta_s[which(dat.all3$Beta_theta_s >1)] <- 1
dat.all3$Beta_theta_ns[which(dat.all3$Beta_theta_ns >1)] <- 1

#- convert values lower than Xl with zero
dat.all3$Beta_theta_s[which(dat.all3$TDR/100 < dat.all3$Xl_theta_s)] <- 0
dat.all3$Beta_theta_ns[which(dat.all3$TDR/100 < dat.all3$Xl_theta_ns)] <- 0

#------------------------------------------------------------------------------------------------------------------




#------------------------------------------------------------------------------------------------------------------
#-- predict photosynthetic rates under each modeling condition
#------------------------------------------------------------------------------------------------------------------

#- Null (no change in physiology with drought)
dat.all3$P_null <- Photosyn(VPD=dat.all3$VpdL,Ca=dat.all3$CO2S,PPFD=dat.all3$PARi,Tleaf=dat.all3$Tleaf,
                            g1=5,Vcmax=80,Jmax=1.6*80)$ALEAF

#- Beta model based on theta
dat.all3$P_theta_s <- Photosyn(VPD=dat.all3$VpdL,Ca=dat.all3$CO2S,PPFD=dat.all3$PARi,Tleaf=dat.all3$Tleaf,
                               g1=5*dat.all3$Beta_theta_s,Vcmax=80,Jmax=1.6*80)$ALEAF
dat.all3$P_theta_ns <- Photosyn(VPD=dat.all3$VpdL,Ca=dat.all3$CO2S,PPFD=dat.all3$PARi,Tleaf=dat.all3$Tleaf,
                                g1=5,Vcmax=80*dat.all3$Beta_theta_ns,Jmax=1.6*80*dat.all3$Beta_theta_ns)$ALEAF
dat.all3$P_theta_sns <- Photosyn(VPD=dat.all3$VpdL,Ca=dat.all3$CO2S,PPFD=dat.all3$PARi,Tleaf=dat.all3$Tleaf,
                                 g1=5*dat.all3$Beta_theta_s,Vcmax=80*dat.all3$Beta_theta_ns,Jmax=1.6*80*dat.all3$Beta_theta_ns)$ALEAF



#------------------------------------------------------------------------------------------------------------------
#- average up the data
dat.m <- summaryBy(Photo+P_null+P_theta_s+P_theta_ns+P_theta_sns~Species+Treat+gxDate,data=dat.all3,
                   FUN=c(mean,standard.error),keep.names=F)


#------------------------------------------------------------------------------------------------------------------
windows(20,20)
par(mfrow=c(1,1),oma=c(1,1,4,0),mar=c(5,6,1,1))

#- plot observed and expected photosynthetic rates during the first two droughts
plot(Photo.mean~P_theta_ns.mean,type="n",ylim=c(0,25),xlim=c(0,25),data=dat.m,axes=F,xlab="",ylab="")

adderrorbars(x=dat.m$P_theta_ns.mean,y=dat.m$Photo.mean,SE=dat.m$Photo.standard.error,direction='updown')
adderrorbars(x=dat.m$P_theta_ns.mean,y=dat.m$Photo.mean,SE=dat.m$P_theta_ns.standard.error,direction='leftright')

plotBy(Photo.mean~P_theta_ns.mean|Species,data=subset(dat.m,gxDate<as.Date("2012-11-15")),legend=F,xlim=c(0,25),ylim=c(0,25),
       legendwhere="bottomright",pch=15,add=T,cex=1.3)

plotBy(Photo.mean~P_theta_ns.mean|Species,data=subset(dat.m,gxDate>as.Date("2012-11-15")),legend=F,add=T,pch=16,cex=1.3)
abline(0,1)
legend(x=6,y=30,legend=c("Cacu1","Cacu2","Eusi1","Eusi2","Eute1","Eute2","Pira1","Pira2"),
       pch=c(15,16),col=c(palette()[c(1,1,2,2,3,3,4,4)]),
       ncol=4,xpd=NA)
magaxis(side=c(1,2),labels=c(1,1),las=1,cex.lab=1.5,frame.plot=T)
title(ylab=expression(Observed~A[sat]~(mu*mol~m^-2~s^-1)),
      xlab=expression(Predicted~A[sat]~(mu*mol~m^-2~s^-1)),cex.lab=2)
#------------------------------------------------------------------------------------------------------------------

