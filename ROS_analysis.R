#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
#- This is the central analysis script for a leaf-level gas exchange analysis of potted trees experiencing
#   controlled droughts under large outdoor rainout shelters (ROS). Most of the actual data manipulation occurs
#   in functions called by this script.
#  
#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

#- load the libraries and scripts that do all of the actual work.
source("R/loadLibraries.R")


#- plot soil moisture over time
plotVWC(ptsize=1.8,output=T,type="4panel") # also can do type="1panel"

  
#- make a big 16-panel plot of photosynthetic variables over time
plotGX(output=T)
  

#- plot the leaf water potential data over time
plotLWP(fillcol="lightgrey",size=1.75,output=T,labsize=1.8)


#- plot the dependence of g1 and non-stomatal limitation relative to VWC and PSIpd

#- fit g1
g1pars <- returng1()
  

#-- plot g1
dat.l <- split(g1pars,g1pars$Species)
#   
#   require(nlme)
#   fit0 <- nlme(g1~SSasymp(TDR, Asym, R0, lrc), # fit null model with nothing but a random pot effect on the asymtote
#                fixed=list(Asym ~1, lrc~Species, R0 ~ 1),
#                random=R0~1|Date,
#                start=list(fixed=c(Asym=rep(4,1),R0=rep(0.5,1),lrc=rep(2.5,4))),
#                data=g1pars2)
#   


windows(16,16)
par(mfrow=c(4,2),mar=c(0,0.25,0,0.25),xpd=FALSE,oma=c(4,5,1,3),cex=1.6,cex.axis=0.9,cex.lab=0.9)
labs <- c("Cacu","Eusi","Eute","Pira")
for (i in 1:length(dat.l)){
  dat.temp <- dat.l[[i]]
  dat.temp$Species <- factor(dat.temp$Species)
  
  #------------------------------------------------------------------------
  #-- plot g1 vs. TDR
  
  startlist <- list(Xlow = 0.0, Xhigh=0.5, q = 0.6)
  dat.temp$g1norm <- dat.temp$g1/max(dat.temp$g1)
  
  fit.sp[[i]] <- nls(g1norm ~ ((TDR-Xlow)/(Xhigh-Xlow))^q,start=startlist,data=dat.temp,algorithm="port",
                        lower=c(0,0.01,0.01),upper=c(0.007,0.6,3))
  
  #Species[i] <- as.character(fit.sp$Species[1])
  #fit.sp <- nls(g1~SSasymp(TDR,Asym,R0,lrc),data=dat.temp,start=list(Asym=4,R0=0.5,lrc=2.5))
  
  plot(g1norm~TDR,data=subset(dat.temp,Treat=="wet"),pch=21,col="black",bg=grey(0.1),axes=F,ylim=c(0,1.05),xlim=c(0,0.4))
  points(g1norm~TDR,data=subset(dat.temp,Treat=="dry"),pch=21,col="black",bg=grey(0.8))
  magaxis(side=c(1:4),labels=c(0,1,0,0),las=1)
  if(i==4)  magaxis(side=c(1:4),labels=c(1,1,0,0),las=1)
  mtext(labs[i],side=2,xpd=T,cex=1.3,line=1.75)
  if(i==4)  mtext(expression(VWC~(m^3~m^-3)),side=1,outer=F,cex=1.5,line=2)
  
  
  # plot model and SE from bootstrapping
  newdat <- expand.grid(Species=levels(dat.temp$Species), TDR=seq(from=0.01,to=max(dat.temp$TDR),length.out=99),lower=NA,upper=NA)
  newdat$wpred <- predict(fit.sp[[i]],newdat,level=0,se.fit=T)
  
  
  b <- bootCase(fit.sp[[i]],B=300)
  for(j in 1:nrow(newdat)){
    #b02 <- SSasymp(newdat$TDR[j],Xlow=b[,"Xlow"],Xhigh=b[,"Xhigh"],q=b[,"q"])
    b02 <- ((newdat$TDR[j]-b[,"Xlow"])/(b[,"Xhigh"]-b[,"Xlow"]))^b[,"q"]
    
    newdat$lower[j] <- unname(quantile(b02,probs=c(0.025,0.975)))[1]
    newdat$upper[j] <- unname(quantile(b02,probs=c(0.025,0.975)))[2]
    
  }
  lines(wpred~TDR,data=newdat)
  polygon(x = c(newdat$TDR, rev(newdat$TDR)), y = c(newdat$lower, rev(newdat$upper)), col = alpha("grey",0.5), border = NA)
  #------------------------------------------------------------------------
  
  
  # 
  # #------------------------------------------------------------------------
  # #-- plot g1 vs. LWP
  # if (i<4)fit.sp2 <- nls(g1~SSlogis(LWP.pd,Asym,xmid,scale),data=dat.temp,start=list(Asym=4,xmid=-2,scale=1))
  # 
  # 
  # plot(g1~LWP.pd,data=subset(dat.temp,Treat=="wet"),pch=21,col="black",bg=grey(0.1),axes=F,ylim=c(0,5.5),xlim=c(-10,0))
  # points(g1~LWP.pd,data=subset(dat.temp,Treat=="dry"),pch=21,col="black",bg=grey(0.8))
  # magaxis(side=c(1:4),labels=c(0,0,0,1),las=1)
  # if(i==4)  magaxis(side=c(1:4),labels=c(1,0,0,1),las=1)
  # #mtext(labs[i],side=2,xpd=T,cex=1.3,line=1.75)
  # if(i==4)  mtext(expression(LWP[pd]~(MPa)),side=1,outer=F,cex=1.5,line=2)
  # 
  # if (i<4){
  #   # Fixed effects predictions
  #   # Make dataframe with all combinations of Species and conductance.
  #   newdat <- expand.grid(Species=levels(dat.temp$Species), LWP.pd=seq(from=-10,to=0,length.out=99,lower=NA,upper=NA))
  #   newdat$wpred <- predict(fit.sp2,newdat,level=0,se.fit=T)
  #   
  #   
  #   b <- bootCase(fit.sp2,B=200)
  #   for(j in 1:nrow(newdat)){
  #     b02 <- SSlogis(newdat$LWP.pd[j],Asym=b[,"Asym"],xmid=b[,"xmid"],scal=b[,"scale"])
  #     newdat$lower[j] <- unname(quantile(b02,probs=c(0.025,0.975)))[1]
  #     newdat$upper[j] <- unname(quantile(b02,probs=c(0.025,0.975)))[2]
  #     
  #   }
  #   lines(wpred~LWP.pd,data=newdat)
  #   polygon(x = c(newdat$LWP.pd, rev(newdat$LWP.pd)), y = c(newdat$lower, rev(newdat$upper)), col = alpha("grey",0.5), border = NA)
  # }
  #------------------------------------------------------------------------
  
  #if(i<4)lines(wpred~LWP.pd,data=newdat)
  
  # 
  #   
  #   plot(g1~LWP.md,data=subset(dat.temp,Treat=="wet"),pch=21,col="black",bg=grey(0.1),axes=F,ylim=c(0,7),xlim=c(-10,0))
  #   points(g1~LWP.md,data=subset(dat.temp,Treat=="dry"),pch=21,col="black",bg=grey(0.8))
  #   magaxis(side=c(1:4),labels=c(0,0,0,4),las=1)
  #   if(i==4)  magaxis(side=c(1:4),labels=c(1,0,0,1),las=1)
  #   if(i==4)  mtext(expression(psi[l-MD]~(MPa)),side=1,outer=F,cex=1.5,line=2)
  
}
mtext(expression(g[1]),side=2,outer=T,cex=3,las=1,line=2.5)
legend(x=-30,y=-0.3,xpd=NA,legend=c("Wet","Dry"),pch=21,pt.bg=c("black","grey"))
