#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
#- This script is an attempt to fit the Tuzets model to the ROS data
#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
source("R/loadLibraries.R")




#------------------------------------------------------------------------------------------------------------------
#- get the gas exchange, VWC, and LWP data
dat.all2 <- return.gx.vwc.lwp()

#- subset to just the variables I want, to make things a little easier
dat.all <- subset(dat.all2[,c("Species","Treat","Pot","Date","Photo","Cond","Trmmol","VpdL","Tleaf","CO2S","PARi","LWP.pd","LWP.md","TDR")],
                  Date > as.POSIXct("2012-11-04")) # the first two dates don't have reliable water potential data (dates merge incorrectly?)

#------------------------------------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------------------------------------
#- wrapper function to return the model-data mismatch for the Tuzets model
#  If fit == 1, the function will return the residual. If fit == 0, the function will return
#   a dataframe of predicted values based on the input parameters and environmental data.
tuzets.cost <- function(pars,dat,fit=1,g1=5,Vcmax=80,adjVcmax=0){
  
  #- pull out the parameters from the pars vector
  psiv <- pars[1]
  SF <- pars[2]
  K <- pars[3]
  
  
  #- adjust Vcmax according to a beta function. 
  if (adjVcmax==1){
    df <- data.frame(species=c("cacu","eusi","eute","pira"),
                     Xl = c(.01,.01,.01,.01),
                     Xh = c(.19,.2,.18,.12),
                     q = c(.55,.38,.34,.48))
    whichline <- which(df$species==dat$Species[1])
    
    Vcmax_a <- Vcmax*((dat$TDR - df[whichline,"Xl"])/(df[whichline,"Xh"]-df[whichline,"Xl"]))^df[whichline,"q"]
    Vcmax_a[which(Vcmax_a>Vcmax)] <- Vcmax
  }
  
  #- or bypass the Vcmax adjustment
  if (adjVcmax==0){
    Vcmax_a <- Vcmax
  }
  
  #- model
  out <- photosyn(SF=SF, PSIV=psiv, G0=0.005, VCMAX=Vcmax_a,G1=g1,K=K,
                 CS=dat$CO2S,WEIGHTEDSWP=dat$LWP.pd, VPD=dat$VpdL,PAR=dat$PARi,TLEAF=dat$Tleaf)
  
  
  
  #- calculate the data-model mismatch. Attempt to weight each equally by dividing by the maximum value
  max.Gs <- max(c(out$GS,dat$Cond))
  max.A <- max(c(out$ALEAF,dat$Photo))
  max.E <- max(c(out$ELEAF,dat$Trmmol))
  max.psi <- max(c(out$PSIL,dat$LWP.md))
  
  resid.gs <- sqrt(((out$GS - dat$Cond)/max.Gs)^2)
  resid.A <- sqrt(((out$ALEAF - dat$Photo)/max.A)^2)
  resid.E <- sqrt(((out$ELEAF - dat$Trmmol)/max.E)^2)
  resid.psi <- sqrt(((out$PSIL - dat$LWP.md)/max.psi)^2)
  
  
  resid.sum <- sum(resid.gs,resid.A,resid.E,resid.psi)
  
  if (fit==1) return(resid.sum)
  if (fit==0) return(out)
  
}
#------------------------------------------------------------------------------------------------------------------




#------------------------------------------------------------------------------------------------------------------
#- set up model estimate for a selected species
dat.all <- dat.all[complete.cases(dat.all),]


#- base model predictions on the mean across dates
dat.m <- summaryBy(.~Species+Treat+Date,FUN=mean,keep.names=T,data=dat.all)

dat.list <- split(dat.m,dat.m$Species)

#- define lower and upper limits of parameter estimates
# pars are psiv, SF,  and K
lower <- c(-10,0.5,0)
upper <- c(0,25,20)
NPmax <- 100
maxiter <- 30


#- set seed for repeatability
set.seed(1234)

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
#- model with no change in photosynthetic capacity
DEfit <- DEfit.best <- DEpred <- list()
for (i in 1:length(dat.list)){
  tofit <- dat.list[[i]]
  
  #- fit the model, extract the best parameters, and rerun the model to get the predicted values  
  DEfit[[i]] <- DEoptim(fn=tuzets.cost,lower=lower,upper=upper,dat=tofit,g1=5,Vcmax=80,adjVcmax=0,
                   fit=1,DEoptim.control(NP = NPmax,itermax=maxiter)) #perhaps adjust reltol, NP=20, itermax = 100
  DEfit.best[[i]] <- unname(DEfit[[i]]$optim$bestmem)
  DEpred[[i]] <- tuzets.cost(pars=DEfit.best[[i]],dat=tofit,fit=0,g1=5,Vcmax=80,adjVcmax=0)
  DEpred[[i]]$Species <- tofit$Species[1]
  DEpred[[i]]$Date <- tofit$Date
  
}






#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
#- this time, model with a change in Vcmax with decreasing TDR

DEfit.nsl <- DEfit.best.nsl <- DEpred.nsl <- list()
for (i in 1:length(dat.list)){
  tofit <- dat.list[[i]]
  
  #- fit the model, extract the best parameters, and rerun the model to get the predicted values  
  DEfit.nsl[[i]] <- DEoptim(fn=tuzets.cost,lower=lower,upper=upper,dat=tofit,g1=5,Vcmax=80,adjVcmax=1,
                        fit=1,DEoptim.control(NP = NPmax,itermax=maxiter)) #perhaps adjust reltol, NP=20, itermax = 100
  DEfit.best.nsl[[i]] <- unname(DEfit.nsl[[i]]$optim$bestmem)
  DEpred.nsl[[i]] <- tuzets.cost(pars=DEfit.best.nsl[[i]],dat=tofit,fit=0,g1=5,Vcmax=80,adjVcmax=1)
  DEpred.nsl[[i]]$Species <- tofit$Species[1]
  DEpred.nsl[[i]]$Date <- tofit$Date
  
}





#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
#- forward prediction of model responses to leaf water potential 

#- params with and without Vcmax change
do.call(rbind,DEfit.best)
# [1,] -0.6571273 0.8781327 5.318193
# [2,] -1.7899833 1.1339467 4.679995
# [3,] -0.5157321 0.5559586 4.596770
# [4,] -1.1086333 8.9338100 5.515464
do.call(rbind,DEfit.best.nsl)
#[1,] -2.671953 24.875104 7.132115
#[2,] -2.574391  1.465136 5.349148
#[3,] -3.033131 24.788989 6.328085
#[4,] -1.282950  2.130503 4.625018

psileaf <- seq(from=-8,to=0,length.out=101)
pred.psileaf <- list()
for(i in 1:length(DEfit.best)){
  
  #- pull out the parameters from the pars vector
  psiv <- DEfit.best[[i]][1]
  SF <- DEfit.best[[i]][2]
  K <- DEfit.best[[i]][3]
  
  pred.psileaf[[i]] <- photosyn(SF=SF, PSIV=psiv, G0=0.0005, VCMAX=80,G1=8,K=K,
                  CS=400,WEIGHTEDSWP=psileaf, VPD=1.2,PAR=1800,TLEAF=25)
  pred.psileaf[[i]]$Species <- levels(dat.all$Species)[i]
}
preddat <- do.call(rbind,pred.psileaf)



#- overlay plots of data for gs vs. psileaf relationship along with model output
dat.m <- summaryBy(Cond+Photo+LWP.pd+LWP.md+TDR~Species+Treat+Date,FUN=c(mean,standard.error),
                   data=subset(dat.all,Treat=="dry"))


windows(25,14)
par(mfrow=c(1,2),oma=c(1,1,0,0),mar=c(6.25,6.25,0.25,0.25))
colors <- brewer.pal(5,"Accent")[c(1,2,3,5)]#brewer.pal(4,"Set1")

#- plot Conductance
plotBy(GS~PSILIN|Species,data=preddat,type="l",legend=F,lwd=3,ylim=c(0,0.6),xlim=c(-8,0),
       col = colors,axes=F,xlab="",ylab="")
adderrorbars(x=dat.m$LWP.md.mean,y=dat.m$Cond.mean,SE=dat.m$Cond.standard.error,direction="updown")
adderrorbars(x=dat.m$LWP.md.mean,y=dat.m$Cond.mean,SE=dat.m$LWP.md.standard.error,direction="leftright")
plotBy(Cond.mean~LWP.md.mean|Species,dat=dat.m,col=colors,pch=16,add=T,legend=F,cex=1.5)
plotBy(GS~PSILIN|Species,data=preddat,type="l",legend=F,lwd=4,ylim=c(0,0.6),xlim=c(-8,0),
       col = colors,axes=F,xlab="",ylab="",add=T)

legend("top",legend=c("Cacu","Eusi","Eute","Pira"),pch=16,col=colors,ncol=2,cex=1.2)
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.2)
title(xlab=expression(Psi[l]~(MPa)),ylab=expression(g[s]~(mol~m^-2~s^-1)),cex.lab=1.5)
legend("topright",letters[1],bty="n",cex=1.1,inset=0.02)

#- plot Photo
plotBy(ALEAF~PSILIN|Species,data=preddat,type="l",legend=F,lwd=3,ylim=c(0,25),xlim=c(-8,0),
       col = colors,axes=F,xlab="",ylab="")
adderrorbars(x=dat.m$LWP.md.mean,y=dat.m$Photo.mean,SE=dat.m$Photo.standard.error,direction="updown")
adderrorbars(x=dat.m$LWP.md.mean,y=dat.m$Photo.mean,SE=dat.m$LWP.md.standard.error,direction="leftright")
plotBy(Photo.mean~LWP.md.mean|Species,dat=dat.m,col=colors,pch=16,add=T,legend=F,cex=1.5)
plotBy(ALEAF~PSILIN|Species,data=preddat,type="l",legend=F,lwd=4,ylim=c(0,25),xlim=c(-8,0),
       col = colors,axes=F,xlab="",ylab="",add=T)
#legend("top",legend=c("Cacu","Eusi","Eute","Pira"),fill=colors,ncol=2,cex=1.2)
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.2)
title(xlab=expression(Psi[l]~(MPa)),ylab=expression(A[sat]~(mu*mol~m^-2~s^-1)),cex.lab=1.5)
legend("topright",letters[2],bty="n",cex=1.1,inset=0.02)

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------






#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
#- write out the parameters for table 2
params1 <- data.frame(Species = c("cacu","eusi","eute","pira"),do.call(rbind,DEfit.best))
names(params1)[2:4] <- c("psiv","SF","K")
write.csv(params1,"Output/table2.csv",row.names=F)
#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
