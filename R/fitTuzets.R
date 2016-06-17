#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
#- This script is an attempt to fit the Tuzets model to the ROS data
#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
source("R/loadLibraries.R")

#- get the soil water curve
soilPars <- plotMoistCurve(output=T)

#------------------------------------------------------------------------------------------------------------------
#- function to predict soil water potentail from VWC measurements
predPsiSoil <- function(VWC,pars){
  thetaSat <- pars[1]
  b <- pars[2]
  PsiE <- pars[3]
  
  PsiSoil <- PsiE*(VWC/thetaSat)^(-1*b)
  return(PsiSoil)
}
#------------------------------------------------------------------------------------------------------------------




#------------------------------------------------------------------------------------------------------------------
#- get the gas exchange, VWC, and LWP data
dat.all2 <- return.gx.vwc.lwp()

#- subset to just the variables I want, to make things a little easier
dat.all <- dat.all2[,c("Species","Treat","Pot","Date","Photo","Cond","VpdL","Tleaf","CO2S","PARi","LWP","LWP.md","TDR")]

#- predict soil water potential
dat.all$SWP <- predPsiSoil(VWC=(dat.all$TDR/100),pars=soilPars)
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
  
  
  #- adjust Vcmax according to a beta function. This uses the pira parameters
  if (adjVcmax==1){
    Vcmax_a <- Vcmax*((dat$TDR/100 - 0.01)/(0.17-0.01))^0.34
    Vcmax_a[which(Vcmax_a>Vcmax)] <- Vcmax
  }
  
  #- or bypass the Vcmax adjustment
  if (adjVcmax==0){
    Vcmax_a <- Vcmax
  }
  
  #- model
  out <- photosyn(SF=SF, PSIV=psiv, G0=0.005, VCMAX=Vcmax_a,G1=g1,K=K,
                 CS=dat$CO2S,WEIGHTEDSWP=dat$LWP, VPD=dat$VpdL,PAR=dat$PARi,TLEAF=dat$Tleaf)
  
  
  
  #out2 <- Photosyn(g0=0.005, Vcmax=Vcmax,g1=g1,
  #                          Ca=dat$CO2S, VPD=dat$VpdL,PPFD=dat$PARi,Tleaf=dat$Tleaf)
  
  #- calculate the data-model mismatch
  #resid.gs <- sqrt((out$GS - dat$Cond)^2)/(abs(out$GS))
  #resid.A <- sqrt((out$ALEAF - dat$Photo)^2)/(abs(out$ALEAF))
  #resid.psi <- sqrt((out$PSIL - dat$LWP.md)^2)/(abs(out$PSIL))
  
  resid.gs <- sqrt((out$GS - dat$Cond)^2)*10
  resid.A <- sqrt((out$ALEAF - dat$Photo)^2)
  resid.psi <- sqrt((out$PSIL - dat$LWP.md)^2)
  
  
  resid.sum <- sum(resid.gs,resid.A,resid.psi)
  
  if (fit==1) return(resid.sum)
  if (fit==0) return(out)
  
}
#------------------------------------------------------------------------------------------------------------------




#------------------------------------------------------------------------------------------------------------------
#- set up model estimate for a selected species
dat <- subset(dat.all,Species=="eusi")
dat <- dat[complete.cases(dat),]

#- define lower and upper limits of parameter estimates
# pars are psiv, SF,  and K
lower <- c(-10,0.5,0)
upper <- c(0,25,20)
NPmax <- 40
maxiter <- 30

#- fit the model, extract the best parameters, and rerun the model to get the predicted values  
DEfit <- DEoptim(fn=tuzets.cost,lower=lower,upper=upper,dat=dat,g1=5,Vcmax=80,adjVcmax=1,
                     fit=1,DEoptim.control(NP = NPmax,itermax=maxiter)) #perhaps adjust reltol, NP=20, itermax = 100
DEfit.best <- unname(DEfit$optim$bestmem)
DEpred <- tuzets.cost(pars=DEfit.best,dat=dat,fit=0,g1=5,Vcmax=80,adjVcmax=1)


#- overlay plot of estimated values and model
windows();par(mfrow=c(2,2))

#- photo vs. cond
plot(Photo~Cond,data=dat)
points(ALEAF~GS,data=DEpred,col="red")
legend("bottomright",pch=1,col=c("black","red"),legend=c("Obs","Pred"))

plot(dat$Photo~DEpred$ALEAF,xlab="Pred Photo",ylab="Photo");abline(0,1)
plot(dat$Cond~DEpred$GS,xlab="Pred cond",ylab="Cond");abline(0,1)
plot(dat$LWP.md~DEpred$PSIL,xlab="Pred Psi leaf",ylab="Psi leaf");abline(0,1)
title(main=dat$Species[1],outer=T,line=-2)

windows();par(mfrow=c(1,1))
plot(Cond~LWP.md,data=dat)
points(DEpred$GS~dat$LWP.md,col="red")
#------------------------------------------------------------------------------------------------------------------
