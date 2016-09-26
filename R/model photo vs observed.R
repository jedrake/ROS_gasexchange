#------------------------------------------------------------------------------------------------------------------
#- Script to make a big figure of predicted vs. observed photosynthetic rates for the six model frameworks
#------------------------------------------------------------------------------------------------------------------




#------------------------------------------------------------------------------------------------------------------
#- get the gas exchange, VWC, and LWP data
dat.all2 <- return.gx.vwc.lwp()

#- subset to just the variables I want, to make things a little easier
dat.all <- dat.all2[,c("Species","Treat","Pot","Date","Photo","Cond","VpdL","Tleaf","CO2S","PARi","LWP.pd","LWP.md","TDR")]
dat.all <- dat.all[complete.cases(dat.all),]
#------------------------------------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------------------------------------
#- read in the parameter values from table1.csv
params1 <- read.csv("Output/table1.csv")
params1$Species <- factor(rep(c("cacu","eusi","eute","pira"),4))
params1$Yvar <- factor(c(rep("s",8),rep("ns",8)))
params1$Xvar <- factor(c(rep("theta",4),rep("lwp",4),rep("theta",4),rep("lwp",4)))
params1$Xo <- as.numeric(substr(params1$Xo,start=1,stop=4))
params1$Xh <- as.numeric(substr(params1$Xh,start=1,stop=4))
params1$q <- as.numeric(substr(params1$q,start=1,stop=4))

#- get the low parameter
params_Xo <- reshape2::dcast(params1,Species~Xvar+Yvar,value.var="Xo")
names(params_Xo)[2:5] <- paste("Xl",names(params_Xo)[2:5],sep="_")

#- get the high parameter
params_Xh <- reshape2::dcast(params1,Species~Xvar+Yvar,value.var="Xh")
names(params_Xh)[2:5] <- paste("Xh",names(params_Xh)[2:5],sep="_")

#- get the q parameter
params_q <- reshape2::dcast(params1,Species~Xvar+Yvar,value.var="q")
names(params_q)[2:5] <- paste("q",names(params_q)[2:5],sep="_")

params2 <- merge(params_Xo,params_Xh,by="Species")
params3 <- merge(params2,params_q,by="Species")

#- manually setup the Tuzet parameters
params.tz <- data.frame(Species=factor(c("cacu","eusi","eute","pira")),
                     #- Tuzet model
                     psiv = c(-1.46,-.02,-.35,-1.16),Sf = c(1.46,.76,.61,3.41),K = c(5.74,2.96,3.5,4.06))

params <- merge(params3,params.tz,by="Species")
#------------------------------------------------------------------------------------------------------------------

# 
# #- parameter estimates for the B and Tuzets models
# params <- data.frame(Species=factor(c("cacu","eusi","eute","pira")),
#                      #- Beta functions, soil water content
#                      Xl_theta_s = c(.01,.01,0,0),Xh_theta_s = c(.33,.42,.6,.37),q_theta_s = c(0.37,.44,.32,.83),
#                      Xl_theta_ns = c(.01,.01,.01,.01),xh_theta_ns = c(.19,.2,.18,.12),q_theta_ns=c(.55,.38,.34,.48),
#                      #- beta functions, leaf water potential
#                      Xl_lwp_s = c(-10,-10,-10,-6),Xh_lwp_s = c(0,0,0,0),q_lwp_s = c(2.49,4.56,3.33,6),
#                      Xl_lwp_ns = c(-9,-10,-9,-2.5),Xh_lwp_ns = c(-0.16,0,0,-.27),q_lwp_ns=c(3.65,2.38,1.84,1.1),
#                      #- Tuzet model
#                      psiv = c(-1.46,-.02,-.35,-1.16),Sf = c(1.46,.76,.61,3.41),K = c(5.74,2.96,3.5,4.06))

dat.all3 <- merge(dat.all,params,by="Species")

#- calculate the beta terms
dat.all3$Beta_theta_s <- with(dat.all3, ((TDR-Xl_theta_s)/(Xh_theta_s-Xl_theta_s))^q_theta_s)
dat.all3$Beta_theta_ns <- with(dat.all3, ((TDR-Xl_theta_ns)/(Xh_theta_ns-Xl_theta_ns))^q_theta_ns)
dat.all3$Beta_lwp_s <- with(dat.all3, ((LWP.pd-Xl_lwp_s)/(Xh_lwp_s-Xl_lwp_s))^q_lwp_s)
dat.all3$Beta_lwp_ns <- with(dat.all3, ((LWP.pd-Xl_lwp_ns)/(Xh_lwp_ns-Xl_lwp_ns))^q_lwp_ns)

#- convert Beta values >1 to 1
dat.all3$Beta_theta_s[which(dat.all3$Beta_theta_s >1)] <- 1
dat.all3$Beta_theta_ns[which(dat.all3$Beta_theta_ns >1)] <- 1
dat.all3$Beta_lwp_s[which(dat.all3$Beta_lwp_s >1)] <- 1
dat.all3$Beta_lwp_ns[which(dat.all3$Beta_lwp_ns >1)] <- 1



#- convert values lower than Xl with zero
dat.all3$Beta_theta_s[which(dat.all3$TDR < dat.all3$Xl_theta_s)] <- 0
dat.all3$Beta_theta_ns[which(dat.all3$TDR < dat.all3$Xl_theta_ns)] <- 0
dat.all3$Beta_lwp_s[which(dat.all3$LWP < dat.all3$Xl_lwp_s )] <- 0
dat.all3$Beta_lwp_ns[which(dat.all3$LWP < dat.all3$Xl_lwp_ns )] <- 0

dat.all3 <- dat.all3[complete.cases(dat.all3),] # remove 10 points with missing data
#- convert a very few NA's to zero
#dat.all3$Beta_lwp_ns[which(is.na(dat.all3$Beta_lwp_ns))]  <- 0
#------------------------------------------------------------------------------------------------------------------







#------------------------------------------------------------------------------------------------------------------
#-- predict photosynthetic rates under each modeling condition
#------------------------------------------------------------------------------------------------------------------

#- subset to dry data only?
dat.all3 <- subset(dat.all3,Treat=="dry")

#- Null (no change in physiology with drought)
dat.all3$P_null <- Photosyn(VPD=dat.all3$VpdL,Ca=dat.all3$CO2S,PPFD=dat.all3$PARi,Tleaf=dat.all3$Tleaf,
                            g1=5,Vcmax=80,Jmax=1.6*80)$ALEAF

#- Beta model based on theta
dat.all3[,c("P_theta_s","Gs_theta_s")] <- Photosyn(VPD=dat.all3$VpdL,Ca=dat.all3$CO2S,PPFD=dat.all3$PARi,Tleaf=dat.all3$Tleaf,
                              g1=5*dat.all3$Beta_theta_s,Vcmax=80,Jmax=1.6*80)[,c("ALEAF","GS")]
dat.all3[,c("P_theta_ns","Gs_theta_ns")] <- Photosyn(VPD=dat.all3$VpdL,Ca=dat.all3$CO2S,PPFD=dat.all3$PARi,Tleaf=dat.all3$Tleaf,
                              g1=5,Vcmax=(79*dat.all3$Beta_theta_ns+1),Jmax=(1.6*79*dat.all3$Beta_theta_ns+1))[,c("ALEAF","GS")]
dat.all3[,c("P_theta_sns","Gs_theta_sns")] <- Photosyn(VPD=dat.all3$VpdL,Ca=dat.all3$CO2S,PPFD=dat.all3$PARi,Tleaf=dat.all3$Tleaf,
                               g1=5*dat.all3$Beta_theta_s,Vcmax=(79*dat.all3$Beta_theta_ns+1),Jmax=(1.6*79*dat.all3$Beta_theta_ns+1))[,c("ALEAF","GS")]

#- Beta models based on lwp
dat.all3[,c("P_lwp_s","Gs_lwp_s")]<- Photosyn(VPD=dat.all3$VpdL,Ca=dat.all3$CO2S,PPFD=dat.all3$PARi,Tleaf=dat.all3$Tleaf,
                              g1=5*dat.all3$Beta_lwp_s,Vcmax=80,Jmax=1.6*80)[,c("ALEAF","GS")]
dat.all3[,c("P_lwp_ns","Gs_lwp_ns")] <- Photosyn(VPD=dat.all3$VpdL,Ca=dat.all3$CO2S,PPFD=dat.all3$PARi,Tleaf=dat.all3$Tleaf,
                               g1=5,Vcmax=(79*dat.all3$Beta_lwp_ns+1),Jmax=(1.6*79*dat.all3$Beta_lwp_ns+1))[,c("ALEAF","GS")]
dat.all3[,c("P_lwp_sns","Gs_lwp_sns")] <- Photosyn(VPD=dat.all3$VpdL,Ca=dat.all3$CO2S,PPFD=dat.all3$PARi,Tleaf=dat.all3$Tleaf,
                           g1=5*dat.all3$Beta_lwp_s,Vcmax=(79*dat.all3$Beta_lwp_ns+1),Jmax=(1.6*79*dat.all3$Beta_lwp_ns+1))[,c("ALEAF","GS")]


#- Tuzet alone
dat.all3[,c("P_tuzet_s","Gs_tuzet_s")] <- photosyn(SF=dat.all3$Sf, PSIV=dat.all3$psiv, G0=0.005, VCMAX=80,G1=5,K=dat.all3$K,
                CS=dat.all3$CO2S,WEIGHTEDSWP=dat.all3$LWP.pd, VPD=dat.all3$VpdL,PAR=dat.all3$PARi,TLEAF=dat.all3$Tleaf)[,c("ALEAF","GS")]
dat.all3[,c("P_tuzet_sns","Gs_tuzet_sns")] <- photosyn(SF=dat.all3$Sf, PSIV=dat.all3$psiv, G0=0.005,G1=5,K=dat.all3$K,
                                 VCMAX=(79*dat.all3$Beta_lwp_ns+1),JMAX=(1.6*79*dat.all3$Beta_lwp_ns+1),
                               CS=dat.all3$CO2S,WEIGHTEDSWP=dat.all3$LWP.pd, VPD=dat.all3$VpdL,PAR=dat.all3$PARi,TLEAF=dat.all3$Tleaf)[,c("ALEAF","GS")]





#--------------------------------
#- plot observed vs. predicted
windows(25,35)
par(mfrow=c(5,2),oma=c(6,9,4,3),mar=c(0,0,0,0),xpd=F)
colors <- brewer.pal(4,"Set1")
fitcol=alpha("darkgrey",0.8)
textsize <- 1.7
xlims <- c(-2,23)
ylims <- c(-2,35)

# #- plot null model (twice)
# plot(dat.all3$Photo~dat.all3$P_null,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
# magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,tcl=0.3)
# title(ylab="Null",xpd=NA,cex.lab=textsize)
# title(xlab=expression(theta~"as"~"predictor"),cex.lab=textsize,xpd=NA,line=-10)
# legend("bottomright",letters[1],bty="n",cex=1.2)
# 
# plot(dat.all3$Photo~dat.all3$P_null,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
# magaxis(side=c(1,2,4),labels=c(0,0,1),frame.plot=T,las=1,tcl=0.3)
# title(xlab=expression(psi[l]~"as"~"predictor"),cex.lab=textsize,xpd=NA,line=-10)
# legend("bottomright",letters[5],bty="n",cex=1.2)

#- plot B models for stomatal limitation (theta and lwp)
par(xpd=F)
plot(dat.all3$Photo~dat.all3$P_theta_s,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="",xpd=F);abline(0,1)
predline(lm(Photo~P_theta_s,data=dat.all3),fittype="confidence",col=fitcol,xpd=F)
magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,tcl=0.3)
title(ylab=expression(beta[s]),xpd=NA,cex.lab=textsize)
legend("bottomright",letters[2],bty="n",cex=1.2)
r2val <-round(summary(lm(Photo~P_theta_s,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")


plot(dat.all3$Photo~dat.all3$P_lwp_s,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(Photo~P_lwp_s,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(0,0,1),frame.plot=T,las=1,tcl=0.3)
legend("bottomright",letters[6],bty="n",cex=1.2)
r2val <-round(summary(lm(Photo~P_lwp_s,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")

#- plot B models for non-stomatal limitation (theta and then lwp)
plot(dat.all3$Photo~dat.all3$P_theta_ns,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(Photo~P_theta_ns,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,tcl=0.3)
title(ylab=expression(beta[ns]),xpd=NA,cex.lab=textsize)
legend("bottomright",letters[3],bty="n",cex=1.2)
r2val <-round(summary(lm(Photo~P_theta_ns,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")

plot(dat.all3$Photo~dat.all3$P_lwp_ns,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(Photo~P_lwp_ns,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(0,0,1),frame.plot=T,las=1,tcl=0.3)
legend("bottomright",letters[7],bty="n",cex=1.2)
r2val <-round(summary(lm(Photo~P_lwp_ns,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")

#- plot B models with both
plot(dat.all3$Photo~dat.all3$P_theta_sns,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(Photo~P_theta_sns,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,tcl=0.3)
title(ylab=expression(beta[s+ns]),xpd=NA,cex.lab=textsize)
legend("bottomright",letters[4],bty="n",cex=1.2)
r2val <-round(summary(lm(Photo~P_theta_sns,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")

plot(dat.all3$Photo~dat.all3$P_lwp_sns,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(Photo~P_lwp_sns,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(0,0,1),frame.plot=T,las=1,tcl=0.3)
legend("bottomright",letters[8],bty="n",cex=1.2)
r2val <-round(summary(lm(Photo~P_lwp_sns,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")

#- plot Tuzet models (stomatal, and both). Note the two empty plots to fill space
plot(Photo~P_theta_sns,data=dat.all3,type="n",axes=F,xlab="",ylab="")
plot(dat.all3$Photo~dat.all3$P_tuzet_s ,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(Photo~P_tuzet_s,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(0,1,1),frame.plot=T,las=1,tcl=0.3)
title(ylab=expression(Tuzet[s]),xpd=NA,cex.lab=textsize)
legend("bottomright",letters[9],bty="n",cex=1.2)
r2val <-round(summary(lm(Photo~P_tuzet_s,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")


plot(Photo~P_theta_sns,data=dat.all3,type="n",axes=F,xlab="",ylab="")
plot(dat.all3$Photo~dat.all3$P_tuzet_sns ,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(Photo~P_tuzet_sns,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(1,1,1),frame.plot=T,las=1,tcl=0.3)
title(ylab=expression(Tuzet[s+ns]),xpd=NA,cex.lab=textsize)
legend("bottomright",letters[10],bty="n",cex=1.2)
r2val <-round(summary(lm(Photo~P_tuzet_sns,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")

title(xlab=expression("Predicted"~A[sat]~(mu*mol~m^-2~s^-1)),outer=T,cex.lab=2)
title(ylab=expression("Observed"~A[sat]~(mu*mol~m^-2~s^-1)),outer=T,cex.lab=2,line=5)

#-- calculate concordance correlation coefficients
null <- summary(agreement(x=dat.all3$Photo,y=dat.all3$P_null,error="constant",TDI_a=20,target="fixed"))
beta_s_theta <- summary(agreement(x=dat.all3$Photo,y=dat.all3$P_theta_s,error="constant",TDI_a=20,target="fixed"))
beta_s_lwp <- summary(agreement(x=dat.all3$Photo,y=dat.all3$P_lwp_s,error="constant",TDI_a=20,target="fixed"))
beta_ns_theta <- summary(agreement(x=dat.all3$Photo,y=dat.all3$P_lwp_ns,error="constant",TDI_a=20,target="fixed"))
beta_ns_lwp <- summary(agreement(x=dat.all3$Photo,y=dat.all3$P_lwp_ns,error="constant",TDI_a=20,target="fixed"))
beta_sns_theta <- summary(agreement(x=dat.all3$Photo,y=dat.all3$P_lwp_sns,error="constant",TDI_a=20,target="fixed"))
beta_sns_lwp <- summary(agreement(x=dat.all3$Photo,y=dat.all3$P_lwp_sns,error="constant",TDI_a=20,target="fixed"))
tuzet_s <- summary(agreement(x=dat.all3$Photo,y=dat.all3$P_tuzet_s,error="constant",TDI_a=20,target="fixed"))
tuzet_sns <- summary(agreement(x=dat.all3$Photo,y=dat.all3$P_tuzet_sns,error="constant",TDI_a=20,target="fixed"))


#-- calculate the residual standard deviation
summary(lm(Photo~P_theta_s,data=dat.all3))$sigma
summary(lm(Photo~P_lwp_s,data=dat.all3))$sigma
summary(lm(Photo~P_theta_ns,data=dat.all3))$sigma
summary(lm(Photo~P_lwp_ns,data=dat.all3))$sigma
summary(lm(Photo~P_theta_sns,data=dat.all3))$sigma
summary(lm(Photo~P_lwp_sns,data=dat.all3))$sigma
summary(lm(Photo~P_theta_sns,data=dat.all3))$sigma
summary(lm(Photo~P_lwp_sns,data=dat.all3))$sigma
summary(lm(Photo~P_tuzet_s,data=dat.all3))$sigma
summary(lm(Photo~P_tuzet_sns,data=dat.all3))$sigma









#--------------------------------
#--------------------------------
#- repeat, but for CONDUCTANCE!
#--------------------------------
#--------------------------------

#--------------------------------
#- plot observed vs. predicted
windows(25,35)
par(mfrow=c(5,2),oma=c(6,9,4,3),mar=c(0,0,0,0),xpd=F)
colors <- brewer.pal(4,"Set1")
fitcol=alpha("darkgrey",0.8)
textsize <- 1.7
xlims <- c(-0.03,0.6)
ylims <- c(-0.03,0.6)

# #- plot null model (twice)
# plot(dat.all3$Photo~dat.all3$Gs_null,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
# magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,tcl=0.3)
# title(ylab="Null",xpd=NA,cex.lab=textsize)
# title(xlab=expression(theta~"as"~"predictor"),cex.lab=textsize,xpd=NA,line=-10)
# legend("bottomright",letters[1],bty="n",cex=1.2)
# 
# plot(dat.all3$Photo~dat.all3$Gs_null,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
# magaxis(side=c(1,2,4),labels=c(0,0,1),frame.plot=T,las=1,tcl=0.3)
# title(xlab=expression(psi[l]~"as"~"predictor"),cex.lab=textsize,xpd=NA,line=-10)
# legend("bottomright",letters[5],bty="n",cex=1.2)

#- plot B models for stomatal limitation (theta and lwp)
par(xpd=F)
plot(dat.all3$Cond~dat.all3$Gs_theta_s,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="",xpd=F);abline(0,1)
predline(lm(Cond~Gs_theta_s,data=dat.all3),fittype="confidence",col=fitcol,xpd=F)
magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,tcl=0.3)
title(ylab=expression(beta[s]),xpd=NA,cex.lab=textsize)
legend("bottomright",letters[2],bty="n",cex=1.2)
r2val <-round(summary(lm(Cond~Gs_theta_s,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")


plot(dat.all3$Cond~dat.all3$Gs_lwp_s,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(Cond~Gs_lwp_s,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(0,0,1),frame.plot=T,las=1,tcl=0.3)
legend("bottomright",letters[6],bty="n",cex=1.2)
r2val <-round(summary(lm(Cond~Gs_lwp_s,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")

#- plot B models for non-stomatal limitation (theta and then lwp)
plot(dat.all3$Cond~dat.all3$Gs_theta_ns,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(Cond~Gs_theta_ns,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,tcl=0.3)
title(ylab=expression(beta[ns]),xpd=NA,cex.lab=textsize)
legend("bottomright",letters[3],bty="n",cex=1.2)
r2val <-round(summary(lm(Cond~Gs_theta_ns,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")

plot(dat.all3$Cond~dat.all3$Gs_lwp_ns,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(Cond~Gs_lwp_ns,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(0,0,1),frame.plot=T,las=1,tcl=0.3)
legend("bottomright",letters[7],bty="n",cex=1.2)
r2val <-round(summary(lm(Cond~Gs_lwp_ns,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")

#- plot B models with both
plot(dat.all3$Cond~dat.all3$Gs_theta_sns,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(Cond~Gs_theta_sns,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,tcl=0.3)
title(ylab=expression(beta[s+ns]),xpd=NA,cex.lab=textsize)
legend("bottomright",letters[4],bty="n",cex=1.2)
r2val <-round(summary(lm(Cond~Gs_theta_sns,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")

plot(dat.all3$Cond~dat.all3$Gs_lwp_sns,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(Cond~Gs_lwp_sns,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(0,0,1),frame.plot=T,las=1,tcl=0.3)
legend("bottomright",letters[8],bty="n",cex=1.2)
r2val <-round(summary(lm(Cond~Gs_lwp_sns,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")

#- plot Tuzet models (stomatal, and both). Note the two empty plots to fill space
plot(Cond~Gs_theta_sns,data=dat.all3,type="n",axes=F,xlab="",ylab="")
plot(dat.all3$Cond~dat.all3$Gs_tuzet_s ,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(Cond~Gs_tuzet_s,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(0,1,1),frame.plot=T,las=1,tcl=0.3)
title(ylab=expression(Tuzet[s]),xpd=NA,cex.lab=textsize)
legend("bottomright",letters[9],bty="n",cex=1.2)
r2val <-round(summary(lm(Cond~Gs_tuzet_s,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")


plot(Cond~Gs_theta_sns,data=dat.all3,type="n",axes=F,xlab="",ylab="")
plot(dat.all3$Cond~dat.all3$Gs_tuzet_sns ,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(Cond~Gs_tuzet_sns,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(1,1,1),frame.plot=T,las=1,tcl=0.3)
title(ylab=expression(Tuzet[s+ns]),xpd=NA,cex.lab=textsize)
legend("bottomright",letters[10],bty="n",cex=1.2)
r2val <-round(summary(lm(Cond~Gs_tuzet_sns,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")

title(xlab=expression("Predicted"~g[s]~(mol~m^-2~s^-1)),outer=T,cex.lab=2)
title(ylab=expression("Observed"~g[s]~(mol~m^-2~s^-1)),outer=T,cex.lab=2,line=5)
