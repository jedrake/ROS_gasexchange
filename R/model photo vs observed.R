#------------------------------------------------------------------------------------------------------------------
#- Script to make a big figure of predicted vs. observed photosynthetic rates for all of the model frameworks
#------------------------------------------------------------------------------------------------------------------




#------------------------------------------------------------------------------------------------------------------
#- get the gas exchange, VWC, and LWP data
dat.all <- return.gx.vwc.lwp()
dat.all2 <- subset(dat.all,!(gxDate %in% as.Date(c("2012-10-04","2012-10-31"))))

#- subset to just the variables I want, to make things a little easier
dat.all <- dat.all2[,c("Species","Treat","Pot","Date","Photo","Cond","VpdL","Tleaf","CO2S","PARi","LWP.pd","LWP.md","TDR")]
dat.all <- dat.all[complete.cases(dat.all),]
#------------------------------------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------------------------------------
#- read in parameters, process them.

#------
#- read in the parameter values from table1.csv
params1 <- read.csv("Output/table1.csv")
params1$Species <- factor(rep(c("cacu","eusi","eute","pira"),2))
params1$Yvar <- factor(c(rep("s",4),rep("ns",4)))
params1$Xvar <- factor(c(rep("theta",8)))
params1$Xo <- as.numeric(substr(params1$Xo,start=1,stop=6))
params1$Xh <- as.numeric(substr(params1$Xh,start=1,stop=6))
params1$q <- as.numeric(substr(params1$q,start=1,stop=6))

#- get the low parameter
params_Xo <- reshape2::dcast(params1,Species~Xvar+Yvar,value.var="Xo")
names(params_Xo)[2:3] <- paste("Xl",names(params_Xo)[2:3],sep="_")

#- get the high parameter
params_Xh <- reshape2::dcast(params1,Species~Xvar+Yvar,value.var="Xh")
names(params_Xh)[2:3] <- paste("Xh",names(params_Xh)[2:3],sep="_")

#- get the q parameter
params_q <- reshape2::dcast(params1,Species~Xvar+Yvar,value.var="q")
names(params_q)[2:3] <- paste("q",names(params_q)[2:3],sep="_")

params2 <- merge(params_Xo,params_Xh,by="Species")
params3 <- merge(params2,params_q,by="Species")
#------



#------
#- do the same, but for the exponential parameters in Table 2
#- read in the parameter values from table1.csv
eparams1 <- read.csv("Output/table2.csv")
eparams1$Species <- factor(rep(c("cacu","eusi","eute","pira"),2))
eparams1$Yvar <- factor(c(rep("s",4),rep("ns",4)))
eparams1$Xvar <- factor(c(rep("lwp",8)))
eparams1$a <- as.numeric(substr(eparams1$a,start=1,stop=6))
eparams1$b <- as.numeric(substr(eparams1$b,start=1,stop=6))

#- get the a parameter
params_a <- reshape2::dcast(eparams1,Species~Xvar+Yvar,value.var="a")
names(params_a)[2:3] <- paste("a",names(params_a)[2:3],sep="_")

#- get the b parameter
params_b <- reshape2::dcast(eparams1,Species~Xvar+Yvar,value.var="b")
names(params_b)[2:3] <- paste("b",names(params_b)[2:3],sep="_")

eparams2 <- merge(params_a,params_b,by="Species")
params4 <- merge(params3,eparams2,by="Species")

#------



# read in the tuzet parameters
params.tz <- read.csv("Output/table3.csv")

#- merge parameters together
params <- merge(params4,params.tz,by="Species")
#------------------------------------------------------------------------------------------------------------------



#- merge parameters with data, predict rates
dat.all3 <- merge(dat.all,params,by="Species")

#- calculate the beta terms
dat.all3$Beta_theta_s <- with(dat.all3, ((TDR-Xl_theta_s)/(Xh_theta_s-Xl_theta_s))^q_theta_s)
dat.all3$Beta_theta_ns <- with(dat.all3, ((TDR-Xl_theta_ns)/(Xh_theta_ns-Xl_theta_ns))^q_theta_ns)
#dat.all3$Beta_lwp_s <- with(dat.all3, ((LWP.pd-Xl_lwp_s)/(Xh_lwp_s-Xl_lwp_s))^q_lwp_s)
#dat.all3$Beta_lwp_ns <- with(dat.all3, ((LWP.pd-Xl_lwp_ns)/(Xh_lwp_ns-Xl_lwp_ns))^q_lwp_ns)
dat.all3$Beta_lwp_s <- with(dat.all3, a_lwp_s*exp(b_lwp_s*LWP.pd)/a_lwp_s)
dat.all3$Beta_lwp_ns <- with(dat.all3, a_lwp_ns*exp(b_lwp_ns*LWP.pd)/a_lwp_ns)


#- convert Beta values >1 to 1
dat.all3$Beta_theta_s[which(dat.all3$Beta_theta_s >1)] <- 1
dat.all3$Beta_theta_ns[which(dat.all3$Beta_theta_ns >1)] <- 1
#dat.all3$Beta_lwp_s[which(dat.all3$Beta_lwp_s >1)] <- 1
#dat.all3$Beta_lwp_ns[which(dat.all3$Beta_lwp_ns >1)] <- 1



#- convert values lower than Xl with zero
dat.all3$Beta_theta_s[which(dat.all3$TDR < dat.all3$Xl_theta_s)] <- 0
dat.all3$Beta_theta_ns[which(dat.all3$TDR < dat.all3$Xl_theta_ns)] <- 0
#dat.all3$Beta_lwp_s[which(dat.all3$LWP < dat.all3$Xl_lwp_s )] <- 0
#dat.all3$Beta_lwp_ns[which(dat.all3$LWP < dat.all3$Xl_lwp_ns )] <- 0

#dat.all3 <- dat.all3[complete.cases(dat.all3),] # remove 10 points with missing data
#- convert a very few NA's to zero
#dat.all3$Beta_lwp_ns[which(is.na(dat.all3$Beta_lwp_ns))]  <- 0
#------------------------------------------------------------------------------------------------------------------







#------------------------------------------------------------------------------------------------------------------
#-- predict photosynthetic rates under each modeling condition
#------------------------------------------------------------------------------------------------------------------

#- subset to dry data only?
dat.all3 <- subset(dat.all3,Treat=="dry")
g1value <- 5
Vcmax_value <- 85
Jmax_multi <- 1.6

#- Null (no change in physiology with drought)
dat.all3$P_null <- Photosyn(VPD=dat.all3$VpdL,Ca=dat.all3$CO2S,PPFD=dat.all3$PARi,Tleaf=dat.all3$Tleaf,
                            g1=g1value,Vcmax=Vcmax_value,Jmax=Jmax_multi*Vcmax_value)$ALEAF

#- Beta model based on theta
dat.all3[,c("P_theta_s","Gs_theta_s")] <- Photosyn(VPD=dat.all3$VpdL,Ca=dat.all3$CO2S,PPFD=dat.all3$PARi,Tleaf=dat.all3$Tleaf,
                              g1=g1value*dat.all3$Beta_theta_s,Vcmax=Vcmax_value,Jmax=Jmax_multi*Vcmax_value)[,c("ALEAF","GS")]
dat.all3[,c("P_theta_ns","Gs_theta_ns")] <- Photosyn(VPD=dat.all3$VpdL,Ca=dat.all3$CO2S,PPFD=dat.all3$PARi,Tleaf=dat.all3$Tleaf,
                              g1=g1value,Vcmax=((Vcmax_value-1)*dat.all3$Beta_theta_ns+1),Jmax=(Jmax_multi*(Vcmax_value-1)*dat.all3$Beta_theta_ns+1))[,c("ALEAF","GS")]
dat.all3[,c("P_theta_sns","Gs_theta_sns")] <- Photosyn(VPD=dat.all3$VpdL,Ca=dat.all3$CO2S,PPFD=dat.all3$PARi,Tleaf=dat.all3$Tleaf,
                               g1=g1value*dat.all3$Beta_theta_s,Vcmax=((Vcmax_value-1)*dat.all3$Beta_theta_ns+1),Jmax=(Jmax_multi*(Vcmax_value-1)*dat.all3$Beta_theta_ns+1))[,c("ALEAF","GS")]

#- Beta models based on lwp
dat.all3[,c("P_lwp_s","Gs_lwp_s")]<- Photosyn(VPD=dat.all3$VpdL,Ca=dat.all3$CO2S,PPFD=dat.all3$PARi,Tleaf=dat.all3$Tleaf,
                              g1=g1value*dat.all3$Beta_lwp_s,Vcmax=Vcmax_value,Jmax=Jmax_multi*Vcmax_value)[,c("ALEAF","GS")]
dat.all3[,c("P_lwp_ns","Gs_lwp_ns")] <- Photosyn(VPD=dat.all3$VpdL,Ca=dat.all3$CO2S,PPFD=dat.all3$PARi,Tleaf=dat.all3$Tleaf,
                               g1=g1value,Vcmax=((Vcmax_value-1)*dat.all3$Beta_lwp_ns+1),Jmax=(Jmax_multi*(Vcmax_value-1)*dat.all3$Beta_lwp_ns+1))[,c("ALEAF","GS")]
dat.all3[,c("P_lwp_sns","Gs_lwp_sns")] <- Photosyn(VPD=dat.all3$VpdL,Ca=dat.all3$CO2S,PPFD=dat.all3$PARi,Tleaf=dat.all3$Tleaf,
                           g1=g1value*dat.all3$Beta_lwp_s,Vcmax=((Vcmax_value-1)*dat.all3$Beta_lwp_ns+1),Jmax=(Jmax_multi*(Vcmax_value-1)*dat.all3$Beta_lwp_ns+1))[,c("ALEAF","GS")]


#- Tuzet alone
dat.all3[,c("P_tuzet_s","Gs_tuzet_s")] <- photosyn(SF=dat.all3$Sf, PSIV=dat.all3$psiv, G0=0.005, VCMAX=Vcmax_value,JMAX=Jmax_multi*Vcmax_value,G1=15,
                                                   K=dat.all3$K*exp(-((-1*dat.all3$LWP.pd/dat.all3$b)^dat.all3$c)),
                CS=dat.all3$CO2S,WEIGHTEDSWP=dat.all3$LWP.pd, VPD=dat.all3$VpdL,PAR=dat.all3$PARi,TLEAF=dat.all3$Tleaf)[,c("ALEAF","GS")]
dat.all3[,c("P_tuzet_sns","Gs_tuzet_sns")] <- photosyn(SF=dat.all3$Sf, PSIV=dat.all3$psiv, G0=0.005,G1=15,K=dat.all3$K,
                                 VCMAX=((Vcmax_value-1)*dat.all3$Beta_lwp_ns+1),JMAX=((Vcmax_value-1)*Jmax_multi*dat.all3$Beta_lwp_ns+1),
                               CS=dat.all3$CO2S,WEIGHTEDSWP=dat.all3$LWP.pd, VPD=dat.all3$VpdL,PAR=dat.all3$PARi,TLEAF=dat.all3$Tleaf)[,c("ALEAF","GS")]





#--------------------------------
#- plot observed vs. predicted
windows(25,35)
par(mfrow=c(5,2),oma=c(6,9,4,3),mar=c(0,0,0,0),xpd=F)
colors <- brewer.pal(4,"Set1")
fitcol=alpha("darkgrey",0.8)
textsize <- 1.7
ylims <- c(-2,31)
xlims <- c(-2,35)

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
plot(dat.all3$Photo,dat.all3$P_theta_s,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="",xpd=F);abline(0,1)
predline(lm(P_theta_s~Photo,data=dat.all3),fittype="confidence",col=fitcol,xpd=F)
magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,tcl=0.3)
title(ylab=expression(beta[s]),xpd=NA,cex.lab=textsize)
legend("bottomright",letters[1],bty="n",cex=1.2)
r2val <-round(summary(lm(P_theta_s~Photo,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")


plot(dat.all3$Photo,dat.all3$P_lwp_s,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(P_lwp_s~Photo,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(0,0,1),frame.plot=T,las=1,tcl=0.3)
legend("bottomright",letters[4],bty="n",cex=1.2)
r2val <-round(summary(lm(P_lwp_s~Photo,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")

#- plot B models for non-stomatal limitation (theta and then lwp)
plot(dat.all3$Photo,dat.all3$P_theta_ns,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(P_theta_ns~Photo,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,tcl=0.3)
title(ylab=expression(beta[ns]),xpd=NA,cex.lab=textsize)
legend("bottomright",letters[2],bty="n",cex=1.2)
r2val <-round(summary(lm(P_theta_ns~Photo,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")

plot(dat.all3$Photo,dat.all3$P_lwp_ns,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(P_lwp_ns~Photo,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(0,0,1),frame.plot=T,las=1,tcl=0.3)
legend("bottomright",letters[5],bty="n",cex=1.2)
r2val <-round(summary(lm(P_lwp_ns~Photo,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")

#- plot B models with both
plot(dat.all3$Photo,dat.all3$P_theta_sns,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(P_theta_sns~Photo,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,tcl=0.3)
title(ylab=expression(beta[s+ns]),xpd=NA,cex.lab=textsize)
legend("bottomright",letters[3],bty="n",cex=1.2)
r2val <-round(summary(lm(P_theta_sns~Photo,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")

plot(dat.all3$Photo,dat.all3$P_lwp_sns,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(P_lwp_sns~Photo,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(0,0,1),frame.plot=T,las=1,tcl=0.3)
legend("bottomright",letters[6],bty="n",cex=1.2)
r2val <-round(summary(lm(P_lwp_sns~Photo,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")

#- plot Tuzet models (stomatal, and both). Note the two empty plots to fill space
plot(P_theta_sns~Photo,data=dat.all3,type="n",axes=F,xlab="",ylab="")
plot(dat.all3$Photo,dat.all3$P_tuzet_s ,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(P_tuzet_s~Photo,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(0,1,1),frame.plot=T,las=1,tcl=0.3)
title(ylab=expression(Tuzet[s]),xpd=NA,cex.lab=textsize)
legend("bottomright",letters[7],bty="n",cex=1.2)
r2val <-round(summary(lm(P_tuzet_s~Photo,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")


plot(P_theta_sns~Photo,data=dat.all3,type="n",axes=F,xlab="",ylab="")
plot(dat.all3$Photo,dat.all3$P_tuzet_sns ,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(P_tuzet_sns~Photo,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(1,1,1),frame.plot=T,las=1,tcl=0.3)
title(ylab=expression(Tuzet[s+ns]),xpd=NA,cex.lab=textsize)
legend("bottomright",letters[8],bty="n",cex=1.2)
r2val <-round(summary(lm(P_tuzet_sns~Photo,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")

title(xlab=expression("Observed"~A[sat]~(mu*mol~m^-2~s^-1)),outer=T,cex.lab=2)
title(ylab=expression("Predicted"~A[sat]~(mu*mol~m^-2~s^-1)),outer=T,cex.lab=2,line=5)

#-- calculate concordance correlation coefficients
#null <- summary(agreement(x=dat.all3$Photo,y=dat.all3$P_null,error="constant",TDI_a=20,target="fixed"))
beta_s_theta <- summary(agreement(y=dat.all3$Photo,x=dat.all3$P_theta_s,error="constant",TDI_a=20,target="fixed"))
beta_s_lwp <- summary(agreement(y=dat.all3$Photo,x=dat.all3$P_lwp_s,error="constant",TDI_a=20,target="fixed"))
beta_ns_theta <- summary(agreement(y=dat.all3$Photo,x=dat.all3$P_theta_ns,error="constant",TDI_a=20,target="fixed"))#### work from here
beta_ns_lwp <- summary(agreement(y=dat.all3$Photo,x=dat.all3$P_lwp_ns,error="constant",TDI_a=20,target="fixed"))
beta_sns_theta <- summary(agreement(y=dat.all3$Photo,x=dat.all3$P_theta_sns,error="constant",TDI_a=20,target="fixed"))
beta_sns_lwp <- summary(agreement(y=dat.all3$Photo,x=dat.all3$P_lwp_sns,error="constant",TDI_a=20,target="fixed"))
tuzet_s <- summary(agreement(y=dat.all3$Photo,x=dat.all3$P_tuzet_s,error="constant",TDI_a=20,target="fixed"))
tuzet_sns <- summary(agreement(y=dat.all3$Photo,x=dat.all3$P_tuzet_sns,error="constant",TDI_a=20,target="fixed"))


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
xlims <- c(-0.03,0.72)
ylims <- c(-0.03,0.72)

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
plot(dat.all3$Cond,dat.all3$Gs_theta_s,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="",xpd=F);abline(0,1)
predline(lm(Gs_theta_s~Cond,data=dat.all3),fittype="confidence",col=fitcol,xpd=F)
magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,tcl=0.3)
title(ylab=expression(beta[s]),xpd=NA,cex.lab=textsize)
legend("bottomright",letters[1],bty="n",cex=1.2)
r2val <-round(summary(lm(Gs_theta_s~Cond,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")


plot(dat.all3$Cond,dat.all3$Gs_lwp_s,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(Gs_lwp_s~Cond,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(0,0,1),frame.plot=T,las=1,tcl=0.3)
legend("bottomright",letters[4],bty="n",cex=1.2)
r2val <-round(summary(lm(Cond~Gs_lwp_s,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")

#- plot B models for non-stomatal limitation (theta and then lwp)
plot(dat.all3$Cond,dat.all3$Gs_theta_ns,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(Gs_theta_ns~Cond,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,tcl=0.3)
title(ylab=expression(beta[ns]),xpd=NA,cex.lab=textsize)
legend("bottomright",letters[2],bty="n",cex=1.2)
r2val <-round(summary(lm(Cond~Gs_theta_ns,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")

plot(dat.all3$Cond,dat.all3$Gs_lwp_ns,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(Gs_lwp_ns~Cond,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(0,0,1),frame.plot=T,las=1,tcl=0.3)
legend("bottomright",letters[5],bty="n",cex=1.2)
r2val <-round(summary(lm(Cond~Gs_lwp_ns,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")

#- plot B models with both
plot(dat.all3$Cond,dat.all3$Gs_theta_sns,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(Gs_theta_sns~Cond,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,tcl=0.3)
axis(side=1,las=1,tcl=0.3,at=c(0,0.1,0.2,0.3,0.4,0.5))
title(ylab=expression(beta[s+ns]),xpd=NA,cex.lab=textsize)
legend("bottomright",letters[3],bty="n",cex=1.2)
r2val <-round(summary(lm(Cond~Gs_theta_sns,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")

plot(dat.all3$Cond,dat.all3$Gs_lwp_sns,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(Gs_lwp_sns~Cond,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(0,0,1),frame.plot=T,las=1,tcl=0.3)
legend("bottomright",letters[6],bty="n",cex=1.2)
r2val <-round(summary(lm(Cond~Gs_lwp_sns,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")

#- plot Tuzet models (stomatal, and both). Note the two empty plots to fill space
plot(Gs_theta_sns~Cond,data=dat.all3,type="n",axes=F,xlab="",ylab="")
plot(dat.all3$Cond,dat.all3$Gs_tuzet_s ,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(Gs_tuzet_s~Cond,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(0,1,1),frame.plot=T,las=1,tcl=0.3)
title(ylab=expression(Tuzet[s]),xpd=NA,cex.lab=textsize)
legend("bottomright",letters[7],bty="n",cex=1.2)
r2val <-round(summary(lm(Cond~Gs_tuzet_s,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")


plot(Gs_theta_sns~Cond,data=dat.all3,type="n",axes=F,xlab="",ylab="")
plot(dat.all3$Cond,dat.all3$Gs_tuzet_sns ,xlim=xlims,ylim=ylims,axes=F,xlab="",ylab="");abline(0,1)
predline(lm(Gs_tuzet_sns~Cond,data=dat.all3),fittype="confidence",col=fitcol,xpd=T)
magaxis(side=c(1,2,4),labels=c(1,1,1),frame.plot=T,las=1,tcl=0.3)
title(ylab=expression(Tuzet[s+ns]),xpd=NA,cex.lab=textsize)
legend("bottomright",letters[8],bty="n",cex=1.2)
r2val <-round(summary(lm(Cond~Gs_tuzet_sns,data=dat.all3))$r.squared,2)
legend("topleft",legend= bquote(r^2 == .(r2val)),bty="n")

title(xlab=expression("Observed"~g[s]~(mol~m^-2~s^-1)),outer=T,cex.lab=2)
title(ylab=expression("Predicted"~g[s]~(mol~m^-2~s^-1)),outer=T,cex.lab=2,line=5)




#-- calculate concordance correlation coefficients. Negative data screw up this code.
dat.all4 <- subset(dat.all3,Cond>0 & Gs_theta_s > 0)

beta_s_theta2 <- summary(agreement(x=dat.all4$Cond*1000,y=dat.all4$Gs_theta_s*1000,error="constant",TDI_a=20,target="fixed"))
beta_s_lwp2 <- summary(agreement(x=dat.all4$Cond*1000,y=dat.all4$Gs_lwp_s*1000,error="constant",TDI_a=20,target="fixed"))
beta_ns_theta2 <- summary(agreement(x=dat.all4$Cond*1000,y=dat.all4$Gs_theta_ns*1000,error="constant",TDI_a=20,target="fixed"))
beta_ns_lwp2 <- summary(agreement(x=dat.all4$Cond*1000,y=dat.all4$Gs_lwp_ns*1000,error="constant",TDI_a=20,target="fixed"))
beta_sns_theta2 <- summary(agreement(x=dat.all4$Cond*1000,y=dat.all4$Gs_theta_sns*1000,error="constant",TDI_a=20,target="fixed"))
beta_sns_lwp2 <- summary(agreement(x=dat.all4$Cond*1000,y=dat.all4$Gs_lwp_sns*1000,error="constant",TDI_a=20,target="fixed"))
tuzet_s2 <- summary(agreement(x=dat.all4$Cond*1000,y=dat.all4$Gs_tuzet_s*1000,error="constant",TDI_a=20,target="fixed"))
tuzet_sns2 <- summary(agreement(x=dat.all4$Cond*1000,y=dat.all4$Gs_tuzet_sns*1000,error="constant",TDI_a=20,target="fixed"))


#-- calculate the residual standard deviation
summary(lm(Cond~Gs_theta_s,data=dat.all4))$sigma
summary(lm(Cond~Gs_lwp_s,data=dat.all4))$sigma
summary(lm(Cond~Gs_theta_ns,data=dat.all4))$sigma
summary(lm(Cond~Gs_lwp_ns,data=dat.all4))$sigma
summary(lm(Cond~Gs_theta_sns,data=dat.all4))$sigma
summary(lm(Cond~Gs_lwp_sns,data=dat.all4))$sigma
summary(lm(Cond~Gs_theta_sns,data=dat.all4))$sigma
summary(lm(Cond~Gs_lwp_sns,data=dat.all4))$sigma
summary(lm(Cond~Gs_tuzet_s,data=dat.all4))$sigma
summary(lm(Cond~Gs_tuzet_sns,data=dat.all4))$sigma


