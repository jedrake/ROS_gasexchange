#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
#- This is the central analysis script for a leaf-level gas exchange analysis of potted trees experiencing
#   controlled droughts under large outdoor rainout shelters (ROS). Most of the actual data manipulation occurs
#   in functions in dataFunction.R and plotFunctions.R. These functions are called by this script.
#  
#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

#- load the libraries and scripts that do all of the actual work.
source("R/loadLibraries.R")


#-------------------------------------------------------------------------------------------------------
#- plot soil moisture over time. Figure 1.
plotVWC(ptsize=1.8,output=T,type="1panel") # can plot type="1panel" or "4panel"
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
#- plot the leaf water potential data over time. Figure 2.
plotLWP(fillcol="lightgrey",size=1.75,output=T,labsize=1.8)
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------  
#- make a big 16-panel plot of photosynthetic variables over time. Figure 3.
plotGX(output=T)
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
#- plot the soil moisture release curve (Figure S2).
plotMoistCurve(output=T)
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
#- plot Photo and Cond vs. leaf water potential, to show iso vs. anisohydry (Figure S3).
plotHydry(output=T)
#-------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------
#- plot normalized g1 and non-stomatal limitation as a function of VWC with beta functions.
#  Figure 4 and S3.

#- get the data, fit the beta functions. For some reason this cannot be dropped into a function.
#    Keep it as a script!
source("R/fitBeta_g1_nsl.R")
plotBetasG1NSL(output=T,g1data=g1values,NSLdata=NSLpars,g1list=g1_TDR_beta,NSLlist=NSL_TDR_beta)

source("R/fitBeta_g1_nsl_LWP.R")
plotBetasG1NSL_LWP(output=T,g1data=g1values2,NSLdata=NSLpars2,g1list=g1_TDR_beta2,NSLlist=NSL_TDR_beta2)
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
#- make a table of the beta parameters and standard errors. This table is written out as a csv in "Output"
writeBetaParams()
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
#- plot the C isotope composition of ROS leaves
plotd13C(export=T)
#-------------------------------------------------------------------------------------------------------




#- put together the met data at high resolution
#there was a large zip file on HIEv containing ROS met data. I downloaded it and placed it in ROS/Data/MET
#these files contain data up until 31 March 2013, which isn't quite recent enough for ROS gax exchange,
# which needs data until 16 May 2013
met1.old <- readTOA5(file="Data/MET/ROSWS_MET_20110619-20130331/data/ROS_WS_Table05min_20130331.dat")

#find the more recent met data from HIEv
met1.new <- downloadTOA5(filename="ROS_WS_Table05min",maxnfiles=100,topath="./Data/MET/ROSWS_MET_20110619-20130331/data/")

#combine the old and new data
met1 <- subset(rbind(met1.old,met1.new),Date>as.Date("2012-09-01") & Date < as.Date("2013-05-01"))

#calculate VPD
met1$VPD <- RHtoVPD(RH=met1$RH,TdegC=met1$AirTC_Avg)

#- get the logged soil moisture dataset, merge with the met data
vwc.loggers <- get.TDR.logger(startdate="2011-11-01")
vwc <- summaryBy(VWC~Date+Species,data=subset(data.frame(vwc.loggers),Treat=="dry"),FUN=mean,keep.names=T,na.rm=T)
met2 <- merge(met1,vwc,by="Date")

#- now we have 15-minutely met data for each species. 

#- maximum Vcmax data (see returnVcmaxa(), but get the 75th (95th?) percentile instead)
maxVcmax <- data.frame(Species=c("cacu","eusi","eute","pira"),Vcmax_max=c(82,84,79,81))
#maxVcmax <- data.frame(Species=c("cacu","eusi","eute","pira"),Vcmax_max=c(146,130,130,127))

# met3 <- merge(met2,maxVcmax,by="Species")
# 
# met4 <- met3
# 
# #- get the g1 fits, merge in the max for each species
# g1values <- returng1()
# maxG1 <- summaryBy(g1~Species,data=g1values,FUN=max,keep.names=T)
# met <- merge(met4,maxG1,by="Species")
# 
# 
# 
# 
# #- reduce Vcmax and g1 for each species
# met$g1_adj <- NA
# met$Vcmax_adj <- NA
# 
# cacus <- which(met$Species=="cacu")
# eusis <- which(met$Species=="eusi")
# eutes <- which(met$Species=="eute")
# piras <- which(met$Species=="pira")
# 
# 
# #- function to return beta values of a thing (g1 or apparent Vcmax). Range from 0 to 1.
# #  Accepts (1) a nls model based on a beta function, (2) a vector of values to predict from,
# #   and (3) a type variable refering to the predictor (TDR or LWP)
# predictBeta <- function(model,data,type="TDR"){
#   params <- coef(model)
#   
#   #- find the data lower than Xlow (to return 0)
#   lows <- which(data < unname(params[1]))
# 
#   #- find the data higher than Xhigh (to return 1.0)
#   highs <- which(data > unname(params[2]))
#   
#   #- predict output value for all x-values
#   if(type=="TDR") yvals <- predict(model,newdat=data.frame(TDR=data))
#   if(type=="LWP") yvals <- predict(model,newdat=data.frame(LWP=data))
#   
#   #- override the low and high values in the predicted value
#   yvals[lows] <- 0.1
#   yvals[highs] <- 1
#   
#   return(yvals)
# }
# 
# #- predict g1 and apparent Vcmax based on the fitted beta models
# met$g1_adj[cacus] <- predictBeta(model=fit.spg1[[1]],data=met$VWC[cacus],type="TDR")*met$g1[cacus]
# met$g1_adj[eusis] <- predictBeta(model=fit.spg1[[2]],data=met$VWC[eusis],type="TDR")*met$g1[eusis]
# met$g1_adj[eutes] <- predictBeta(model=fit.spg1[[3]],data=met$VWC[eutes],type="TDR")*met$g1[eutes]
# met$g1_adj[piras] <- predictBeta(model=fit.spg1[[4]],data=met$VWC[piras],type="TDR")*met$g1[piras]
# 
# met$Vcmax_adj[cacus] <- predictBeta(model=fit.spNSL[[1]],data=met$VWC[cacus],type="TDR")*met$Vcmax_max[cacus]
# met$Vcmax_adj[eusis] <- predictBeta(model=fit.spNSL[[2]],data=met$VWC[eusis],type="TDR")*met$Vcmax_max[eusis]
# met$Vcmax_adj[eutes] <- predictBeta(model=fit.spNSL[[3]],data=met$VWC[eutes],type="TDR")*met$Vcmax_max[eutes]
# met$Vcmax_adj[piras] <- predictBeta(model=fit.spNSL[[4]],data=met$VWC[piras],type="TDR")*met$Vcmax_max[piras]
# 
# 
# 
# #- model photo in the absence of any drought effect
# met$ALEAF <- Photosyn(VPD=met$VPD,PPFD=1800,Tleaf=met$AirTC_Avg,Vcmax=met$Vcmax_max,Jmax=0.8*met$Vcmax_max,
#                       g1=met$g1)$ALEAF
# 
# #- model photo with a drought effect on g1
# met$ALEAF_g1 <- Photosyn(VPD=met$VPD,PPFD=1800,Tleaf=met$AirTC_Avg,Vcmax=met$Vcmax_max,Jmax=0.8*met$Vcmax_max,
#                       g1=met$g1_adj)$ALEAF
# #- model photo with a drought effect on Vcmax
# met$ALEAF_NSL <- Photosyn(VPD=met$VPD,PPFD=1800,Tleaf=met$AirTC_Avg,Vcmax=met$Vcmax_adj,Jmax=0.8*met$Vcmax_adj,
#                          g1=met$g1)$ALEAF
# #- model photo with a drought effect on both
# met$ALEAF_both <- Photosyn(VPD=met$VPD,PPFD=1800,Tleaf=met$AirTC_Avg,Vcmax=met$Vcmax_adj,Jmax=0.8*met$Vcmax_adj,
#                          g1=met$g1_adj)$ALEAF
# 
# #- plot daily averages
# 
# met.d <- summaryBy(PPFD_Avg+AirTC_Avg+VPD+VWC+ALEAF+ALEAF_g1+ALEAF_NSL+ALEAF_both~Date+Species,FUN=mean,data=met,keep.names=T)
# met.d2 <- summaryBy(VWC~Date,FUN=mean,data=met,keep.names=T)
# met.d.l <- split(met.d,met.d$Species)
# windows();par(mfrow=c(5,1),mar=c(0,0,0,0),oma=c(5,5,1,1))
# 
# colors=c("black","forestgreen","blue","red")
# startdate <- as.Date("2013-01-18")
# plot(VWC~Date,data=subset(met.d2,Date>startdate),type="l",col="black",legend=F)
# 
# for (i in 1:length(met.d.l)){
#   toplot <- subset(met.d.l[[i]],Date>startdate)
#   plot(ALEAF~Date,data=toplot,type="l",col=colors[1],legend=F,ylim=c(5,15))
#   lines(ALEAF_g1~Date,data=toplot,type="l",col=colors[2],legend=F)
#   lines(ALEAF_NSL~Date,data=toplot,type="l",col=colors[3],legend=F)
#   lines(ALEAF_both~Date,data=toplot,type="l",col=colors[4],legend=F)
#   
# }
# 
# 
# legend("top",)
# 
# #- model photo with a drought effect on both








#- take an alternative approach, and predict values across a wide range of VWC values, along with Tair of 25, VPD of 1.5
predData.1 <- expand.grid(VWC=seq(0.001,0.4,length=101),Species=levels(met.d$Species))
predData.2 <- merge(predData.1,maxVcmax)
predData <- merge(predData.2,maxG1)

cacus <- which(predData$Species=="cacu")
eusis <- which(predData$Species=="eusi")
eutes <- which(predData$Species=="eute")
piras <- which(predData$Species=="pira")


#- predict g1 and apparent Vcmax based on the fitted beta models
predData$g1_adj[cacus] <- predictBeta(model=fit.spg1[[1]],data=predData$VWC[cacus],type="TDR")*predData$g1[cacus]
predData$g1_adj[eusis] <- predictBeta(model=fit.spg1[[2]],data=predData$VWC[eusis],type="TDR")*predData$g1[eusis]
predData$g1_adj[eutes] <- predictBeta(model=fit.spg1[[3]],data=predData$VWC[eutes],type="TDR")*predData$g1[eutes]
predData$g1_adj[piras] <- predictBeta(model=fit.spg1[[4]],data=predData$VWC[piras],type="TDR")*predData$g1[piras]

predData$Vcmax_adj[cacus] <- predictBeta(model=fit.spNSL[[1]],data=predData$VWC[cacus],type="TDR")*predData$Vcmax_max[cacus]
predData$Vcmax_adj[eusis] <- predictBeta(model=fit.spNSL[[2]],data=predData$VWC[eusis],type="TDR")*predData$Vcmax_max[eusis]
predData$Vcmax_adj[eutes] <- predictBeta(model=fit.spNSL[[3]],data=predData$VWC[eutes],type="TDR")*predData$Vcmax_max[eutes]
predData$Vcmax_adj[piras] <- predictBeta(model=fit.spNSL[[4]],data=predData$VWC[piras],type="TDR")*predData$Vcmax_max[piras]


#- model photo in the absence of any drought effect
predData$ALEAF <- Photosyn(VPD=1.5,PPFD=1800,Tleaf=25,Ca=400,Vcmax=predData$Vcmax_max,Jmax=1.6*predData$Vcmax_max,
                      g1=predData$g1)$ALEAF

#- model photo with a drought effect on g1
predData$ALEAF_g1 <- Photosyn(VPD=1.5,PPFD=1800,Tleaf=25,Ca=400,Vcmax=predData$Vcmax_max,Jmax=1.6*predData$Vcmax_max,
                         g1=predData$g1_adj)$ALEAF
#- model photo with a drought effect on Vcmax
predData$ALEAF_NSL <- Photosyn(VPD=1.5,PPFD=1800,Tleaf=25,Ca=400,Vcmax=predData$Vcmax_adj,Jmax=1.6*predData$Vcmax_adj,
                          g1=predData$g1)$ALEAF
predData$ALEAF_both <- Photosyn(VPD=1.5,PPFD=1800,Tleaf=25,Ca=400,Vcmax=predData$Vcmax_adj,Jmax=1.6*predData$Vcmax_adj,
                           g1=predData$g1_adj)$ALEAF

#- get the data to overlay
Adat <- return.gx.vwc()
Adat.m <- summaryBy(Photo+TDR~Date+Species+Treat,data=Adat,keep.names=T)
Adat.m$VWC <- Adat.m$TDR/100


# windows()
# #plotBy(ALEAF~VWC|Species,data=predData,type="l",lwd=2,lty=1,ylim=c(-1,15),legend=F)
# plotBy(ALEAF_NSL~VWC|Species,data=predData,type="l",lwd=2,lty=2,add=F,ylim=c(-1,20),legend=F)
# plotBy(ALEAF_g1~VWC|Species,data=predData,type="l",lwd=2,lty=3,add=T,legend=F)
# plotBy(ALEAF_both~VWC|Species,data=predData,type="l",lwd=2,lty=1,add=T,legend=F)
# 

#- plot data with predictions
windows();par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(7,7,2,2))
predData.l <- split(predData,predData$Species)
for (i in 1:length(predData.l)){
  toplot <- predData.l[[i]]
  plotBy(ALEAF_NSL~VWC,data=toplot,type="l",lwd=2,lty=2,add=F,ylim=c(-1,25),legend=F,axes=F)
  plotBy(ALEAF_g1~VWC,data=toplot,type="l",lwd=2,lty=3,add=T,legend=F)
  plotBy(ALEAF_both~VWC,data=toplot,type="l",lwd=2,lty=1,add=T,legend=F)
  points(Photo~VWC,data=subset(Adat.m,Species==as.character(toplot$Species[1])),cex=1.5,col="black",pch=16)
  
  title(main=toplot$Species[1],line=-1.5)
  
  if(i==1) magaxis(side=c(1:4),labels=c(0,1,0,0),frame.plot=T,las=1)
  if(i==2) magaxis(side=c(1:4),labels=c(0,0,0,2),frame.plot=T,las=1)
  if(i==3) magaxis(side=c(1:4),labels=c(1,1,0,0),frame.plot=T,las=1)
  if(i==4) magaxis(side=c(1:4),labels=c(1,0,0,1),frame.plot=T,las=1)
  
}
title(ylab=expression(A[sat]~(mu*mol~m^-2~s^-1)),outer=T,cex.lab=2)
title(xlab=expression(Soil~volumetric~water~content~(theta~";"~m^3~m^-3)),outer=T,cex.lab=2)
dev.copy2pdf(file="Output/Asat_VWC_modelPredictions.pdf")
