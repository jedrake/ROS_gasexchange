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
#- plot soil moisture over time
plotVWC(ptsize=1.8,output=T,type="1panel") # can plot type="1panel" or "4panel"
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------  
#- make a big 16-panel plot of photosynthetic variables over time
plotGX(output=T)
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
#- plot the leaf water potential data over time
plotLWP(fillcol="lightgrey",size=1.75,output=T,labsize=1.8)
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
#- plot normalized g1 and non-stomatal limitation as a function of VWC with beta functions

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

#- maximum Vcmax data (see returnVcmaxa(), but get the 75 percentile instead)
maxVcmax <- data.frame(Species=c("cacu","eusi","eute","pira"),Vcmax_max=c(106.3,97.3,104.7,97.9))
met3 <- merge(met2,maxVcmax,by="Species")
met4 <- met3

#- get the g1 fits, merge in the max for each species
g1values <- returng1()
maxG1 <- summaryBy(g1~Species,data=g1values,FUN=max,keep.names=T)
met <- merge(met4,maxG1,by="Species")

#- model photo in the absence of any drought effect
met$ALEAF <- Photosyn(VPD=met$VPD,PPFD=met$PPFD_Avg,Tleaf=met$AirTC_Avg,Vcmax=met$Vcmax_max,Jmax=0.8*met$Vcmax_max,
                      g1=met$g1)$ALEAF


#- reduce Vcmax and g1 for each species
met$g1_adj <- NA

cacus <- which(met$Species=="cacu")
eusis <- which(met$Species=="eusi")
eutes <- which(met$Species=="eute")
piras <- which(met$Species=="pira")

met$g1_adj[cacus] <- predict(fit.spg1[[1]],newdat=data.frame(TDR=met$VWC[cacus]))*met$g1[cacus]
met$g1_adj[eusis] <- predict(fit.spg1[[2]],newdat=data.frame(TDR=met$VWC[eusis]))*met$g1[eusis]
met$g1_adj[eutes] <- predict(fit.spg1[[3]],newdat=data.frame(TDR=met$VWC[eutes]))*met$g1[eutes]
met$g1_adj[piras] <- predict(fit.spg1[[4]],newdat=data.frame(TDR=met$VWC[piras]))*met$g1[piras]

#- model photo with a drought effect on g1
met$ALEAF_g1 <- Photosyn(VPD=met$VPD,PPFD=met$PPFD_Avg,Tleaf=met$AirTC_Avg,Vcmax=met$Vcmax_max,Jmax=0.8*met$Vcmax_max,
                      g1=met$g1_adj)$ALEAF

#- plot daily averages
windows()
met.d <- summaryBy(PPFD_Avg+AirTC_Avg+VPD+ALEAF+ALEAF_g1~Date+Species,FUN=mean,data=subset(met,PPFD_Avg>0),keep.names=T)
plotBy(ALEAF~Date|Species,data=met.d,type="l",col=c("black","red","blue","forestgreen"),legend=F)
plotBy(ALEAF_g1~Date|Species,data=met.d,type="l",col=c("black","red","blue","forestgreen"),legend=F,lty=2,add=T)

legend("top",)
      