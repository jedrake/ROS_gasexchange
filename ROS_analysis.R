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
plotLWP(fillcol="lightgrey",size=1.75,output=T,labsize=1.8,output=T)


#- plot normalized g1 and non-stomatal limitation as a function of VWC with beta functions
plotBetasG1NSL(output=T)




#-------------------------------------------------------------------------------------------------------------------------------------
#-- isotopic analysis of plants at the ROS.

#-- notes: As a reminder, the samples ending -21 were collected on 21/12/13, after the second dry down, 
#        while those ending 26 were collected at the end of the final dry down (26/3/14). If you recall, 
#        we tagged the twigs just before the first dry down began, so the first set of “leaves” sampled had formed 
#        at some point during - or between - the first two dry downs, while those collected later had developed 
#        under even more varied (and extreme) conditions.

# I found two possible mislabellings: One of the wet cacu in shelter 3 was labeled 318 on the first sampling and 378 on the second.
# One of the wet Eute in shelter 5 was labeled 505 on the first sampling and 503 on the second. These were confirmed with Sally
# and changed in the original dataset (the correct labels were 318 and 503)

#---- read in the data, do a little processing
d1 <- get_d13C()

#--------------------------------------------------------------
# plots
d1$bigDelta <- d1$bigDelta*1000
d1.m <- summaryBy(deltaC+bigDelta~Species+Treat,data=d1,FUN=c(mean,standard.error))

#-- barplots of bigdelta
windows(18,12);par(mfrow=c(2,1),mar=c(3,5,1,1),oma=c(1,1,3,1))

bp2 <- barplot(height=d1.m$deltaC.mean[1:8],col=c("darkgrey","white"),ylim=c(-35,-25),xpd=F,axes=F)
adderrorbars(x=bp2,y=d1.m$deltaC.mean[1:8],SE=d1.m$deltaC.standard.error[1:8],direction="updown")
magaxis(side=c(2,4),labels=c(1,0),frame.plot=T,las=1)
text(x=c(1.25,3.7,6.1,8.5),y=17.5,xpd=T,labels=c("Cacu","Eusi","Eute","Pira"),cex=1.5)
title(ylab=expression(delta^13~C),cex.lab=1.5,line=1.8)
text(x=c(1.25,3.7,6.1,8.5),y=-35.8,xpd=T,labels=c("Cacu","Eusi","Eute","Pira"),cex=1.5)
legend("bottomright",xpd=NA,legend=c("Wet","Dry"),fill=c("darkgrey","white"),bty="n",ncol=1,cex=1.5)


bp1 <- barplot(height=d1.m$bigDelta.mean[1:8],col=c("darkgrey","white"),ylim=c(18,25),xpd=F,axes=F)
adderrorbars(x=bp1,y=d1.m$bigDelta.mean[1:8],SE=d1.m$bigDelta.standard.error[1:8],direction="updown")
magaxis(side=c(2,4),labels=c(1,0),frame.plot=T,las=1)
text(x=c(1.25,3.7,6.1,8.5),y=17.5,xpd=T,labels=c("Cacu","Eusi","Eute","Pira"),cex=1.5)
title(ylab=expression(Delta~"*"~10^3),cex.lab=1.5,line=1.8)

#dev.copy2pdf(file="Output/ROS_bigdelta_bars.pdf")

      