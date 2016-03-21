source("R/loadLibraries.R")



#------------------------------------------------------------------------------------------------------------------
#- make the soil moisture plot
#------------------------------------------------------------------------------------------------------------------
plotVWC <- function(ptsize=1.5,output=T,type="4panel"){
  
  #- get the VWC data from the loggers and handheld
  vwc.loggers <- get.TDR.logger(startdate="2011-11-01")
  vwc.hand <- get.TDR.handheld(filterBlanks=1)
  vwc.hand$Species <- droplevels(vwc.hand$Species)
  
  #----
  # remove the handheld tdr measurements with <3 observations per treatment (some dry pots were mislabeled as wet)
  #writing my own length function that omits NA's
  length2 <- function(x,...){
    length(na.omit(x))
  }
  vwc.hand0 <- summaryBy(TDR~Date+Species+Treat,FUN=c(mean,length2,standard.error),dat=vwc.hand,keep.names=FALSE,na.rm=TRUE)
  vwc.hand1 <- subset(vwc.hand0,TDR.length2>3 | Date>as.Date("2013-05-01") ) #remove dates that have <3 observations (some pots were miscoded)
  names(vwc.hand1)[4] <- "TDR"
  names(vwc.hand1)[6] <- "TDR.se"
  
    #-----
  # average across dates for dry pots only and have a look. Try to identify the dates where plants were re-watered
  vwc.dates.dry <- summaryBy(TDR~Date,data=subset(vwc.hand1,Treat=="dry"))
  
  # 
  # #------------------------------------------------------------------------------------------------------------------
  # #Write out a csv of the TDR data for sharing
  # vwc.hand.out <- vwc.hand[,c(3,7,8,4,5,6)];names(vwc.hand.out)[5] <- "VWC"
  # vwc.hand.out$Treat <- as.factor(ifelse(vwc.hand.out$Treat=="dry","Dry","Wet"))
  # vwc.hand.out$VWC <- vwc.hand.out$VWC/100
  # write.csv(vwc.hand.out,file="Output/ROS_MD_PM_SOILMOIST-HANDHELD_L1.csv",row.names=F)
  
  names(vwc.loggers)[3] <- "Pot"
  vwc.loggers.out <- vwc.loggers[,c(6,5,4,3,1,2,7)]
  vwc.loggers.out$Treat <- as.factor(ifelse(vwc.loggers.out$Treat=="dry","Dry","Wet"))
  #write.csv(vwc.loggers.out,file="./Output/ROS_MD_PM_SOILMOIST-LOGGERS_MA.csv",row.names=FALSE)
  # #------------------------------------------------------------------------------------------------------------------
  
  
  
  
  #------------------------------------------------------------------------------------------------------------------
  #get gas exchange data for the gasex dates
  gx <- read.gx.ROS()
  gxdates <- as.Date(levels(as.factor(gx$Date)))
  
  
  #------------------------------------------------------------------------------------------------------------------
  #calculate treatment averages in soil moisture, plot 4-panel plot
  
  vwc.loggers2 <- summaryBy(VWC~Date+Species+Treat,data=as.data.frame(vwc.loggers),FUN=c(mean,standard.error), na.rm=TRUE)
  #vwc.hand2 <- summaryBy(TDR~Date+Species+Treat,data=subset(vwc.hand1,Date!=as.Date("2013-01-18")),FUN=c(mean), na.rm=TRUE)
  vwc.hand2 <- subset(vwc.hand1,Date!=as.Date("2013-01-18") & is.na(TDR) ==F)
  
  vwc.log.list <- split(vwc.loggers2,vwc.loggers2$Species)
  vwc.hand.list <- split(vwc.hand2,vwc.hand2$Species)
  
  
  #---- make the plot, if it's the 4-panel option
  if(type=="4panel"){
    windows(20,20)
    par(mfrow=c(2,2),oma=c(5,5,0.5,2),mar=c(0,0,0,0))
    palette(c("white","black"))
    colors <- c("white","black")
    
    startdate <- as.Date(x="2012-10-1",format="%Y-%m-%d")
    enddate <- as.Date(x="2013-04-30",format="%Y-%m-%d")
    dates <- as.Date(c("2012-10-8","2012-11-15","2012-11-25","2012-12-27","2013-1-26","2013-5-17"),format="%Y-%m-%d")
    for (i in 1:length(vwc.log.list)){
      dat1 <- vwc.log.list[[i]]
      dat2 <- vwc.hand.list[[i]]
      
      plot.new()
      plot.window(xlim=c(startdate,enddate),ylim=c(0,0.4))
      
      rect(xleft=dates[1],ybottom=-0.2,xright=dates[2],ytop=0.42,col="lightgrey") #add rectangles for droughts
      rect(xleft=dates[3],ybottom=-0.2,xright=dates[4],ytop=0.42,col="lightgrey") #add rectangles for droughts
      rect(xleft=dates[5],ybottom=-0.2,xright=enddate+13,ytop=0.42,col="lightgrey") #add rectangles for droughts
      
      #if (i==1){
      #  points(VWC.mean~Date,data=subset(dat1,Treat=="dry"),type="b",pch=21,cex=0.5,ylim=c(0,0.4),col="black",bg="white",legend=F,axes=FALSE)
      # points(VWC.mean~Date,data=subset(dat1,Treat=="wet"),type="b",pch=21,cex=0.5,ylim=c(0,0.4),col="black",bg="black",legend=F,axes=FALSE)
      lines(VWC.mean~Date,data=subset(dat1,Treat=="dry"),lwd=1,lty=3,col="black")
      lines(VWC.mean~Date,data=subset(dat1,Treat=="wet"),lwd=1,lty=1,col="black")
      
      #}   
      #if (i > 1)plotBy(VWC.mean~Date|Treat,data=dat1,type="b",cex=0.5,ylim=c(0,0.4),legend=FALSE,axes=FALSE,add=T)
      if (i%%2==1) magaxis(c(2,4),labels=c(1,0),frame.plot=TRUE)
      if (i%%2==0) magaxis(c(2,4),labels=c(0,1),frame.plot=TRUE)
      points(TDR/100~Date,col="black",data=subset(dat2,Treat=="dry"),pch=21,bg="white",cex=ptsize)
      points(TDR/100~Date,col="black",data=subset(dat2,Treat=="wet"),pch=21,bg="black",cex=ptsize)
      
      rug(gxdates,lwd=3,line=-0.5)
      
      #labeling and legends
      legend("topright",legend=dat1$Species[1],bty="n",cex=1.3)
      if (i<3)axis.Date(side=1,at=seq(startdate,enddate,by="month"),format="%m/%y",labels=FALSE)
      if (i>2)axis.Date(side=1,at=seq(startdate,enddate,by="month"),format="%m/%y")
      
      
    }
    title(xlab="Date",outer=TRUE,ylab=expression(Volumetric~water~content~(m^3~m^-3)),cex.lab=2,line=2.5)
    if(output==T) dev.copy2pdf(file="./Output/VWC_ROS_jed.pdf")
  }
  
  
  
  
  
  
  #------------------------------------------------------------------------------------------------------------------
  #--- make another plot, with just the handheld data on one panel
  if(type=="1panel"){
    windows(20,20)
    par(mfrow=c(1,1),oma=c(1,3,0,0),mar=c(4,4,1,3))
    palette(c("white","black"))
    colors <- c("white","black")
    startdate <- as.Date(x="2012-10-1",format="%Y-%m-%d")
    enddate <- as.Date(x="2013-04-30",format="%Y-%m-%d")
    dates <- as.Date(c("2012-10-8","2012-11-8","2012-11-21","2012-12-27","2013-1-26","2013-5-17"),format="%Y-%m-%d")
    pchs=c(21,22,23,24)
    
    plot.new()
    plot.window(xlim=c(startdate,enddate),ylim=c(0,0.4))
    rect(xleft=dates[1],ybottom=-0.2,xright=dates[2],ytop=0.42,col="lightgrey") #add rectangles for droughts
    rect(xleft=dates[3],ybottom=-0.2,xright=dates[4],ytop=0.42,col="lightgrey") #add rectangles for droughts
    rect(xleft=dates[5],ybottom=-0.2,xright=enddate+13,ytop=0.42,col="lightgrey") #add rectangles for droughts
    
    for (i in 1:length(vwc.hand.list)){
      dat2 <- vwc.hand.list[[i]]
      
      adderrorbars(x=dat2$Date,y=dat2$TDR/100,SE=dat2$TDR.se/100,direction="updown")
      
      points(TDR/100~Date,col="black",data=subset(dat2,Treat=="dry"),pch=pchs[i],bg="white",cex=ptsize)
      points(TDR/100~Date,col="black",data=subset(dat2,Treat=="wet"),pch=pchs[i],bg="black",cex=ptsize)
      
      
    }
    magaxis(c(2,4),labels=c(1,1),las=1,cex.axis=1.5,frame.plot=TRUE)
    #rug(gxdates,lwd=3,line=-0.5)
    axis.Date(side=1,at=seq(startdate,enddate,by="month"),format="%m/%y",cex.axis=1.5,tcl=0.5)
    abline(h=0.1,lty=3)
    title(xlab="Date",outer=TRUE,ylab=expression(Volumetric~water~content~(m^3~m^-3)),cex.lab=2,line=-1)
    legend(x=as.Date("2013-1-1"),y=0.416,legend=c("  Cacu","  Eusi","  Eute","  Pira"),pch=c(pchs),
           col=c(rep("black",4)),ncol=2,pt.bg="black",cex=1.4,bg="white")
    legend(x=as.Date("2013-1-8"),y=0.416,legend=c("        ","        ","          ","           "),pch=c(pchs),
           col=c(rep("black",4)),ncol=2,pt.bg="white",cex=1.4,bty="n")
    rug(gxdates,lwd=3,line=-0.5)
    
    if(output==T) dev.copy2pdf(file="./Output/VWC_ROS_jed_1panel.pdf")
  }
  #------------------------------------------------------------------------------------------------------------------
  #------------------------------------------------------------------------------------------------------------------
  
}



