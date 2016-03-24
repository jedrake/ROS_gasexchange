
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
  
}
#------------------------------------------------------------------------------------------------------------------








#------------------------------------------------------------------------------------------------------------------
#- function to make the plot of gas exchange data over time. This function is called by plotGX().
#------------------------------------------------------------------------------------------------------------------
plot.gs.ros <- function(dat,toplot,toplotse,ylims,ylabs,xlab){ #function needs arguments for the dataframe, and
  #the COLUMN NUMBERS for the variables to plot
  startdate <- as.Date(x="2012-09-30",format="%Y-%m-%d")
  enddate <- as.Date(x="2013-05-30",format="%Y-%m-%d")
  dates <- as.Date(c("2012-10-18","2012-11-15","2012-11-25","2012-12-27","2013-1-26","2013-5-17"),format="%Y-%m-%d")
  
  dat2 <- subset(dat,is.na(dat[,toplot])==F)
  dat.wet <- subset(dat2,Treat=="wet")
  dat.dry <- subset(dat2,Treat=="dry")
  
  
  
  plot.new()
  plot.window(xlim=c(startdate,enddate),ylim=ylims)
  
  rect(xleft=dates[1],ybottom=-200,xright=dates[2],ytop=4000,col="lightgrey") #add rectangles for droughts
  rect(xleft=dates[3],ybottom=-200,xright=dates[4],ytop=4000,col="lightgrey") #add rectangles for droughts
  rect(xleft=dates[5],ybottom=-200,xright=enddate+13,ytop=4000,col="lightgrey") #add rectangles for droughts
  
  if(ylims[1] < 0) abline(h=0,lty=2)
  
  #plot error bars
  plotCI(x=dat.wet$gxDate,y=dat.wet[,toplot],uiw=dat.wet[,toplotse],add=T)
  plotCI(x=dat.dry$gxDate,y=dat.dry[,toplot],uiw=dat.dry[,toplotse],add=T)
  
  #plot big points over the error bars
  points(dat.wet[,toplot]~dat.wet$gxDate,type="b",pch=21,col="black",bg="black",lty=1,add=T,cex=2)
  points(dat.dry[,toplot]~dat.dry$gxDate,type="b",pch=21,col="black",bg="white",lty=2,add=T,cex=2)
  
  
  
  #plotBy(toplot~gxDate|Treat,data=dat,type="b",xlim=c(startdate,enddate),legend=F,
  #       pch=16,cex=1.2,add=T,col=colors,legendwhere="topleft",axes=FALSE,ylim=ylims)
  #magaxis(c(2,4),labels=ylabs,box=TRUE,cex.axis=1.5)
  box()
  axis(side=2,labels=ylabs[1],cex.axis=1.5,tck=0.025,las=2)
  axis(side=4,labels=ylabs[2],cex.axis=1.5,tck=0.025,las=2)
  
  
  #points(TDR.mean/100~Date,col=Treat,data=dat2,pch=16,cex=1.3)
  #rug(gxdates,lwd=3,line=-0.5)
  #labeling and legends
  #legend("topright",legend=dat$Species[1],bty="n",cex=1.3)
  if (xlab==T)axis.Date(side=1,at=seq(startdate,enddate,by=60),tck=0.025,format="%m/%y",labels=T,cex.axis=1.7,las=2)
  if (xlab==F)axis.Date(side=1,at=seq(startdate,enddate,by=60),tck=0.025,format="%m/%y",labels=F,cex.axis=1.7,las=3)
  axis.Date(side=3,at=seq(startdate,enddate,by=60),tck=0.025,format="%m/%y",labels=F,cex.axis=1.7)
}
#------------------------------------------------------------------------------------------------------------------






#------------------------------------------------------------------------------------------------------------------
#- main plotting script for Asat, gs, WUE, and Ci/Ca
#------------------------------------------------------------------------------------------------------------------
plotGX <- function(output=F){
  #- get the gas exchange and handheld TDR data
  ros <- return.gx.vwc()
  
  
  #- get the lwp data and process a bit
  lwp <- read.lwp()
  
  names(lwp)[3] <- "LWPdate"
  lwp.pd <- subset(lwp,Type=="PD" & is.na(LWP)==FALSE)
  lwp.pd$LWP <- -1*lwp.pd$LWP
  lwp.md <- subset(lwp,Type=="MD" & is.na(LWP)==FALSE)
  lwp.md$LWP <- -1*lwp.md$LWP
  
  names(lwp.md)[2] <- "LWP.md"
  lwp2 <- merge(lwp.pd,lwp.md,by=c("Pot","LWPdate","Treat","Species"),all=FALSE)
  lwp2$diff <- with(lwp2,abs(LWP.md)-abs(LWP))
  
  lwp.trt<- summaryBy(LWP+LWP.md+diff~LWPdate+Species+Treat,data=lwp2,FUN=c(mean,standard.error), na.rm=TRUE)
  lwp.trt.list <- split(lwp.trt,lwp.trt$Species)
  #------------------------------------------------------------------------------------------------------------------
  
  
  
  #do some data manipulation, set graphing parameters 
  ros$ITE <- with(ros,Photo/Trmmol)
  ros$WUE <- with(ros,Photo/Cond)
  
  ITE.trt <- summaryBy(Photo+Cond+ITE+WUE+TDR+Ci.Ca~gxDate+Species+Treat,data=ros,FUN=c(mean,standard.error), na.rm=TRUE)
  ITE.trt.list <- split(ITE.trt,ITE.trt$Species)
  
  
  #------------------------------------------------------------------------------------------------------------------
  #set up 16-panel plot, to plot Asat, Cond, WUE, and Ci/Ca over time for all 4 species
  dat.trt <- ITE.trt[,c(1:5,7,9,10,11,13,15)]
  
  
  
  
  dat.trt.list <- split(dat.trt,dat.trt$Species)
  ylabs.list <- list(c(T,F),c(F,F),c(F,F),c(F,T))
  ylims.list <- list(c(0,27),c(0,0.6),c(-20,250),c(0,1.5))
  xlabs <- list(F,F,F,T)
  toplot <- c(4,5,6,7)
  ses <- c(8:11)
  
  labels <- list(expression(A[sat]),expression(g[s]),
                 expression(WUE),expression(C[i]~"/"~C[a]))
  
  #-- get the average WUE for the first two droughts
  firstdrs <- subset(ITE.trt,gxDate==as.Date("2012-10-31") | gxDate==as.Date("2012-11-6") | gxDate==as.Date("2012-12-19"))
  summaryBy(WUE.mean~Species+Treat,data=firstdrs)
  
  windows(14,12)
  par(mfrow=c(4,4),oma=c(8,7,2,5),mar=c(0.25,0.25,0.25,0.25))
  for (i in 1:4){
    #print(i)
    #plot photosynthesis
    for (j in 1:4){
      #print(j)
      plot.gs.ros(dat=dat.trt.list[[j]],toplot=toplot[i],toplotse=ses[i],
                  ylims=ylims.list[[i]],ylabs=ylabs.list[[j]],xlab=xlabs[i])
      
    }
    title(xlab="",outer=TRUE,ylab=labels[[i]],cex.lab=2,line=4,adj=((4-i)/4+0.13))
    if (i==1) title(main=expression(Cacu~~~~~~~~~~~~~~~~~~~~~~Eusi~~~~~~~~~~~~~~~~~~~~~~Eute~~~~~~~~~~~~~~~~~~~~~Pira),outer=T,cex.main=2)
  }
  title(xlab="Date",outer=TRUE,ylab="",cex.lab=2,line=5)
  if(output==T) dev.copy2pdf(file="Output/Asat_Gs_WUE_CiCa_ROS.pdf")
}
#---------------------------------------------------------------------------------------------------------------------






#------------------------------------------------------------------------------------------------------------------
#- plot leaf water potentials over time
#------------------------------------------------------------------------------------------------------------------
plotLWP <- function(fillcol="lightgrey",size=1.75,output=F,labsize=1.8){
  
  ros2 <- return.gx.vwc.lwp()
  
  
  
  
  #---------------------------------------------------------------------------------------------------------------------
  #--- plot pre-dawn and mid-day leaf water potentials over time
  lwp.m <- summaryBy(LWP+LWP.md+diff~gxDate+Treat+Species,FUN=c(mean,standard.error),data=subset(ros2,gxDate>as.Date("2012-11-1")))
  lwp.m.l <- split(lwp.m,lwp.m$Species)
  
  
  windows(20,20)
  par(mfrow=c(2,1),oma=c(4,4,0.5,2),mar=c(2,3,1,1),cex.axis=1.5,cex.lab=1.8)
  
  startdate <- as.Date(x="2012-10-1",format="%Y-%m-%d")
  enddate <- as.Date(x="2013-06-1",format="%Y-%m-%d")
  dates <- as.Date(c("2012-10-18","2012-11-15","2012-11-25","2012-12-27","2013-1-26","2013-5-17"),format="%Y-%m-%d")
  #colors <- rev(c("black","darkgrey"))
  colors <- c("white","black")
  #pchs <- c(15,16,17,18)
  pchs <- c(21,22,23,24)
  
  
  #- plot pre-dawn LWP over time for each species
  plot.new()
  plot.window(xlim=c(startdate,enddate),ylim=c(-10,0))
  rect(xleft=dates[1],ybottom=-45,xright=dates[2],ytop=270,col=fillcol) #add rectangles for droughts
  rect(xleft=dates[3],ybottom=-45,xright=dates[4],ytop=270,col=fillcol) #add rectangles for droughts
  rect(xleft=dates[5],ybottom=-45,xright=enddate+13,ytop=270,col=fillcol) #add rectangles for droughts
  
  #plot pre-dawn
  for (i in 1:length(lwp.m.l)){
    dat1 <- lwp.m.l[[i]]
    
    adderrorbars(x=dat1$gxDate,y=dat1$LWP.mean,SE=dat1$LWP.standard.error,direction="updown")
    points(LWP.mean~gxDate,col="black",data=subset(dat1,Treat=="dry"),pch=pchs[i],bg=colors[1],cex=size,type="o")
    points(LWP.mean~gxDate,col="black",data=subset(dat1,Treat=="wet"),pch=pchs[i],bg=colors[2],cex=size,type="o")
    
    
    #if (i==1) plotBy(LWP.mean~gxDate|Treat,data=dat1,type="b",xlim=c(startdate,enddate),legend=F,ylim=c(-10,0),pch=pchs[i],cex=size,add=T,bg=colors[i],
    #                 panel.first=adderrorbars(x=dat1$gxDate,y=dat1$LWP.mean,SE=dat1$LWP.standard.error,direction="updown"))
    #if (i > 1)plotBy(LWP.mean~gxDate|Treat,data=dat1,type="b",xlim=c(startdate,enddate),legend=F,ylim=c(-10,0),pch=pchs[i],cex=size,add=T,bg=colors,
    #                 panel.first=adderrorbars(x=dat1$gxDate,y=dat1$LWP.mean,SE=dat1$LWP.standard.error,direction="updown"))
    if (i%%2==1) magaxis(c(2,4),labels=c(1,1),frame.plot=TRUE,las=1,ylab="",tline=0.2,majorn=4,minorn=0)
    if (i<3)axis.Date(side=1,at=seq(startdate,enddate,by="month"),format="%m/%y",labels=FALSE)
    if (i>2)axis.Date(side=1,at=seq(startdate,enddate,by="month"),format="%m/%y")
    if (i==1) title(ylab="Pre-dawn",cex.lab=labsize,line=2.3,xpd=NA)
  }
  
  legend(x=as.Date("2012-10-1"),y=-7.1,legend=c("  Cacu","  Eusi","  Eute","  Pira"),pch=c(pchs),
         col=c(rep("black",4)),ncol=2,pt.bg="black",cex=1.4,bg="white")
  legend(x=as.Date("2012-10-8"),y=-7.1,legend=c("        ","        ","          ","           "),pch=c(pchs),
         col=c(rep("black",4)),ncol=2,pt.bg="white",cex=1.4,bty="n")
  #legend("bottomleft",c("Cacu-wet","Eusi-wet","Eute-wet","Pira-wet","Cacu-dry","Eusi-dry","Eute-dry","Pira-dry"),pch=c(pchs,pchs),
  #       col=c("black","black","black","black","darkgrey","darkgrey","darkgrey","darkgrey"),ncol=2,bg="white",cex=1.4)
  #plot mid-day
  plot.new()
  plot.window(xlim=c(startdate,enddate),ylim=c(-10,0))
  rect(xleft=dates[1],ybottom=-45,xright=dates[2],ytop=270,col=fillcol) #add rectangles for droughts
  rect(xleft=dates[3],ybottom=-45,xright=dates[4],ytop=270,col=fillcol) #add rectangles for droughts
  rect(xleft=dates[5],ybottom=-45,xright=enddate+13,ytop=270,col=fillcol) #add rectangles for droughts
  for (i in 1:length(lwp.m.l)){
    dat1 <- lwp.m.l[[i]]
    
    adderrorbars(x=dat1$gxDate,y=dat1$LWP.md.mean,SE=dat1$LWP.md.standard.error,direction="updown")
    points(LWP.md.mean~gxDate,col="black",data=subset(dat1,Treat=="dry"),pch=pchs[i],bg=colors[1],cex=size,type="o")
    points(LWP.md.mean~gxDate,col="black",data=subset(dat1,Treat=="wet"),pch=pchs[i],bg=colors[2],cex=size,type="o")
    
    
    #if (i==1) plotBy(LWP.md.mean~gxDate|Treat,data=dat1,type="b",xlim=c(startdate,enddate),legend=F,ylim=c(-10,0),pch=pchs[i],cex=size,add=T,col=colors,
    #                 panel.first=adderrorbars(x=dat1$gxDate,y=dat1$LWP.md.mean,SE=dat1$LWP.md.standard.error,direction="updown"))
    #if (i > 1)plotBy(LWP.md.mean~gxDate|Treat,data=dat1,type="b",xlim=c(startdate,enddate),legend=F,ylim=c(-10,0),pch=pchs[i],cex=size,add=T,col=colors,
    #                 panel.first=adderrorbars(x=dat1$gxDate,y=dat1$LWP.md.mean,SE=dat1$LWP.md.standard.error,direction="updown"))
    if (i%%2==1) magaxis(c(2,4),labels=c(1,1),frame.plot=TRUE,las=1,ylab="",tline=0.2,majorn=4,minorn=0,cex.lab=1.5)
    
    if (i<3)axis.Date(side=1,at=seq(startdate,enddate,by="month"),format="%m/%y",labels=FALSE)
    if (i>2)axis.Date(side=1,at=seq(startdate,enddate,by="month"),format="%m/%y")
    if (i==1) title(ylab="Mid-day",cex.lab=labsize,line=2.3,xpd=NA)
    
  }
  
  title(xlab="Date",outer=TRUE,ylab=expression(Leaf~water~potential~(MPa)),cex.lab=2,line=1)
  
  if(output==T) dev.copy2pdf(file="Output/LWP_over_time.pdf")
  
}
#---------------------------------------------------------------------------------------------------------------------






#---------------------------------------------------------------------------------------------------------------------
#- plot the dependence of g1 and non-stomatal limitation relative to VWC.
#   The fitting of NSL might need some work...
#---------------------------------------------------------------------------------------------------------------------
plotBetasG1NSL <- function(output=F,g1data,NSLdata,g1list,NSLlist){
  #- list for g1 data
  dat.l <- split(g1data,g1data$Species)
  
  #- fit NSL
  dat.l2 <- split(NSLdata,NSLdata$Species)
  
  #- plot normalized g1 and non-stomatal limitation as a function of VWC
  windows(16,16)
  par(mfrow=c(4,2),mar=c(0,0.25,0,0.25),xpd=FALSE,oma=c(4,5,1,5),cex=1.6,cex.axis=0.9,cex.lab=0.9)
  labs <- c("Cacu","Eusi","Eute","Pira")
  for (i in 1:length(dat.l)){
    newdat <- newdat2 <- data.frame()
    dat.temp <- dat.l[[i]]
    dat.temp$g1norm <- dat.temp$g1/max(dat.temp$g1)
    
    dat.temp$Species <- factor(dat.temp$Species)
    #------------------------------------------------------------------------
    #-- plot g1 vs. TDR
    
    plot(g1norm~TDR,data=subset(dat.temp,Treat=="wet"),pch=21,col="black",bg=grey(0.1),axes=F,ylim=c(0,1.05),xlim=c(0,0.4))
    points(g1norm~TDR,data=subset(dat.temp,Treat=="dry"),pch=21,col="black",bg=grey(0.8))
    magaxis(side=c(1:4),labels=c(0,1,0,0),las=1)
    if(i==4)  magaxis(side=c(1:4),labels=c(1,1,0,0),las=1)
    mtext(labs[i],side=2,xpd=T,cex=1.3,line=1.75)
    if(i==4)  mtext(expression(VWC~(m^3~m^-3)),side=1,outer=F,cex=1.5,line=2)
    if(i==1) legend("bottomright",xpd=NA,legend=c("Wet","Dry"),pch=21,pt.bg=c("black","grey"),ncol=2,cex=0.75)
    
    
    # plot model and SE from bootstrapping
    rm(newdat)
    newdat <- g1list[[2]][[i]]
    lines(wpred~TDR,data=newdat)
    polygon(x = c(newdat$TDR, rev(newdat$TDR)), y = c(newdat$lower, rev(newdat$upper)),
          col = alpha("grey",0.5), border = NA,xpd=F)
    

    #------------------------------------------------------------------------
    #-- repeat, but for NSL
    
    dat.temp2 <- dat.l2[[i]]
    dat.temp2$Species <- factor(dat.temp2$Species)
    
    plot(NSL~TDR,data=subset(dat.temp2,Treat=="wet"),pch=21,col="black",bg=grey(0.1),axes=F,ylim=c(0,1.5),xlim=c(0,0.4))
    points(NSL~TDR,data=subset(dat.temp2,Treat=="dry"),pch=21,col="black",bg=grey(0.8))
    magaxis(side=c(1:4),labels=c(0,0,0,1),las=1)
    if(i==4)  magaxis(side=c(1:4),labels=c(1,0,0,0),las=1)
    if(i==4)  mtext(expression(VWC~(m^3~m^-3)),side=1,outer=F,cex=1.5,line=2)
    
    newdat2 <- NSLlist[[2]][[i]]
    lines(wpred~TDR,data=newdat2)
    polygon(x = c(newdat2$TDR, rev(newdat2$TDR)), y = c(newdat2$lower, rev(newdat2$upper)), 
            col = alpha("grey",0.5), border = NA, xpd=F)
    lines(x=c(max(newdat2$TDR),max(dat.temp2$TDR)),y=c(1,1))
    
    
  }
  mtext(expression(normalized~g[1]),side=2,outer=T,cex=2.5,las=0,line=2.5)
  mtext(expression(Nonstomatal~limitation~(A/A[e])),side=4,outer=T,cex=2.5,las=0,line=2.5)
  
  if(output==T) dev.copy2pdf(file="Output/Beta_g1andNSL_VWC.pdf")
  
}
#---------------------------------------------------------------------------------------------------------------------










#---------------------------------------------------------------------------------------------------------------------
#- plot the dependence of g1 and non-stomatal limitation relative to LWP!
#---------------------------------------------------------------------------------------------------------------------
plotBetasG1NSL_LWP <- function(output=F,g1data,NSLdata,g1list,NSLlist){
  #- list for g1 data
  dat.l <- split(g1data,g1data$Species)
  
  #- fit NSL
  dat.l2 <- split(NSLdata,NSLdata$Species)
  
  #- plot normalized g1 and non-stomatal limitation as a function of LWP
  windows(16,16)
  par(mfrow=c(4,2),mar=c(0,0.25,0,0.25),xpd=FALSE,oma=c(4,5,1,5),cex=1.6,cex.axis=0.9,cex.lab=0.9)
  labs <- c("Cacu","Eusi","Eute","Pira")
  for (i in 1:length(dat.l)){
    newdat <- newdat2 <- data.frame()
    dat.temp <- dat.l[[i]]
    dat.temp$g1norm <- dat.temp$g1/max(dat.temp$g1)
    dat.temp$LWPpos <- dat.temp$LWP+11
    dat.temp$Species <- factor(dat.temp$Species)
    #------------------------------------------------------------------------
    #-- plot g1 vs. LWP
    
    plot(g1norm~LWP,data=subset(dat.temp,Treat=="wet"),pch=21,col="black",bg=grey(0.1),axes=F,ylim=c(0,1.05),xlim=c(-10,0))
    points(g1norm~LWP,data=subset(dat.temp,Treat=="dry"),pch=21,col="black",bg=grey(0.8))
    magaxis(side=c(1:4),labels=c(0,1,0,0),las=1)
    if(i==4)  magaxis(side=c(1:4),labels=c(1,1,0,0),las=1)
    mtext(labs[i],side=2,xpd=T,cex=1.3,line=1.75)
    if(i==4)  mtext(expression(psi[pd]~(mPa)),side=1,outer=F,cex=1.5,line=2)
    if(i==1) legend("topleft",xpd=NA,legend=c("Wet","Dry"),pch=21,pt.bg=c("black","grey"),ncol=1,cex=0.75)
    
    
    # plot model and SE from bootstrapping
    rm(newdat)
    newdat <- g1list[[2]][[i]]
    newdat$LWP <- newdat$LWPpos-11
    lines(wpred~LWP,data=newdat)
    polygon(x = c(newdat$LWP, rev(newdat$LWP)), y = c(newdat$lower, rev(newdat$upper)),
            col = alpha("grey",0.5), border = NA,xpd=F)
    
    
    #------------------------------------------------------------------------
    #-- repeat, but for NSL
    
    dat.temp2 <- dat.l2[[i]]
    dat.temp2$Species <- factor(dat.temp2$Species)
    
    plot(NSL~LWP,data=subset(dat.temp2,Treat=="wet"),pch=21,col="black",bg=grey(0.1),axes=F,ylim=c(0,1.5),xlim=c(-10,0))
    points(NSL~LWP,data=subset(dat.temp2,Treat=="dry"),pch=21,col="black",bg=grey(0.8))
    magaxis(side=c(1:4),labels=c(0,0,0,1),las=1)
    if(i==4)  magaxis(side=c(1:4),labels=c(1,0,0,0),las=1)
    if(i==4)  mtext(expression(psi[pd]~(mPa)),side=1,outer=F,cex=1.5,line=2)
    
    newdat2 <- NSLlist[[2]][[i]]
    newdat2$LWP <- newdat2$LWPpos-11
    
    lines(wpred~LWP,data=newdat2)
    polygon(x = c(newdat2$LWP, rev(newdat2$LWP)), y = c(newdat2$lower, rev(newdat2$upper)), 
            col = alpha("grey",0.5), border = NA, xpd=F)
    #lines(x=c(max(newdat2$TDR),max(dat.temp2$TDR)),y=c(1,1))
    
    
  }
  mtext(expression(normalized~g[1]),side=2,outer=T,cex=2.5,las=0,line=2.5)
  mtext(expression(Nonstomatal~limitation~(A/A[e])),side=4,outer=T,cex=2.5,las=0,line=2.5)
  
  if(output==T) dev.copy2pdf(file="Output/Beta_g1andNSL_LWP.pdf")
  
}
#---------------------------------------------------------------------------------------------------------------------












#-------------------------------------------------------------------------------------------------------------------------------------
#-- plot barplots of d13C and bigdelta 
#-------------------------------------------------------------------------------------------------------------------------------------
plotd13C <- function(export=F){
  
  #---- read in the data, do a little processing
  d1 <- get_d13C()
  
  #--------------------------------------------------------------
  # plots
  d1$bigDelta <- d1$bigDelta*1000
  d1.m <- summaryBy(deltaC+bigDelta~Species+Treat,data=d1,FUN=c(mean,standard.error))
  
  #-- barplots of d13C and bigdelta
  windows(18,12);par(mfrow=c(2,1),mar=c(3,5,1,1),oma=c(1,1,3,1))
  
  bp2 <- barplot(height=d1.m$deltaC.mean[1:8],col=c("darkgrey","white"),ylim=c(-35,-25),xpd=F,axes=F)
  adderrorbars(x=bp2,y=d1.m$deltaC.mean[1:8],SE=d1.m$deltaC.standard.error[1:8],direction="updown")
  magaxis(side=c(2,4),labels=c(1,0),frame.plot=T,las=1)
  text(x=c(1.25,3.7,6.1,8.5),y=17.5,xpd=T,labels=c("Cacu","Eusi","Eute","Pira"),cex=1.5)
  title(ylab=expression(paste(delta^{13}, "C (\u2030)")),cex.lab=1.5,line=1.8)
  text(x=c(1.25,3.7,6.1,8.5),y=-35.8,xpd=T,labels=c("Cacu","Eusi","Eute","Pira"),cex=1.5)
  legend("bottomright",xpd=NA,legend=c("Wet","Dry"),fill=c("darkgrey","white"),bty="n",ncol=1,cex=1.5)
  
  
  bp1 <- barplot(height=d1.m$bigDelta.mean[1:8],col=c("darkgrey","white"),ylim=c(18,25),xpd=F,axes=F)
  adderrorbars(x=bp1,y=d1.m$bigDelta.mean[1:8],SE=d1.m$bigDelta.standard.error[1:8],direction="updown")
  magaxis(side=c(2,4),labels=c(1,0),frame.plot=T,las=1)
  text(x=c(1.25,3.7,6.1,8.5),y=17.5,xpd=T,labels=c("Cacu","Eusi","Eute","Pira"),cex=1.5)
  title(ylab=expression(Delta~"*"~10^3),cex.lab=1.5,line=1.8)
  
  if(export==T) dev.copy2pdf(file="Output/ROS_d13C_bars.pdf")
}