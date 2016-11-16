
#------------------------------------------------------------------------------------------------------------------
#- make the soil moisture plot
#------------------------------------------------------------------------------------------------------------------
plotVWC <- function(ptsize=1.5,output=T,type="4panel"){
  
  #- get the VWC data from the loggers and handheld
  vwc.loggers <- get.TDR.logger(startdate="2011-11-01",enddate="2013-12-1")
  vwc.hand <- get.TDR.handheld(filterBlanks=1)
  vwc.hand$Species <- droplevels(vwc.hand$Species)
  
  #----
  # remove the handheld tdr measurements with <3 observations per treatment (some dry pots were mislabeled as wet)
  #writing my own length function that omits NA's
  length2 <- function(x,...){
    length(na.omit(x))
  }
  names(vwc.hand)[which(names(vwc.hand)=="VWC")] <- "TDR"
  vwc.hand0 <- summaryBy(TDR~Date+Species+Treat,FUN=c(mean,length2,standard.error),dat=vwc.hand,keep.names=FALSE,na.rm=TRUE)
  vwc.hand1 <- subset(vwc.hand0,TDR.length2>3 | Date>as.Date("2013-05-01") ) #remove dates that have <3 observations (some pots were miscoded)
  names(vwc.hand1)[4] <- "TDR"
  names(vwc.hand1)[6] <- "TDR.se"
  
    #-----
  # average across dates for dry pots only and have a look. Try to identify the dates where plants were re-watered
  vwc.dates.dry <- summaryBy(TDR~Date,data=subset(vwc.hand1,Treat=="dry"))
  

  names(vwc.loggers)[3] <- "Pot"
  vwc.loggers.out <- vwc.loggers[,c(6,5,4,3,1,2,7)]
  vwc.loggers.out$Treat <- as.factor(ifelse(vwc.loggers.out$Treat=="dry","Dry","Wet"))
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
      points(TDR~Date,col="black",data=subset(dat2,Treat=="dry"),pch=21,bg="white",cex=ptsize)
      points(TDR~Date,col="black",data=subset(dat2,Treat=="wet"),pch=21,bg="black",cex=ptsize)
      
      rug(gxdates,lwd=3,line=-0.5)
      
      #labeling and legends
      legend("topright",legend=dat1$Species[1],bty="n",cex=1.3)
      if (i<3)axis.Date(side=1,at=seq(startdate,enddate,by="month"),format="%m/%y",labels=FALSE)
      if (i>2)axis.Date(side=1,at=seq(startdate,enddate,by="month"),format="%m/%y")
      
      
    }
    title(xlab="Date",outer=TRUE,ylab=expression(Volumetric~water~content~(theta~m^3~m^-3)),cex.lab=2,line=2)
    title(xlab="Date",outer=TRUE,ylab=expression(Volumetric~water~content~(theta~m^3~m^-3)),cex.lab=2,line=2)
    if(output==T) dev.copy2pdf(file="Output/Figure1_VWC_fourPanels.pdf")
    
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
      
      adderrorbars(x=dat2$Date,y=dat2$TDR,SE=dat2$TDR.se,direction="updown")
      
      points(TDR~Date,col="black",data=subset(dat2,Treat=="dry"),pch=pchs[i],bg="white",cex=ptsize)
      points(TDR~Date,col="black",data=subset(dat2,Treat=="wet"),pch=pchs[i],bg="black",cex=ptsize)
      
      
    }
    magaxis(c(2,4),labels=c(1,1),las=1,cex.axis=1.5,frame.plot=TRUE)
    #rug(gxdates,lwd=3,line=-0.5)
    axis.Date(side=1,at=seq(startdate,enddate,by="month"),format="%m/%y",cex.axis=1.5,tcl=0.5)
    abline(h=0.1,lty=3)
    title(xlab="Date",outer=TRUE,ylab=expression(Volumetric~water~content~(theta~";"~m^3~m^-3)),cex.lab=2,line=-1)
    legend(x=as.Date("2013-1-1"),y=0.416,legend=c("  Cacu","  Eusi","  Eute","  Pira"),pch=c(pchs),
           col=c(rep("black",4)),ncol=2,pt.bg="black",cex=1.4,bg="white")
    legend(x=as.Date("2013-1-8"),y=0.416,legend=c("        ","        ","          ","           "),pch=c(pchs),
           col=c(rep("black",4)),ncol=2,pt.bg="white",cex=1.4,bty="n")
    rug(gxdates,lwd=3,line=-0.5)
    
  }
  if(output==T) dev.copy2pdf(file="./Output/Figure1_VWC_onePanel.pdf")
  
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
  
  # names(lwp)[3] <- "LWPdate"
  # lwp.pd <- subset(lwp,Type=="PD" & is.na(LWP)==FALSE)
  # lwp.pd$LWP <- -1*lwp.pd$LWP
  # lwp.md <- subset(lwp,Type=="MD" & is.na(LWP)==FALSE)
  # lwp.md$LWP <- -1*lwp.md$LWP
  # 
  # names(lwp.md)[2] <- "LWP.md"
  # lwp2 <- merge(lwp.pd,lwp.md,by=c("Pot","LWPdate","Treat","Species"),all=FALSE)
  # lwp2$diff <- with(lwp2,abs(LWP.md)-abs(LWP))
  # 
  lwp.trt<- summaryBy(LWP.pd+LWP.md+diff~LWPdate+Species+Treat,data=lwp,FUN=c(mean,standard.error), na.rm=TRUE)
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
                 expression(WUE[i]),expression(C[i]~"/"~C[a]))
  
  #-- get the average WUE for the first two droughts
  firstdrs <- subset(ITE.trt,gxDate==as.Date("2012-10-31") | gxDate==as.Date("2012-11-6") | gxDate==as.Date("2012-12-19"))
  summaryBy(WUE.mean~Species+Treat,data=firstdrs)
  
  windows(14,12)
  par(mfrow=c(4,4),oma=c(8,7,2,5),mar=c(0.25,0.25,0.25,0.25))
  count <- 0
  for (i in 1:4){
    #print(i)
    #plot photosynthesis
    for (j in 1:4){
      #print(j)
      count <- count +1
      plot.gs.ros(dat=dat.trt.list[[j]],toplot=toplot[i],toplotse=ses[i],
                  ylims=ylims.list[[i]],ylabs=ylabs.list[[j]],xlab=xlabs[i])
      legend("topright",letters[count],bty="n",cex=1.5,inset=-0.02)
    }
    title(xlab="",outer=TRUE,ylab=labels[[i]],cex.lab=2,line=4,adj=((4-i)/4+0.13))
    if (i==1) title(main=expression(Cacu~~~~~~~~~~~~~~~~~~~~~~Eusi~~~~~~~~~~~~~~~~~~~~~~Eute~~~~~~~~~~~~~~~~~~~~~Pira),outer=T,cex.main=2)
  }
  title(xlab="Date",outer=TRUE,ylab="",cex.lab=2,line=5)
  if(output==T) dev.copy2pdf(file="Output/Figure3_Asat_Gs_WUE_CiCa_ROS.pdf")
}
#---------------------------------------------------------------------------------------------------------------------






#------------------------------------------------------------------------------------------------------------------
#- plot leaf water potentials over time
#------------------------------------------------------------------------------------------------------------------
plotLWP <- function(fillcol="lightgrey",size=1.75,output=F,labsize=1.8){
  
  ros2 <- return.gx.vwc.lwp()
  
  
  
  
  #---------------------------------------------------------------------------------------------------------------------
  #--- plot pre-dawn and mid-day leaf water potentials over time
  lwp.m <- summaryBy(LWP.pd+LWP.md+diff~gxDate+Treat+Species,FUN=c(mean,standard.error),na.rm=T,data=subset(ros2,gxDate>as.Date("2012-11-1")))
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
    
    adderrorbars(x=dat1$gxDate,y=dat1$LWP.pd.mean,SE=dat1$LWP.pd.standard.error,direction="updown")
    points(LWP.pd.mean~gxDate,col="black",data=subset(dat1,Treat=="dry"),pch=pchs[i],bg=colors[1],cex=size,type="o")
    points(LWP.pd.mean~gxDate,col="black",data=subset(dat1,Treat=="wet"),pch=pchs[i],bg=colors[2],cex=size,type="o")
    
    
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
  legend("topleft","a",bty="n",inset=-0.05,cex=1.5)
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
  legend("topleft","b",bty="n",inset=-0.05,cex=1.5)
  
  title(xlab="Date",outer=TRUE,ylab=expression(Leaf~water~potential~(MPa)),cex.lab=2,line=1)
  
  if(output==T) dev.copy2pdf(file="Output/Figure2_LWP_over_time.pdf")
  
}
#---------------------------------------------------------------------------------------------------------------------






#---------------------------------------------------------------------------------------------------------------------
#- plot the dependence of g1 and non-stomatal limitation relative to VWC.
#   The fitting of NSL might need some work...
#---------------------------------------------------------------------------------------------------------------------
plotBetasG1NSL <- function(output=F,g1data,NSLdata,g1list,NSLlist,xlims=c(0,0.33)){
  #- list for g1 data
  dat.l <- split(g1data,g1data$Species)
  
  #- fit NSL
  dat.l2 <- split(NSLdata,NSLdata$Species)
  
  #- plot normalized g1 and non-stomatal limitation as a function of VWC
  windows(16,16)
  par(mfrow=c(4,2),mar=c(0,0.25,0,0.25),xpd=FALSE,oma=c(4,5,1,5),cex=1.6,cex.axis=0.9,cex.lab=0.9)
  labs <- c("Cacu","Eusi","Eute","Pira")
  count <- 0
  for (i in 1:length(dat.l)){
    newdat <- newdat2 <- data.frame()
    dat.temp <- dat.l[[i]]
    dat.temp$g1norm <- dat.temp$g1/max(dat.temp$g1)
    
    dat.temp$Species <- factor(dat.temp$Species)
    #------------------------------------------------------------------------
    #-- plot g1 vs. TDR
    count <- count+1
    
    plot(g1norm~TDR,data=subset(dat.temp,Treat=="wet"),pch=21,col="black",bg=grey(0.1),axes=F,ylim=c(0,1.18),xlim=xlims)
    points(g1norm~TDR,data=subset(dat.temp,Treat=="dry"),pch=21,col="black",bg=grey(0.8))
    magaxis(side=c(1:4),labels=c(0,1,0,0),las=1,tcl=0.3,ratio=0.25,majorn=3,cex.axis=0.8)
    if(i==4)  magaxis(side=c(1:4),labels=c(1,1,0,0),las=1,tcl=0.3,ratio=0.25,majorn=3,cex.axis=0.8)
    mtext(labs[i],side=2,xpd=T,cex=1.3,line=1.75)
    if(i==1) legend("bottomright",xpd=NA,legend=c("Wet","Dry"),pch=21,pt.bg=c("black","grey"),ncol=2,cex=0.75)
    if(i==1) title(main="Stomatal",cex.main=0.75,line=0.5,xpd=NA)
    legend("topleft",letters[count],cex=1,inset=-0.1,bty="n")
    
    
    # plot model and SE from bootstrapping
    rm(newdat)
    newdat <- g1list[[2]][[i]]
    lines(wpred~TDR,data=newdat,xpd=F)
    par(xpd=F)
    polygon(x = c(newdat$TDR, rev(newdat$TDR)), y = c(newdat$lower, rev(newdat$upper)),
          col = alpha("grey",0.5),border=NA)
    

    #------------------------------------------------------------------------
    #-- repeat, but for NSL
    
    dat.temp2 <- dat.l2[[i]]
    dat.temp2$Species <- factor(dat.temp2$Species)
    count <- count+1
    
    plot(NSL~TDR,data=subset(dat.temp2,Treat=="wet"),pch=21,col="black",bg=grey(0.1),axes=F,ylim=c(0,1.18),xlim=xlims)
    points(NSL~TDR,data=subset(dat.temp2,Treat=="dry"),pch=21,col="black",bg=grey(0.8))
    magaxis(side=c(1:4),labels=c(0,0,0,1),las=1,tcl=0.3,ratio=0.25,majorn=3,cex.axis=0.8)
    if(i==4)  magaxis(side=c(1:4),labels=c(1,0,0,0),las=1,tcl=0.3,ratio=0.25,majorn=3,cex.axis=0.8)
    if(i==1) title(main="Non-stomatal",cex.main=0.75,line=0.5,xpd=NA)
    
    legend("topleft",letters[count],cex=1,inset=-0.1,bty="n")
    
    newdat2 <- NSLlist[[2]][[i]]
    lines(wpred~TDR,data=newdat2)
    polygon(x = c(newdat2$TDR, rev(newdat2$TDR)), y = c(newdat2$lower, rev(newdat2$upper)), 
            col = alpha("grey",0.5),border=NA)
    lines(x=c(max(newdat2$TDR),max(dat.temp2$TDR)),y=c(1,1))
    
    
  }
  mtext(expression(Normalized~g[1]),side=2,outer=T,cex=2,las=0,line=3)
  mtext(expression(A/A[e]),side=4,outer=T,cex=2,las=0,line=2.5)
  mtext(expression(theta~(m^3~m^-3)),side=1,outer=T,cex=1.5,line=2)
  
  if(output==T) dev.copy2pdf(file="Output/Figure4_Beta_g1andNSL_VWC.pdf")
  
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
  count <- 0
  #- plot normalized g1 and non-stomatal limitation as a function of LWP
  windows(16,16)
  par(mfrow=c(4,2),mar=c(0,0.25,0,0.25),xpd=FALSE,oma=c(4,5,1,5),cex=1.6,cex.axis=0.9,cex.lab=0.9)
  labs <- c("Cacu","Eusi","Eute","Pira")
  for (i in 1:length(dat.l)){
    newdat <- newdat2 <- data.frame()
    dat.temp <- dat.l[[i]]
    dat.temp$g1norm <- dat.temp$g1/max(dat.temp$g1)
    dat.temp$LWPpos <- dat.temp$LWP.pd+11
    dat.temp$Species <- factor(dat.temp$Species)
    #------------------------------------------------------------------------
    #-- plot g1 vs. LWP    
    count <- count+1

    
    plot(g1norm~LWP.pd,data=subset(dat.temp,Treat=="wet"),pch=21,col="black",bg=grey(0.1),axes=F,ylim=c(0,1.05),xlim=c(-10,0))
    points(g1norm~LWP.pd,data=subset(dat.temp,Treat=="dry"),pch=21,col="black",bg=grey(0.8))
    magaxis(side=c(1:4),labels=c(0,1,0,0),las=1,tcl=0.3,ratio=0.25,majorn=3,cex.axis=0.8)
    if(i==4)  magaxis(side=c(1:4),labels=c(1,1,0,0),las=1,tcl=0.3,ratio=0.25,majorn=3,cex.axis=0.8)
    mtext(labs[i],side=2,xpd=T,cex=1.3,line=1.75)
    if(i==4)  mtext(expression(psi[pd]~(MPa)),side=1,outer=F,cex=1.5,line=2)
    if(i==1) legend("topleft",xpd=NA,legend=c("Wet","Dry"),pch=21,pt.bg=c("black","grey"),ncol=1,cex=0.75)
    legend("bottomleft",letters[count],cex=1,inset=-0.05,bty="n")
    if(i==1) title(main="Stomatal",cex.main=0.75,line=0.5,xpd=NA)
    
    
    # plot model and SE from bootstrapping
    rm(newdat)
    newdat <- g1list[[2]][[i]]
    #newdat$LWP <- newdat$LWPpos-11
    lines(wpred~LWP,data=newdat)
    polygon(x = c(newdat$LWP, rev(newdat$LWP)), y = c(newdat$lower, rev(newdat$upper)),
            col = alpha("grey",0.5), border = NA,xpd=F)
    
    
    #------------------------------------------------------------------------
    #-- repeat, but for NSL
    count <- count+1
    
    dat.temp2 <- dat.l2[[i]]
    dat.temp2$Species <- factor(dat.temp2$Species)
    
    plot(NSL~LWP.pd,data=subset(dat.temp2,Treat=="wet"),pch=21,col="black",bg=grey(0.1),axes=F,ylim=c(0,1.09),xlim=c(-10,0))
    points(NSL~LWP.pd,data=subset(dat.temp2,Treat=="dry"),pch=21,col="black",bg=grey(0.8))
    magaxis(side=c(1:4),labels=c(0,0,0,1),las=1,tcl=0.3,ratio=0.25,majorn=3,cex.axis=0.8)
    if(i==4)  magaxis(side=c(1:4),labels=c(1,0,0,0),las=1,tcl=0.3,ratio=0.25,majorn=3,cex.axis=0.8)
    if(i==4)  mtext(expression(psi[pd]~(MPa)),side=1,outer=F,cex=1.5,line=2)
    legend("bottomleft",letters[count],cex=1,inset=-0.05,bty="n")
    if(i==1) title(main="Non-stomatal",cex.main=0.75,line=0.5,xpd=NA)
    
    newdat2 <- NSLlist[[2]][[i]]
    #newdat2$LWP <- newdat2$LWPpos-11
    
    if (count %in% c(7,8)){
      newdat2 <- subset(newdat2,LWP> -3)
    }
    lines(wpred~LWP,data=newdat2)
    polygon(x = c(newdat2$LWP, rev(newdat2$LWP)), y = c(newdat2$lower, rev(newdat2$upper)), 
            col = alpha("grey",0.5), border = NA, xpd=F)
    #lines(x=c(max(newdat2$TDR),max(dat.temp2$TDR)),y=c(1,1))
    
    
  }
  mtext(expression(normalized~g[1]),side=2,outer=T,cex=2.5,las=0,line=2.5)
  mtext(expression(A/A[e]),side=4,outer=T,cex=2.5,las=0,line=2.5)
  
  if(output==T) dev.copy2pdf(file="Output/FigureS4_Beta_g1andNSL_LWP.pdf")
  
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
  legend("topright",xpd=NA,legend=c("Wet","Dry"),fill=c("darkgrey","white"),bty="n",ncol=1,cex=1.5)
  
  if(export==T) dev.copy2pdf(file="Output/ROS_d13C_bars.pdf")
}








#-------------------------------------------------------------------------------------------------------------------------------------
#- make publication-ready plot of soil moisture release curve
#-------------------------------------------------------------------------------------------------------------------------------------
plotMoistCurve <- function(output=F){
  
  
  
  #-----
  #- get the soil data
  curve <- getMoistCurve()
  curve$pressure_MPa_neg <- -1*curve$pressure_MPa
  #-----
  
  
  
  #-----
  #- get the leaf water potential data
  lwp1 <- return.gx.vwc.lwp()
  lwp <- subset(lwp1,!(gxDate %in% as.Date(c("2012-10-04","2012-10-31"))))
  lwp.mean <- summaryBy(LWP.pd+TDR.mean~gxDate+Species,data=subset(lwp,Treat=="dry"),FUN=c(mean,standard.error))
  #-----
  
  
  
  #----
  #- predict soil matrix potential via Campbell 1974 and Cosby et al. 1984
  thetaSat <- 0.09
  #- values  Duursma 2008, loamy sand 
  b <- 4.26
  PsiE <- -0.36
  Ksat <- 5.8
  
  VWC <- seq(0,0.3,length=101)
  PsiSoil <- PsiE*(VWC/thetaSat)^(-1*b)
  
  Cosby_optim <- function(par,VWC,dat){
    #- get the parameters
    thetaSat <- par[1]
    b <- par[2]
    PsiE <- par[3]

    #- predict matrix potential
    pred <- PsiE*(VWC/thetaSat)^(-1*b)
    
    #- 
    cost <- ((dat-pred)^2)/abs(dat)
    return(sum(cost))
  }
  
  set.seed(1234)
  upper <- c(0.2,15,-0.01)
  lower <- c(0.1,1,-10)
  
  #--- fit and predict off of the soil data
  DEoutput <- DEoptim(Cosby_optim,lower=lower,upper=upper,VWC=curve$VWC,dat=curve$pressure_MPa_neg,
                    DEoptim.control(NP=400,itermax=200,trace=F))
  bestpars <- DEoutput$optim$bestmem
  names(bestpars) <- c("thetaSat","b","PsiE")
  
  VWC <- seq(0,0.3,length=1001)
  PsiSoil_best <- bestpars[3]*(VWC/bestpars[1])^(-1*bestpars[2])
  #---
  
  
  
  #--- Try to fit Cosby in a log-transformed framework (See equation 10 in Duursma 2008)
  thetasat <- 0.14
  curve$Xval <- log10(curve$VWC/thetasat)
  curve$Yval <- log10(curve$pressure_MPa)
  Cosby_loglin <- lm(Yval~Xval,data=subset(curve,VWC<0.1))
  plot(Yval~Xval,dat=curve)
  abline(Cosby_loglin)
  b <- coef(Cosby_loglin)[[2]]
  PsiE <- -10^coef(Cosby_loglin)[[1]]
  PsiSoil_loglin <- PsiE*(VWC/thetasat)^(b)
  
  
  
  #--- fit the van Genuchten model via the package "HydroMe"
  vn.ns=nlsLM(VWC~SSvgm(actual_pressure,thr,ths,alp,nscal,mscal),data=curve,trace=T,
              control=nls.lm.control(maxiter=300,options(warn=1,trace=T)))
  Psi <- seq(from=0,to=1000,length.out=1001)
  newdata <- data.frame(actual_pressure=Psi)
  newdata$pred <- predict(vn.ns,newdata=newdata)
  newdata$psi <- newdata$actual_pressure/-10
  #---
  
  
  #---
  #- fit and predict off of the LWP data
  dat2 <- subset(lwp.m,Species != "Pira")
  set.seed(1234)
  upper2 <- c(0.3,15,-0.01)
  lower2 <- c(0.001,0.1,-10)
  
  DEoutput2 <- DEoptim(Cosby_optim,lower=lower2,upper=upper2,VWC=dat2$TDR,dat=dat2$LWP.pd,
                      DEoptim.control(NP=400,itermax=200,trace=T))
  bestpars2 <- DEoutput2$optim$bestmem
  names(bestpars2) <- c("thetaSat","b","PsiE")
  PsiSoil_best2 <- bestpars2[3]*(VWC/bestpars2[1])^(-1*bestpars2[2])
  #---
  
  
  #- plot
  windows(12,12);par(cex.lab=1.5,cex.axis=1.5,mar=c(5,5,1,1))
  plot(pressure_MPa_neg~VWC,data=curve,axes=F,pch=1,cex=2,ylim=c(-10,0),xlim=c(0,0.2),
       ylab=expression(Pre-dawn~leaf~or~soil~water~potential~(Psi~";"~MPa)),xlab=expression(Volumetric~water~content~(theta~";"~m^3~m^-3)))
  plotBy(LWP.pd.mean~TDR.mean.mean|Species,data=lwp.mean,legend=F,add=T,pch=16,cex=2,
         panel.first=adderrorbars(x=lwp.mean$TDR.mean.mean,y=lwp.mean$LWP.pd.mean,
                                  SE=lwp.mean$LWP.pd.standard.error,direction="updown"))
  
  magaxis(c(1:4),labels=c(1,1,0,0),las=1)
  
  #lines(PsiSoil_loglin~VWC,lty=1) # overlay fit of log-linearized soil data
  lines(PsiSoil_best~VWC,lty=1)
  #lines(PsiSoil_best2~VWC,lty=2) #overlay fit of the pre-dawn data, excluding Pira?
  lines(psi~pred,data=newdata,lty=2)
  abline(h=0,lty=3)
  
  #lines(PsiSoil~VWC,lty=2)
  
  legend(x=0.125,y=-7,lty=c(1,2),col="black",legend=c("Cosby","van Genutchen"),bty="n",cex=1.2)
  legend(x=0.137,y=-7.8,pch=c(1,16,16,16,16),col=c("black","black","red","green3","blue")
         ,legend=c("Soil measurement","Cacu","Eusi","Eute","Pira"),bty="n",cex=1.2)
  
  
  if(output==T) dev.copy2pdf(file="Output/FigureS5_LWP_moistureReleaseCurve.pdf")
  
  return(bestpars)
}
#-------------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------------
#- plot Asat and gs relative to pre-dawn leaf water potential for the four species. Isohydry vs. anisohydry
#-------------------------------------------------------------------------------------------------------------------------------------
plotHydry <- function(output=F){
  dat <- return.gx.vwc.lwp()
  
  #- remove a troublesome Cacu point
  dat[which(dat$Photo>30 & dat$Species=="cacu" & dat$LWP< -2),c("Photo")] <- NA
  #dat.m <- summaryBy(Photo+Cond+LWP+LWP.md ~ Species+Treat+Date,data=dat,FUN=c(mean,standard.error))
  dat.l <- split(dat,dat$Species)
  
  windows(16,16)
  par(mfrow=c(4,2),mar=c(0,0.25,0,0.25),xpd=FALSE,oma=c(4,5,1,5),cex=1.6,cex.axis=0.9,cex.lab=0.9)
  
  count <- 0
  for (i in 1:length(dat.l)){
    toplot <- dat.l[[i]]
    
    count <- count +1
    plot(Photo~LWP,data=toplot,ylim=c(-2,35),xlim=c(-10,0),axes=F)
    magaxis(side=c(1,2,3,4),labels=c(0,1,0,0),las=1,frame.plot=T,tcl=0.2)
    title(ylab=toplot$Species[1],cex.lab=1,xpd=NA,line=1.5)
    if (i ==4) magaxis(side=c(1,2,3,4),labels=c(1,0,0,0),las=1,frame.plot=T,tcl=0.2)
    abline(h=0,lty=2)
    legend("topleft",letters[count],bty="n",inset=-0.1)
    
    count <- count +1
    plot(Cond~LWP,data=toplot,ylim=c(-0.05,1),xlim=c(-10,0),axes=F)
    magaxis(side=c(1,2,3,4),labels=c(0,0,0,1),las=1,frame.plot=T,tcl=0.2)
    if (i ==4) magaxis(side=c(1,2,3,4),labels=c(1,0,0,0),las=1,frame.plot=T,tcl=0.2)
    abline(h=0,lty=2)
    legend("topleft",letters[count],bty="n",inset=-0.1)
    
    
  }  
  title(xlab=expression(psi[l-PD]~(MPa)),outer=TRUE,cex.lab=1.5,line=2.5)
  title(ylab=expression(A[sat]~(mu*mol~m^-2~s^-1)),outer=TRUE,cex.lab=1.5,line=2.5)
  title(ylab=expression(g[s]~(mol~m^-2~s^-1)),outer=TRUE,cex.lab=1.5,line=-18,xpd=NA)
  if (output==T) dev.copy2pdf(file="Output/FigureS3_isohydry.pdf")
}

#-------------------------------------------------------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------------------------------------------------------
#- compare the isotopes and the long-term gas exchange data
#-------------------------------------------------------------------------------------------------------------------------------------
plotd13C_gx <- function(output=F,ptsize=1.8){
  #---- read in the isotope data, do a little processing
  d1 <- get_d13C()
  
  d1$bigDelta <- d1$bigDelta*1000
  d1.m <- summaryBy(deltaC+bigDelta~Species+Treat,data=d1,FUN=c(mean,standard.error))
  d1.m$Species <- tolower(d1.m$Species)
  
  #--- read in the leaf gas exchange data, do a little processing
  gx <- return.gx.vwc()
  gx$WUEi <- with(gx,Photo/Cond)
  
  #- try to merge individual pots?
  gx.m2 <- summaryBy(Photo+Cond+WUEi~Species+Treat+Pot,data=subset(gx,gxDate < as.Date("2013-3-01")),FUN=mean,keep.names=T)
  dat2 <- merge(gx.m2,subset(d1,Date==as.Date("2013-03-26")),by=c("Species","Treat","Pot"))
  
  
  #- plot isotopes relative to gas exchange,
  #  overlay fits
  windows(18,12);par(mfrow=c(1,2),mar=c(5,5,1,1),oma=c(1,1,3,1))
  
  lm1 <- list()
  dat.l <- split(dat2,dat2$Species)
  plot(bigDelta~WUEi,data=subset(dat2,Treat=="wet" & Species=="cacu"),pch=16,ylim=c(17,27),xlim=c(40,120),
       axes=F,ylab="",xlab="",col="white")
  pointstyles <- 21:24
  #pointstyles <- c(15,16,17,18)
  #pointstyles2 <- c(0,1,2,5)
  for (i in 1:length(dat.l)){
    
    lm1[[i]] <- lm(bigDelta~WUEi,data=dat.l[[i]])
    predline(lm1[[i]])
  }
  for (i in 1:length(dat.l)){
    
    points(bigDelta~WUEi,data=subset(dat.l[[i]],Treat=="wet"),col="black",pch=pointstyles[i],bg="black",cex=ptsize)
    points(bigDelta~WUEi,data=subset(dat.l[[i]],Treat=="dry"),col="black",pch=pointstyles[i],bg="white",cex=ptsize)
    
  }
  #legend("topright",pch=15:18,legend=levels(dat2$Species),cex=1.2)
  magaxis(side=1:4,labels=c(1,1,0,0),frame.plot=T,las=1)
  title(ylab=expression(Delta~"*"~10^3),cex.lab=2)
  title(xlab=expression(WUE[i]~(A[sat]/g[s]~";"~mu*mol~mol^-1)),cex.lab=2)
  legend("topleft","a",bty="n",cex=1.5)
  legend(x=75,y=27,legend=c("  Cacu","  Eusi","  Eute","  Pira"),pch=c(pointstyles),
         col=c(rep("black",4)),ncol=2,pt.bg="black",cex=1.4,bg="white")
  legend(x=78,y=27,legend=c("        ","        ","          ","           "),pch=c(pointstyles),
         col=c(rep("black",4)),ncol=2,pt.bg="white",cex=1.4,bty="n")
  
  #- fit one big model instead
  lm.all <- lm(bigDelta~WUEi+Species+WUEi:Species,data=dat2)
  anova(lm.all)
  summary(lm.all)
  
  #-- barplot of bigdelta
  d1.m <- summaryBy(deltaC+bigDelta~Species+Treat,data=d1,FUN=c(mean,standard.error))
  
  bp1 <- barplot(height=d1.m$bigDelta.mean[1:8],col=c("darkgrey","white"),ylim=c(18,25),xpd=F,axes=F)
  adderrorbars(x=bp1,y=d1.m$bigDelta.mean[1:8],SE=d1.m$bigDelta.standard.error[1:8],direction="updown")
  magaxis(side=c(2,4),labels=c(1,0),frame.plot=T,las=1)
  text(x=c(1.25,3.7,6.1,8.5),y=17.5,xpd=T,labels=c("Cacu","Eusi","Eute","Pira"),cex=1.5)
  title(ylab=expression(Delta~"*"~10^3),cex.lab=1.5,line=1.8)
  legend("topright",xpd=NA,legend=c("Wet","Dry"),fill=c("darkgrey","white"),bty="n",ncol=1,cex=1.5)
  legend("topleft","b",bty="n",cex=1.5)
  title(xlab="Species",cex.lab=2,line=3.5)
  if(output==T) dev.copy2pdf(file="Output/Figure6_d13C_WUE.pdf")
}
#-------------------------------------------------------------------------------------------------------------------------------------








#-------------------------------------------------------------------------------------------------------
#- Function to model Asat's response to drought given changes in g1, photosynthetic capacity, or both
#-------------------------------------------------------------------------------------------------------
modelAsatVWC <- function(output=F,fit.spg1,fit.spNSL){
  
  #- maximum Vcmax data (see returnVcmaxa(), but get the 75th (95th?) percentile instead)
  maxVcmax <- data.frame(Species=c("cacu","eusi","eute","pira"),Vcmax_max=c(82,84,79,81))
  
  # #- get the g1 fits, merge in the max for each species
  g1values <- returng1()
  maxG1 <- summaryBy(g1~Species,data=g1values,FUN=max,keep.names=T)
  
  #- take an alternative approach, and predict values across a wide range of VWC values, along with Tair of 25, VPD of 1.5
  predData.1 <- expand.grid(VWC=seq(0.001,0.4,length=101),Species=levels(maxVcmax$Species))
  predData.2 <- merge(predData.1,maxVcmax,by="Species")
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
    plotBy(ALEAF_NSL~VWC,data=toplot,type="l",lwd=2,lty=2,add=F,ylim=c(-1,25),legend=F,axes=F,col="black")
    plotBy(ALEAF_g1~VWC,data=toplot,type="l",lwd=2,lty=3,add=T,legend=F,col="black")
    plotBy(ALEAF_both~VWC,data=toplot,type="l",lwd=2,lty=1,add=T,legend=F,col="black")
    points(Photo~VWC,data=subset(Adat.m,Species==as.character(toplot$Species[1])),cex=1.5,col="black",pch=16)
    
    title(main=toplot$Species[1],line=-1.5)
    
    if(i==1) magaxis(side=c(1:4),labels=c(0,1,0,0),frame.plot=T,las=1)
    if(i==2) magaxis(side=c(1:4),labels=c(0,0,0,2),frame.plot=T,las=1)
    if(i==3) magaxis(side=c(1:4),labels=c(1,1,0,0),frame.plot=T,las=1)
    if(i==4) magaxis(side=c(1:4),labels=c(1,0,0,1),frame.plot=T,las=1)
    legend("topleft",letters[i],bty="n",cex=1.5,inset=0)
    
    if (i==1) legend("bottomright",legend=c("Stomatal","Non-stomatal","Both"),lty=c(3,2,1),cex=1.5,bty="n",lwd=2)
  }
  title(ylab=expression(A[sat]~(mu*mol~m^-2~s^-1)),outer=T,cex.lab=2)
  title(xlab=expression(theta~(m^3~m^-3)),outer=T,cex.lab=2)
  if (output==T) dev.copy2pdf(file="Output/Figure5_Asat_VWC_modelPredictions.pdf")
}
#-------------------------------------------------------------------------------------------------------





#------------------------------------------------------------------------------------------------------------------
#- A new function that plots A, gs, WUE, Ci/Ca, and LWP-PD an LWP-MD as a function of volumetric water content
#------------------------------------------------------------------------------------------------------------------
plotGX_theta <- function(output=F,xlims = c(0,0.33),colors= brewer.pal(4,"Set1")){
  #- get the gas exchange and handheld TDR data
  ros <- return.gx.vwc()
  
  
  #- get the lwp data and process a bit
  #lwp <- read.lwp()
  lwp1 <- return.gx.vwc.lwp()
  lwp <- subset(lwp1,!(gxDate %in% as.Date(c("2012-10-04","2012-10-31"))))
  
  # names(lwp)[3] <- "LWPdate"
  # lwp.pd <- subset(lwp,Type=="PD" & is.na(LWP)==FALSE)
  # lwp.pd$LWP <- -1*lwp.pd$LWP
  # lwp.md <- subset(lwp,Type=="MD" & is.na(LWP)==FALSE)
  # lwp.md$LWP <- -1*lwp.md$LWP
  # 
  # names(lwp.md)[2] <- "LWP.md"
  # lwp2 <- merge(lwp.pd,lwp.md,by=c("Pot","LWPdate","Treat","Species"),all=FALSE)
  # lwp2$diff <- with(lwp2,abs(LWP.md)-abs(LWP))
  
  #- log-transform LWP number
  lwp$LWP.pd.trans <- -log(-1*lwp$LWP.pd)
  lwp$LWP.md.trans <- -log(-1*lwp$LWP.md)
  
  lwp.trt<- summaryBy(LWP.pd+LWP.md+LWP.pd.trans+LWP.pd.trans+diff+TDR~LWPdate+Species+Treat,data=lwp,FUN=c(mean,standard.error), na.rm=TRUE)
  lwp.trt$TDR <- lwp.trt$TDR.mean
  lwp.trt.list <- split(lwp.trt,lwp.trt$Species)
  #------------------------------------------------------------------------------------------------------------------
  
  
  
  #do some data manipulation, set graphing parameters 
  ros$ITE <- with(ros,Photo/Trmmol)
  ros$WUE <- with(ros,Photo/Cond)
  
  ITE.trt <- summaryBy(Photo+Cond+ITE+WUE+TDR+Ci.Ca~gxDate+Species+Treat,data=ros,FUN=c(mean,standard.error), na.rm=TRUE)
  ITE.trt$TDR <- ITE.trt$TDR.mean
  ITE.trt.list <- split(ITE.trt,ITE.trt$Species)
  
  
  
  windows(25,14)
  par(mfrow=c(2,3),oma=c(8,7,2,5),mar=c(0.25,6.5,0.25,0.5))
  symbols <- c(16,17,18,15)
  colors <- colors
  
  
  
  #- plot Photo
  for(i in 1:4){
    toplot <- ITE.trt.list[[i]]
    if(i ==1) plot.new();plot.window(xlim=xlims,ylim=c(0,25));abline(h=0)
    
    
    #- smoothplot
    smoothplot(TDR, Photo.mean, polycolor=alpha(colors[i],0.3),linecols=colors[i],
               ylab="",cex=0.5,add=T,
               data=toplot, kgam=5, axes=F)
    #- overlay points and error bars
    plotBy(Photo.mean~TDR|Treat,data=toplot,add=T,pch=symbols[i],col=c(colors[i]),cex=2,legend=F,
           panel.first=adderrorbars(x=toplot$TDR,y=toplot$Photo.mean,SE=toplot$Photo.standard.error,direction="updown"))
    
  }
  magaxis(side=1:4,labels=c(0,1,0,0),las=2,cex.axis=1.5)
  title(ylab=expression(A[sat]~(mu*mol~m^-2~s^-1)),xpd=NA,cex.lab=3)
  legend("bottomright",letters[1],bty="n",cex=2)
  
  #- plot WUE
  for(i in 1:4){
    toplot <- ITE.trt.list[[i]]
    if(i ==1) plot.new();plot.window(xlim=xlims,ylim=c(-30,250));abline(h=0)
    
    #- smoothplot
    smoothplot(TDR, WUE.mean, polycolor=alpha(colors[i],0.3),linecols=colors[i],
               ylab="",cex=0.5,add=T,
               data=toplot, kgam=7, axes=F)
    #- overlay points and error bars
    
    plotBy(WUE.mean~TDR|Treat,data=toplot,add=T,pch=symbols[i],col=colors[i],cex=2,legend=F,
           panel.first=adderrorbars(x=toplot$TDR,y=toplot$WUE.mean,SE=toplot$WUE.standard.error,direction="updown"))
    
  }
  magaxis(side=1:4,labels=c(0,1,0,0),las=2,cex.axis=1.5)
  title(ylab=expression(WUE[i]~(mu*mol~mol^-1)),xpd=NA,cex.lab=3)
  legend("bottomright",letters[3],bty="n",cex=2)
  legend("topright",legend=c("Cacu","Eusi","Eute","Pira"),pch=symbols,col=colors,ncol=2,cex=2)
  
  #- plot pre-dawn LWP
  for(i in 1:4){
    toplot <- lwp.trt.list[[i]]
    if(i ==1) plot.new();plot.window(xlim=xlims,ylim=c(-10,0));abline(h=0)
    
    #- smoothplot
    smoothplot(TDR, LWP.pd.mean, polycolor=alpha(colors[i],0.3),linecols=colors[i],
               ylab="",cex=0.5,add=T,
               data=toplot, kgam=5, axes=F)
    #- overlay points and error bars
    
    plotBy(LWP.pd.mean~TDR|Treat,data=toplot,add=T,pch=symbols[i],col=colors[i],cex=2,legend=F,
           panel.first=adderrorbars(x=toplot$TDR,y=toplot$LWP.pd.mean,SE=toplot$LWP.pd.standard.error,direction="updown"))
    
  }
  magaxis(side=1:4,labels=c(0,1,0,0),las=2,cex.axis=1.5)
  title(ylab=expression(psi[pd]~(MPa)),xpd=NA,cex.lab=3)
  legend("bottomright",letters[5],bty="n",cex=2)
  
  #- plot Cond
  for(i in 1:4){
    toplot <- ITE.trt.list[[i]]
    if(i ==1) plot.new();plot.window(xlim=xlims,ylim=c(0,0.6));abline(h=0)
    
    #- smoothplot
    smoothplot(TDR, Cond.mean, polycolor=alpha(colors[i],0.3),linecols=colors[i],
               ylab="",cex=0.5,add=T,
               data=toplot, kgam=4, axes=F)
    
    #- overlay points and error bars
    plotBy(Cond.mean~TDR|Treat,data=toplot,add=T,pch=symbols[i],col=c(colors[i]),cex=2,legend=F,
           panel.first=adderrorbars(x=toplot$TDR,y=toplot$Cond.mean,SE=toplot$Cond.standard.error,direction="updown"))
    
  }
  magaxis(side=1:4,labels=c(1,1,0,0),las=1,cex.axis=1.5)
  title(ylab=expression(g[s]~(mol~m^-2~s^-1)),xpd=NA,cex.lab=3)
  title(xlab=expression(theta~(m^3~m^-3)),xpd=NA,cex.lab=3,line=4)
  legend("bottomright",letters[2],bty="n",cex=2)
  
  #- plot Ci/Ca
  for(i in 1:4){
    toplot <- ITE.trt.list[[i]]
    if(i ==1) plot.new();plot.window(xlim=xlims,ylim=c(0,1.5));abline(h=0)
    
    #- smoothplot
    smoothplot(TDR, Ci.Ca.mean, polycolor=alpha(colors[i],0.3),linecols=colors[i],
               ylab="",cex=0.5,add=T,
               data=toplot, kgam=5, axes=F)
    
    #- overlay points and error bars
    plotBy(Ci.Ca.mean~TDR|Treat,data=toplot,add=T,pch=symbols[i],col=colors[i],cex=2,legend=F,
           panel.first=adderrorbars(x=toplot$TDR,y=toplot$Ci.Ca.mean,SE=toplot$Ci.Ca.standard.error,direction="updown"))
    
  }
  magaxis(side=1:4,labels=c(1,1,0,0),las=1,cex.axis=1.5)
  title(ylab=expression(C[i]*"/"*C[a]),xpd=NA,cex.lab=3)
  title(xlab=expression(theta~(m^3~m^-3)),xpd=NA,cex.lab=3,line=4)
  legend("bottomright",letters[4],bty="n",cex=2)
  
  #- plot mid-day LWP
  for(i in 1:4){
    toplot <- lwp.trt.list[[i]]
    if(i ==1) plot.new();plot.window(xlim=xlims,ylim=c(-10,0));abline(h=0)
    
    #- smoothplot
    smoothplot(TDR, LWP.md.mean, polycolor=alpha(colors[i],0.3),linecols=colors[i],
               ylab="",cex=0.5,add=T,
               data=toplot, kgam=5, axes=F)
    
    #- overlay points and error bars
    
    plotBy(LWP.md.mean~TDR|Treat,data=toplot,add=T,pch=symbols[i],col=colors[i],cex=2,legend=F,
           panel.first=adderrorbars(x=toplot$TDR,y=toplot$LWP.md.mean,SE=toplot$LWP.md.standard.error,direction="updown"))
    
  }
  magaxis(side=1:4,labels=c(1,1,0,0),las=1,cex.axis=1.5)
  title(ylab=expression(psi[md]~(MPa)),xpd=NA,cex.lab=3)
  title(xlab=expression(theta~(m^3~m^-3)),xpd=NA,cex.lab=3,line=4)
  legend("bottomright",letters[6],bty="n",cex=2)
  
  if(output==T) dev.copy2pdf(file="Output/Figure3_gx_vs_VWC.pdf")
}
#---------------------------------------------------------------------------------------------------------------------










#---------------------------------------------------------------------------------------------------------------------
#- plot the dependence of g1 and non-stomatal limitation relative to TSW.
#---------------------------------------------------------------------------------------------------------------------
plotBetasG1NSL.TSW <- function(output=F,g1data,NSLdata,g1list,NSLlist){
  #- list for g1 data
  dat.l <- split(g1data,g1data$Species)
  
  #- fit NSL
  NSLdata$TSW <- (NSLdata$TDR)/(max(NSLdata$TDR)-min(NSLdata$TDR)) # normalize TDR data to estimate the transpirable soil water
  dat.l2 <- split(NSLdata,NSLdata$Species)
  
  #- plot normalized g1 and non-stomatal limitation as a function of TSW
  windows(16,16)
  par(mfrow=c(4,2),mar=c(0,0.25,0,0.25),xpd=FALSE,oma=c(4,5,1,5),cex=1.6,cex.axis=0.9,cex.lab=0.9)
  labs <- c("Cacu","Eusi","Eute","Pira")
  count <- 0
  for (i in 1:length(dat.l)){
    newdat <- newdat2 <- data.frame()
    dat.temp <- dat.l[[i]]
    dat.temp$g1norm <- dat.temp$g1/max(dat.temp$g1)
    dat.temp$TSW <- (dat.temp$TDR)/(max(dat.temp$TDR)-min(dat.temp$TDR))
    
    dat.temp$Species <- factor(dat.temp$Species)
    #------------------------------------------------------------------------
    #-- plot g1 vs. TSW
    count <- count+1
    
    plot(g1norm~TSW,data=subset(dat.temp,Treat=="wet"),pch=21,col="black",bg=grey(0.1),axes=F,ylim=c(0,1.18),xlim=c(0,1.05))
    points(g1norm~TSW,data=subset(dat.temp,Treat=="dry"),pch=21,col="black",bg=grey(0.8))
    magaxis(side=c(1:4),labels=c(0,1,0,0),las=1,tcl=0.3,ratio=0.25,majorn=3,cex.axis=0.8)
    if(i==4)  magaxis(side=c(1:4),labels=c(1,1,0,0),las=1,tcl=0.3,ratio=0.25,majorn=3,cex.axis=0.8)
    mtext(labs[i],side=2,xpd=T,cex=1.3,line=1.75)
    if(i==1) legend("bottomright",xpd=NA,legend=c("Wet","Dry"),pch=21,pt.bg=c("black","grey"),ncol=2,cex=0.75)
    if(i==1) title(main="Stomatal",cex.main=0.75,line=0.5,xpd=NA)
    legend("topleft",letters[count],cex=1,inset=-0.1,bty="n")
    
    
    # plot model and SE from bootstrapping
    rm(newdat)
    newdat <- g1list[[2]][[i]]
    lines(wpred~Xval,data=newdat,xpd=F)
    par(xpd=F)
    polygon(x = c(newdat$Xval, rev(newdat$Xval)), y = c(newdat$lower, rev(newdat$upper)),
            col = alpha("grey",0.5),border=NA)
    
    
    #------------------------------------------------------------------------
    #-- repeat, but for NSL
    
    dat.temp2 <- dat.l2[[i]]
    dat.temp2$Species <- factor(dat.temp2$Species)
    count <- count+1
    
    plot(NSL~TSW,data=subset(dat.temp2,Treat=="wet"),pch=21,col="black",bg=grey(0.1),axes=F,ylim=c(0,1.18),xlim=c(0,1))
    points(NSL~TSW,data=subset(dat.temp2,Treat=="dry"),pch=21,col="black",bg=grey(0.8))
    magaxis(side=c(1:4),labels=c(0,0,0,1),las=1,tcl=0.3,ratio=0.25,majorn=3,cex.axis=0.8)
    if(i==4)  magaxis(side=c(1:4),labels=c(1,0,0,0),las=1,tcl=0.3,ratio=0.25,majorn=3,cex.axis=0.8)
    if(i==1) title(main="Non-stomatal",cex.main=0.75,line=0.5,xpd=NA)
    
    legend("topleft",letters[count],cex=1,inset=-0.1,bty="n")
    
    newdat2 <- NSLlist[[2]][[i]]
    lines(wpred~Xval,data=subset(newdat2,wpred<=1))
    polygon(x = c(newdat2$Xval, rev(newdat2$Xval)), y = c(newdat2$lower, rev(newdat2$upper)), 
            col = alpha("grey",0.5),border=NA)
    lines(x=c(max(subset(newdat2,wpred<=1)$Xval),1),y=c(1,1))
    
    
  }
  mtext(expression(Normalized~g[1]),side=2,outer=T,cex=2,las=0,line=3)
  mtext(expression(A/A[e]),side=4,outer=T,cex=2,las=0,line=2.5)
  mtext(expression(TSW~(proportion)),side=1,outer=T,cex=1.5,line=2)
  
  if(output==T) dev.copy2pdf(file="Output/Figure4_Beta_g1andNSL_TSW.pdf")
  
}
#---------------------------------------------------------------------------------------------------------------------

