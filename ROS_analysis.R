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
  



#- fix up the code below to plot leaf-water potential over time
ros2 <- return.gx.vwc.lwp()




#---------------------------------------------------------------------------------------------------------------------
#--- plot pre-dawn and mid-day leaf water potentials over time
lwp.m <- summaryBy(LWP+LWP.md+diff~gxDate+Treat+Species,FUN=c(mean,standard.error),data=ros2)
lwp.m.l <- split(lwp.m,lwp.m$Species)


windows(20,20)
par(mfrow=c(2,1),oma=c(4,4,0.5,2),mar=c(2,3,1,1),cex.axis=1.5,cex.lab=1.8)

startdate <- as.Date(x="2012-10-1",format="%Y-%m-%d")
enddate <- as.Date(x="2013-06-1",format="%Y-%m-%d")
dates <- as.Date(c("2012-10-18","2012-11-15","2012-11-25","2012-12-27","2013-1-26","2013-5-17"),format="%Y-%m-%d")
colors <- rev(c("black","darkgrey"))
pchs <- c(15,16,17,18)
size=1.75
fillcol <- gray(0.95)
#- plot pre-dawn LWP over time for each species
plot.new()
plot.window(xlim=c(startdate,enddate),ylim=c(-10,0))
rect(xleft=dates[1],ybottom=-45,xright=dates[2],ytop=270,col=fillcol) #add rectangles for droughts
rect(xleft=dates[3],ybottom=-45,xright=dates[4],ytop=270,col=fillcol) #add rectangles for droughts
rect(xleft=dates[5],ybottom=-45,xright=enddate+13,ytop=270,col=fillcol) #add rectangles for droughts

#plot pre-dawn
for (i in 1:length(lwp.m.l)){
  dat1 <- lwp.m.l[[i]]
  
  if (i==1) plotBy(LWP.mean~gxDate|Treat,data=dat1,type="b",xlim=c(startdate,enddate),legend=F,ylim=c(-10,0),pch=pchs[i],cex=size,add=T,col=colors,
                   panel.first=adderrorbars(x=dat1$gxDate,y=dat1$LWP.mean,SE=dat1$LWP.standard.error,direction="updown"))
  if (i > 1)plotBy(LWP.mean~gxDate|Treat,data=dat1,type="b",xlim=c(startdate,enddate),legend=F,ylim=c(-10,0),pch=pchs[i],cex=size,add=T,col=colors,
                   panel.first=adderrorbars(x=dat1$gxDate,y=dat1$LWP.mean,SE=dat1$LWP.standard.error,direction="updown"))
  if (i%%2==1) magaxis(c(2,4),labels=c(1,1),frame.plot=TRUE,las=1,ylab="Pre-dawn",tline=0.2,majorn=4,minorn=0)
  if (i<3)axis.Date(side=1,at=seq(startdate,enddate,by="month"),format="%m/%y",labels=FALSE)
  if (i>2)axis.Date(side=1,at=seq(startdate,enddate,by="month"),format="%m/%y")
}
legend("bottomleft",c("Cacu-wet","Eusi-wet","Eute-wet","Pira-wet","Cacu-dry","Eusi-dry","Eute-dry","Pira-dry"),pch=c(pchs,pchs),
       col=c("black","black","black","black","darkgrey","darkgrey","darkgrey","darkgrey"),ncol=2,bg="white",cex=1.4)
#plot mid-day
plot.new()
plot.window(xlim=c(startdate,enddate),ylim=c(-10,0))
rect(xleft=dates[1],ybottom=-45,xright=dates[2],ytop=270,col=fillcol) #add rectangles for droughts
rect(xleft=dates[3],ybottom=-45,xright=dates[4],ytop=270,col=fillcol) #add rectangles for droughts
rect(xleft=dates[5],ybottom=-45,xright=enddate+13,ytop=270,col=fillcol) #add rectangles for droughts
for (i in 1:length(lwp.m.l)){
  dat1 <- lwp.m.l[[i]]
  
  if (i==1) plotBy(LWP.md.mean~gxDate|Treat,data=dat1,type="b",xlim=c(startdate,enddate),legend=F,ylim=c(-10,0),pch=pchs[i],cex=size,add=T,col=colors,
                   panel.first=adderrorbars(x=dat1$gxDate,y=dat1$LWP.md.mean,SE=dat1$LWP.md.standard.error,direction="updown"))
  if (i > 1)plotBy(LWP.md.mean~gxDate|Treat,data=dat1,type="b",xlim=c(startdate,enddate),legend=F,ylim=c(-10,0),pch=pchs[i],cex=size,add=T,col=colors,
                   panel.first=adderrorbars(x=dat1$gxDate,y=dat1$LWP.md.mean,SE=dat1$LWP.md.standard.error,direction="updown"))
  if (i%%2==1) magaxis(c(2,4),labels=c(1,1),frame.plot=TRUE,las=1,ylab="Mid-day",tline=0.2,majorn=4,minorn=0,cex.lab=1.5)
  
  if (i<3)axis.Date(side=1,at=seq(startdate,enddate,by="month"),format="%m/%y",labels=FALSE)
  if (i>2)axis.Date(side=1,at=seq(startdate,enddate,by="month"),format="%m/%y")
}

title(xlab="Date",outer=TRUE,ylab=expression(Leaf~water~potential~(MPa)),cex.lab=2,line=1)
dev.copy2pdf(file="Output/LWP_over_time.pdf")
#---------------------------------------------------------------------------------------------------------------------


#ros2$SWP <- PSIsoil(ros2$TDR,units="percent") #i'm not sure this is the best function to use here
#plotBy(TDR~LWP|sp,data=ros2)

#------------------------------------------------------------------------------------------------------- 
#------------------------------------------------------------------------------------------------------- 
#Great, now let's bin those data and create some averages
#average data into VWC and LWP bins

#LWP ~ soil moisture
ros2$TDR_bin <- cut(ros2$TDR,breaks=c(0,2,4,6,8,12,15,20,25,32,42))
ros2$TDR_bin_mid <- sapply(strsplit(gsub("^\\W|\\W$", "", ros2$TDR_bin), ","), function(x)sum(as.numeric(x))/2) 
plot1.means <- summaryBy(Photo+Cond+LWP+VpdL~Species+TDR_bin_mid,data=ros2,FUN=c(mean,standard.error),keep.names=FALSE)

#Photo ~ LWP
ros2$LWP_bin <- cut(ros2$LWP,breaks=c(-10,-8,-6,-4,-2,-1,-0.8,-0.6,-0.3,0))
ros2$LWP_bin_mid <- sapply(strsplit(gsub("^\\W|\\W$", "", ros2$LWP_bin), ","), function(x)sum(as.numeric(x))/2) 
plot2.means <- summaryBy(Photo+Cond~Species+LWP_bin_mid,data=ros2,FUN=c(mean,standard.error),keep.names=FALSE)



#split into lists for each species
plot1.list <- split(plot1.means,plot1.means$Species)
plot2.list <- split(plot2.means,plot2.means$Species)

#also split the raw data
ros2.list <- split(ros2,ros2$Species)



#---------------------------------------------------------
#--- make 3-panel plot, showing interdependence of cond, VWC, and LWP

bgcols <- c("black",grey(0.5),grey(0.85),"white")
windows(16,12)
par(mfrow=c(3,1),mar=c(5,7,1,1))
size <- 2.5
labsize <- 1.3



#photo vs. VWC
plot.new()
plot.window(xlim=c(0,40),ylim=c(0,1.0))
for (i in 1:length(plot1.list)){
  dat <- plot1.list[[i]]
  dat.all <- ros2.list[[i]]
  adderrorbars(x=dat$TDR_bin_mid,y=dat$Cond.mean,SE=dat$Cond.standard.error,direction="updown",col="black",add=T)
  
  points(Cond~TDR,data=dat.all,pch=21,col="black",bg=bgcols[i],cex=size/2)
  points(Cond.mean~TDR_bin_mid,data=dat,pch=21,type="b",col="black",bg=bgcols[i],cex=size+1,lty=i)
  
}

box()
magaxis(side=c(1,2,3,4),labels=c(1,1,0,0),cex.axis=2,las=1)
title(xlab="VWC (%)",ylab=expression(g[s]~(mol~m^-2~s^-1)),cex.lab=2)

#pdlwp vs. VWC
plot.new()
plot.window(xlim=c(0,40),ylim=c(-10,0))
for (i in 1:length(plot1.list)){
  dat <- plot1.list[[i]]
  dat.all <- ros2.list[[i]]
  
  adderrorbars(x=dat$TDR_bin_mid,y=dat$LWP.mean,SE=dat$LWP.standard.error,direction="updown",col="black",
               add=T)
  points(LWP~TDR,data=dat.all,pch=21,type="p",col="black",bg=bgcols[i],cex=size/2)
  points(LWP.mean~TDR_bin_mid,data=dat,pch=21,type="b",col="black",bg=bgcols[i],cex=size+1,lty=i)
}
box()
magaxis(side=c(1,2,3,4),labels=c(1,1,0,0),cex.axis=2,las=1)
title(xlab="VWC (%)",ylab=expression(psi[l-pd]~("MPa")),cex.lab=2)
legend("bottomright",legend=c("Cacu","Eusi","Eute","Pira"),pch=21,col="black",pt.bg=bgcols,pt.cex=size+1,cex=1.8,lty=c(1:4))


#Cond vs. pdlwp
plot.new()
plot.window(xlim=c(-10,0),ylim=c(0,1))
for (i in 1:length(plot2.list)){
  dat <- plot2.list[[i]]
  dat.all <- ros2.list[[i]]
  
  adderrorbars(x=dat$LWP_bin_mid,y=dat$Cond.mean,SE=dat$Cond.standard.error,direction="updown",col="black",add=T)
  
  points(Cond~LWP,data=dat.all,pch=21,type="p",col="black",bg=bgcols[i],cex=size/2)
  points(Cond.mean~LWP_bin_mid,data=dat,pch=21,type="b",col="black",bg=bgcols[i],cex=size+1,lty=i)
  
}
box()
magaxis(side=c(1,2,3,4),labels=c(1,1,0,0),cex.axis=2,las=1)

title(xlab=expression(psi[l-pd]~("MPa")),ylab=expression(g[s]~(mol~m^-2~s^-1)),cex.lab=2)
#dev.copy2pdf(file="./Output/Cond_VWC_PD-LWP_3-panel.pdf")
