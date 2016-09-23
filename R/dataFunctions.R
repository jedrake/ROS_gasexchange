#-----------------------------------------------------------------------------------------
#- Functions to read, process, and return data in the ROS project
#-----------------------------------------------------------------------------------------





#-----------------------------------------------------------------------------------------
#read and process handheld TDR data
get.TDR.handheld <- function(filterBlanks=1){
  
  #- read in the data from HIEv
  TDR1 <- read.csv("Data/ROS_MD_PM_SOILMOIST-HANDHELD_L2.csv")
  
  #- Filter out the blanks and non-shelter measurements
  if (filterBlanks==1){
    TDRdat <- subset(TDR1, Species !="blank" & Pot < 700)
    TDRdat$Species <- droplevels(TDRdat$Species)} #remove data from blank pots, keep only data from shelters (pots 100-699)
  if (filterBlanks!=1){
    TDRdat <- TDR1
  }
  
  TDRdat$Treat <- tolower(TDRdat$Treat)
  TDRdat$Date <- as.Date(TDRdat$Date)
  
  #- convert from % to m3 m-3
  TDRdat$TDR <- TDRdat$TDR/100
  
  return(TDRdat)
}
#-----------------------------------------------------------------------------------------






#-----------------------------------------------------------------------------------------
#- download the automated TDR measures of soil moisture.
#-----------------------------------------------------------------------------------------
get.TDR.logger <- function(plotson=0,startdate="2012-08-01",enddate="2013-12-1"){


  
  #find all of the files with soil moisture at ROS
  dfr <- downloadTOA5(filename="ROS_[0-9]_ROS[0-9][0-9]Table15min",maxnfiles=120,topath="Data/VWC/",
                      startDate=startdate,endDate=enddate)
  #cachefile="TDR_ROS_Logger.RData")
  dfr$Plot <- as.numeric((substr(dfr$Source,start=10,stop=11)))
  attr(dfr$DateTime,"tzone") <- "GMT"
  
  
  #reshape data into long format with VWC data from shelters only
  dfr2 <- subset(dfr,Plot<7,select=c("DateTime","Date","Plot","VW.1.","VW.2.","VW.3.","VW.4.","VW.5.",
                                     "VW.6.","VW.7.","VW.8."))
  dfr3 <- melt(dfr2,id.vars=c("DateTime","Date","Plot"))
  names(dfr3)[5] <- "VWC"
  dfr3$sensor <- as.numeric(substr(dfr3$variable,start=4,stop=4))
  
  
  #read in the identifying data as to which sensor is paired with which pot
  ids <- read.csv("./Data/ROS_MD_PM_SOILMOIST-CODES.csv")
  
  
  #merge data with ids. Note that there are ~2 million rows in this dataframe. Only take dates range that I want
  dfr4.p <- plyr::join(dfr3,ids,type="full")
  
  #dfr4 <- merge(dfr3,ids,by=c("Plot","sensor")) #removed, plyr join is faster
  maxdate <- as.POSIXct(x="2013-05-20",format="%Y-%m-%d",tz="GMT")
  mindate <- as.POSIXct(x="2012-5-01",format="%Y-%m-%d",tz="GMT")
  dfr5 <- subset(dfr4.p,DateTime>mindate & DateTime<maxdate)
  
  #make daily averages
  #VWC.day <- summaryBy(VWC~Plot+sensor+pot+Treat+Species+Date,FUN=mean,data=dfr5,keep.names=TRUE) # removed. dplyr method is so much faster
  VWC.day <- dplyr::summarize(dplyr::group_by(dfr5,Plot,sensor,pot,Treat,Species,Date),
                              VWC = mean(VWC))
  if(plotson==1){
    #make a plot of VWC for each plot
    pdf("./Ouput/VWC by species.pdf",onefile=TRUE)
    VWC.day.list <- split(VWC.day,VWC.day$Species)
    startdate <- as.Date(x="2012-07-01",format="%Y-%m-%d")
    enddate <- as.Date(x="2013-10-01",format="%Y-%m-%d")
    colors <- c("red","blue")
    for (i in 1:length(VWC.day.list)){
      
      dat <- VWC.day.list[[i]]
      plotBy(VWC~Date|Treat,ylim=c(0,0.4),axes=FALSE,col=colors,
             data=dat,type="p",legend=TRUE,main=VWC.day.list[[i]]$Species[1])
      box()
      axis(2)
      axis(4)
      axis.Date(side=1,at=seq(startdate,enddate,by="month"),format="%d-%m-%y")
      
    }
    dev.off()
  }
  return(VWC.day)
}
#-----------------------------------------------------------------------------------------







#-----------------------------------------------------------------------------------------
#- standard error, while removing NA's by default
standard.error <- function(dat,na.rm=F,...){
  if(na.rm==T){
    dat <- subset(dat,is.na(dat)==F)
  }
  std <- sd(dat)
  n <- length(dat)
  se <- std/sqrt(n)
  return(se)
}
#-----------------------------------------------------------------------------------------






#-----------------------------------------------------------------------------------------
#- Adds error bars to a plot
#-----------------------------------------------------------------------------------------
adderrorbars <- function(x,y,SE,direction,barlen=0.04,...){
  
  if(length(direction)>1)stop("direction must be of length one.")
  if(direction == "updown")
    direction <- c("up","down")
  else if(direction == "rightleft" | direction == "leftright")direction <- c("left","right")
  
  if("up" %in% direction)
    arrows(x0=x, x1=x, y0=y, y1=y+SE, code=3, angle=90, length=barlen,...)
  if("down" %in% direction) 
    arrows(x0=x, x1=x, y0=y, y1=y-SE, code=3, angle=90, length=barlen,...)
  if("left" %in% direction) 
    arrows(x0=x, x1=x-SE, y0=y, y1=y, code=3, angle=90, length=barlen,...)
  if("right" %in% direction)
    arrows(x0=x, x1=x+SE, y0=y, y1=y, code=3, angle=90, length=barlen,...)  
  
}
#-----------------------------------------------------------------------------------------






#-----------------------------------------------------------------------------------------
#--- This functions "does the Joey" and calculates the non-stomatal limitation of 
#         photosyntheis as in Zhou et al. 2013 Tree Phys.
#-----------------------------------------------------------------------------------------
returnVcmaxa <- function(){
  #--- get the leaf-level gas exchange data
  #ros2 <- return.gx.vwc.lwp()
  ros2 <- return.gx.vwc()
  
  #- remove gx data in ros2 where Ci < 5 or Ci > 800
  ros3 <- subset(ros2,Ci>5 & Ci<1000)
  
  #- define function to return Vcmax based on a single point A:Ci curve as in Zhou et al. 2013 Tree Phys.
  returnVcmax <- function(A,Ci,Tleaf){
    Tk <- Tleaf + 273.15
    R <- 8.314                                                       #universal gas constant in J mol-1 K-1
    gammastar <- 42.75*exp( (37830*(Tk-298.15))/(298.15*R*Tk) )            #calculate gamma star as in bernacchi et al. 2001
    #Ko <- exp(20.3-36.38/(R/1000*Tk))                                     #calculate Ko as in Bernacchi et al. 2001
    Ko <- 278.4*exp((36380*(Tk-298.15)/(298.15*R*Tk)))
    #Kc <- exp(38.05-79.43/(R/1000*Tk))                                    #calculate Kc as in Barnacchi et al. 2001
    Kc <- 404.9*exp((79403*(Tk-298.15)/(298.15*R*Tk)))
    Km <- Kc*(1+210/Ko)                                              # calculate Km. Assumes [O2] = 210 mmol mol-1
    
    #Vcmax <- A*(Ci+Km)/(Ci-gammastar - 0.015)     #calculate apparent Vcmax. De Kauwe (2016) eq. 3
    Vcmax <- A/((Ci-gammastar)/(Ci+Km) - 0.015)     #calculate apparent Vcmax. De Kauwe (2016) eq. 3
    
    return(Vcmax)
  }
  
  #- get apparent Vcmax for each observation
  ros3$Vcmax_a <- returnVcmax(A=ros3$Photo,Ci=ros3$Ci,Tleaf=ros3$Tleaf)
  #ros3 <- subset(ros3,Vcmax_a>5 & Vcmax_a<600)
  
  #- get the "maximum" Vcmax as the 50th percentile of the apparent Vcmax of the well-watered treatments
  ros3.list <- split(ros3,ros3$Species)
  Species <- c()
  Vcmax_max <- c()
  for(i in 1:length(ros3.list)){
    dat <- subset(ros3.list[[i]],Treat=="wet" & TDR > 15 & TDR < 25 & Ci > 100)
    dat <- ros3.list[[i]]
    Species[i] <- as.character(dat$Species[1])
    Vcmax_max[i] <- unname(quantile(dat$Vcmax_a,probs=0.5))
    #Vcmax_max[i] <- max(dat$Vcmax_a)
  }
  df2 <- data.frame(Species=Species,Vcmax_max=Vcmax_max)
  #df2$Vcmax_max <- max(df2$Vcmax_max)
  
  ros4 <- merge(ros3,df2,by=c("Species"))
  ros4$Photo_a <- Photosyn(VPD=ros4$VpdL,Ca=ros4$CO2S,PPFD=ros4$PARi,Tleaf=ros4$Tleaf,Ci=ros4$Ci,
                           Vcmax=ros4$Vcmax_max,Jmax=1.6*ros4$Vcmax_max,Tcorrect=T)$ALEAF
  ros4$NSL <- with(ros4,Photo/Photo_a)
  
  #- remove a really troublesome eute point and a cacu point
  #ros4[258,] <- NA #- these were the bad points when using return.vwc.lwp()
  #ros4[42,] <- NA
  ros4[1003:1004,] <- NA
  ros4[629,] <- NA
  ros4[479,] <- NA
  ros4[975,] <- NA
  
  ros5 <- ros4[complete.cases(ros4),]
  ros5$TDR <- ros5$TDR/100
  
  boxplot(NSL~Species,data=ros5)
  
  
  return(ros5)
}
#-----------------------------------------------------------------------------------------







#-----------------------------------------------------------------------------------------
#- Read in the leaf water potential data
#-----------------------------------------------------------------------------------------
read.lwp <- function(){
  #LWP data from ROS
  library(reshape)
  library(plotBy)
  
  #read in Danielle's LWP data. Note that they are in the "wide" format, and some dates have missing data
  lwp <- read.csv("./Data/LWP/ROS Water potential compiled_JED.csv")
  names(lwp)[2] <- "sp"
  
  #create species lookup table
  sp <- c("ES","ET","P","C")
  Species <- c("eusi","eute","pira","cacu")
  spconv <- data.frame(sp,Species)
  
  #create treatment lookup table
  Treatment <- c("D","W")
  Treat <- c("dry","wet")
  treatconv <- data.frame(Treatment,Treat)
  
  #melt into "long" data format
  lwp.m <- melt(lwp,c("Shelter","sp","Treatment","Pot"))
  names(lwp.m)[6] <- "LWP"
  lwp.m$Date <- as.Date(substr(lwp.m$variable,start=3,stop=8),format="%d%m%y")
  lwp.m$Type <- factor(substr(lwp.m$variable,start=1,stop=2))
  lwp.m <- merge(lwp.m,spconv,by="sp")
  lwp.m <- merge(lwp.m,treatconv,by="Treatment")
  
  
  #read in the LWP data measured for the DRUID campaign, make formatting the same as the other data
  lwp_druid <- read.csv("./Data/LWP/DruidLWP_26Sept2013.csv")
  lwp_druid$Treatment <- toupper(lwp_druid$trt)
  names(lwp_druid)[3:7] <- c("Species","trt","Pot","PD","MD")
  lwp_druid$Date <- as.Date(lwp_druid$date,format="%d/%m/%Y")
  lwp_druid2 <- subset(lwp_druid,select=c("Species","Pot","Treatment","Date","PD","MD"))
  lwp_druid3 <- melt(lwp_druid2,c("Species","Pot","Treatment","Date"),variable_name="Type")
  names(lwp_druid3)[6] <- "LWP"
  lwp_druid3$LWP <- lwp_druid3$LWP*-1
  
  #merge the LWP dataframes
  commoncols <- intersect(colnames(lwp.m),colnames(lwp_druid3))
  lwp2 <- rbind(subset(lwp.m,select=commoncols),subset(lwp_druid3,select=commoncols))
  lwp3 <- merge(lwp2,treatconv)[,2:7]
  
  
  #read in the extra LWP data that was measured on a different set of pots, on 4.10.12
  extras <- read.csv("./Data/LWP/LWP_10April12_ROS.csv")
  extras$Date <- as.Date(extras$Date,format="%m/%d/%Y")
  
  #merge the extra dataframes
  commoncols <- intersect(colnames(lwp3),colnames(extras))
  lwp4 <- rbind(subset(lwp3,select=commoncols),subset(extras,select=commoncols))
  
  return(lwp4)
}
#-----------------------------------------------------------------------------------------






#-----------------------------------------------------------------------------------------
#-read in the processed gas exchange data, combine with soil moisture
#-----------------------------------------------------------------------------------------
return.gx.vwc <- function(){

  gx <- read.gx.ROS()
  vwc <- na.omit(get.TDR.handheld())
  
  gx.dry <- subset(gx,Treat=="dry")
  vwc.dry <- subset(vwc,Treat=="dry")
  
  gx.wet <- subset(gx,Treat=="wet")
  vwc.wet <- subset(vwc,Treat=="wet")
  
  #-------------------------------------------------------------------------------------
  # Find the closest date when VWC was measured, to match these obsevations with the gx measurements.

  vwcdates.wet <- summaryBy(TDR~Date,FUN=length,data=vwc.wet,keep.names=TRUE)
  vwcdates.wet <- subset(vwcdates.wet,TDR>15) #get rid of the dates with few observaitons
  #get the dates gx and vwc was measured in the dry treatments
  gxdates.dry <- levels(as.factor(gx.dry$Date))
  vwcdates.dry <- levels(as.factor(vwc.dry$Date))
  #get the dates gx and vwc was measured in the wet treatments
  gxdates.wet <- levels(as.factor(gx.wet$Date))
  vwcdates.wet <- levels(as.factor(vwcdates.wet$Date))
  
  
  #this function returns the closest date for a searchDate, given a vector of dates to search (dateList)
  # I use this function to help merge the TDR data with the GX data, as these were often measured on adjacent days
  closestDate <- function(searchDate, dateList){
    dist2date <- abs(as.Date(dateList) - as.Date(searchDate))
    closest <- which(min(dist2date)==dist2date)
    return(dateList[closest])
  }
  
  #find the closest date of TDR measurements to the gas exchange date.
  # Note that this is done separately for wet and dry, as vwc was measured less
  # often in the wet treatment
  gx.dry$VWCdate <- NA
  for (i in 1:nrow(gx.dry)){
    gx.dry$VWCdate[i] <- closestDate(searchDate=gx.dry$Date[i],dateList=vwcdates.dry)
  }
  gx.wet$VWCdate <- NA
  for (i in 1:nrow(gx.wet)){
    gx.wet$VWCdate[i] <- closestDate(searchDate=gx.wet$Date[i],dateList=vwcdates.wet)
  }
  #put the wet and dry subsets back together
  gx2 <- rbind(gx.wet,gx.dry)
  gx2$VWCdate <- as.POSIXct(gx2$VWCdate,format="%Y-%m-%d",tz="GMT")
  
  
  #merge vwc and gx datasets together, using a date by species by treatment average for the vwc data
  names(gx2)[4] <- "gxDate"
  names(gx2)[57] <- "Date"
  
  #average the VWC data
  vwc.treat.avg <- summaryBy(TDR+Ladd~Date+Species+Treat,
                             FUN=mean,keep.names=TRUE,data=vwc)
  
  #- merge the gas-exchange data (gx2) with a treatment by date by species mean of the TDR measurements
  gx3 <- merge(gx2,vwc.treat.avg,by=c("Date","Species","Treat"))
  
  return(gx3)
}
#-----------------------------------------------------------------------------------------






#-----------------------------------------------------------------------------------------
#- read in the L1 ROS gas exchange data (manually cleaned by RS, but not averaged up by species or anything)
#-----------------------------------------------------------------------------------------
read.gx.ROS <- function(){

  
  #-- the old way of doing things...
  
  
  # #find the right files
  # datadir <- "oldData/ROS gas exchange data/Cleaned CSV & Text files/"
  # datapath <- file.path(getwd(),datadir)
  # files <- list.files(datapath,pattern="*.csv")
  # 
  # #read them into a dataframe "dat"
  # dat <- list()
  # for (i in 1:length(files)){
  #   dat[[i]] <- read.csv(file.path(datapath,files[i])) # read in each file and extract the date from the file name
  #   dat[[i]]$Date <- str_extract(str_extract(files[i],"[0-9]+[A-z]+[0-9]+.csv"),"[0-9]+[A-z]+[0-9]+") #- had to remove "*" at beginning of first regex, 15/5/15
  # }
  # names(dat[[4]]) <- names(dat[[3]]) #for some reason, "pot." was not named quite the same thing in the fourth file
  # names(dat[[6]]) <- names(dat[[3]])
  # dat <- do.call(rbind,dat)
  # 
  # #do some renaming and creating the right factor levels
  # names(dat)[2:3] <- c("Treat","Pot")
  # dat$Species <- factor(tolower(substr(dat$Species,start=1,stop=4)))
  # dat$Treat <- factor(tolower(dat$Treat))
  # #   dat$Pot <- factor(dat$Pot)
  # dat$Date <- as.Date(dat$Date,format="%d%b%y")
  # 
  # #average sub-replicate logs
  # dat2 <- summaryBy(.~Species+Treat+Pot+Date,data=dat,FUN=mean,keep.names=TRUE)
  # toremove214 <- which(dat2$Pot == 214 & dat2$Date == "2013-03-21") #these two observations are crazy outliers and will be removed
  # toremove219 <- which(dat2$Pot == 219 & dat2$Date == "2013-03-21")
  # dat2[toremove214,] <- NA
  # dat2[toremove219,] <- NA
  # 
  # dat3 <- subset(dat2,is.na(Date)==FALSE)
  # 
  
  #--- alternatively, read in the data on HIEv
  dat.hiev <- read.csv("Data/ROS_MD_GX-Asat_120410-130516_L1.csv")
  
  #- remove first column (names)
  dat.hiev$X <- NULL
  
  #- remove a few outliers
  toremove214 <- which(dat.hiev$Pot == 214 & dat.hiev$Date == "2013-03-21") #these two observations are crazy outliers and will be removed
  toremove219 <- which(dat.hiev$Pot == 219 & dat.hiev$Date == "2013-03-21")
  dat.hiev[toremove214,] <- NA
  dat.hiev[toremove219,] <- NA
  dat3 <- subset(dat.hiev,is.na(Date)==FALSE)
  
  #- reformat variables
  dat3$Date <- as.Date(dat3$Date)
   
  return(dat3)
}
#-----------------------------------------------------------------------------------------






#-----------------------------------------------------------------------------------------
#- merge leaf water potential data with the gas exchange and VWC data
#-----------------------------------------------------------------------------------------
return.gx.vwc.lwp <- function(){


  #get the raw data
  ros <- return.gx.vwc()
  lwp <- read.lwp()
  
  #--------------- process LWP data a bit
  names(lwp)[3] <- "LWPdate"
  lwp.pd <- subset(lwp,Type=="PD" & is.na(LWP)==FALSE)
  lwp.pd$LWP <- -1*lwp.pd$LWP
  lwp.md <- subset(lwp,Type=="MD" & is.na(LWP)==FALSE)
  lwp.md$LWP <- -1*lwp.md$LWP
  
  names(lwp.md)[2] <- "LWP.md"
  lwp2 <- merge(lwp.pd,lwp.md,by=c("Pot","LWPdate","Treat","Species"),all=FALSE)
  lwp2$diff <- with(lwp2,abs(LWP.md)-abs(LWP))
  lwp2$Type.x <- lwp2$Type.y <- NULL
  
  #----------------
  
  
  
  #---------------- Match Dates
  #the dates aren't exactly matched, BUT they're damn close
  #with(ros,plot(CO2S~gxDate))
  #rug(lwp.pd$LWPDate) 
  
  #so, I need to match dates using the same idea as I did in "Merge_gx_TDR.R"
  #this function returns the closest date for a searchDate, given a vector of dates to search (dateList)
  closestDate <- function(searchDate, dateList){
    dist2date <- abs(as.Date(dateList) - as.Date(searchDate))
    closest <- which(min(dist2date)==dist2date)
    return(dateList[closest])
  }
  #grab the dates that leaf water potential was measured
  LWPdates <- levels(as.factor(lwp.pd$LWPdate))
  
  #find the closest date of LWP measurements to the gas exchange date.
  ros$LWPdate <- NA
  for (i in 1:nrow(ros)){
    ros$LWPdate[i] <- closestDate(searchDate=ros$gxDate[i],dateList=LWPdates)
  }
  ros$LWPdate <- as.Date(ros$LWPdate)
  
  # merge the dataframes
  ros2 <- merge(ros,lwp2,by=c("LWPdate","Pot","Species","Treat"))
  ros2$LWP[which(ros2$diff >= 4)] <- NA
  ros2$LWP.md[which(ros2$diff <= -2)] <- NA
  ros2$diff[which(ros2$diff <= -2)] <- NA
  return(ros2)
}
#-----------------------------------------------------------------------------------------






#-----------------------------------------------------------------------------------------
#- fits Belinda's optimal stomatal model to the data, returns g1 parameters with confidence intervals
#-----------------------------------------------------------------------------------------
returng1 <- function(){
  
  #get the gas exchange data
  ros <- return.gx.vwc()
  #ros <- return.gx.vwc.lwp() # note, this cuts out more than half the data.
  
  #--- na-fill some crazy cacu data from the dry-down
  ros$Cond[which(ros$Species=="cacu" & ros$gxDate==as.Date("2013-3-26") & ros$Treat=="dry" & ros$Cond >0.3)] <- NA
  ros$Cond[which(ros$Species=="cacu" & ros$gxDate==as.Date("2013-3-21") & ros$Treat=="dry" & ros$Cond >0.3)] <- NA
  ros$Cond[which(ros$Species=="cacu" & ros$gxDate==as.Date("2013-4-4") & ros$Treat=="dry" & ros$Cond >0.2)] <- NA
  
  
  #split data into a list, for each species and treatment on each date
  ros.date.list <- split(ros,paste(ros$Species,ros$gxDate,ros$Treat))
  
  #fit the optimal stomatal model for each set of data in the list
  fits.list <- list()
  Species <- vector(mode="character")
  Date <- as.Date(NA)
  Treat <- vector(mode="character")
  TDR <- c()
  pdf("Output/ROS_g1plots.pdf")
  for (i in 1:length(ros.date.list)){
    dat <- ros.date.list[[i]]
    fits.list[[i]] <- nls(Cond ~ 1.6*(1+g1/sqrt(VpdL))*(Photo/CO2S),start=list(g1=4),data=dat,algorithm="port",
                          lower=c(0),upper=(10))
    Date[i] <- as.Date(dat$gxDate[1])
    Species[i] <- as.character(dat$Species[1])
    Treat[i] <- as.character(dat$Treat[1])
    TDR[i] <- mean(dat$TDR)
    
    # plot diagnostics
    dat$opti <- with(dat,Photo/(CO2S*sqrt(VpdL)))
    plot(Cond~opti,data=dat,ylim=c(-0.05,0.5),xlim=c(-0.05,0.15))
    lm1 <- lm(Cond~0+opti,data=dat)
    abline(h=0);abline(v=0);abline(lm1)
    title(main=paste(dat$Date[1],dat$Species[1],dat$Treat[1],round(coef(fits.list[[i]]),2)))
    
  }
  dev.off()
  
  #extract the model fits
  g1pars <-as.data.frame(do.call(rbind,(lapply(fits.list,FUN=coefficients))))
  g1pars$Species <- factor(Species)
  g1pars$Treat <- factor(Treat)
  g1pars$Treat <- relevel(g1pars$Treat,ref=2)
  g1pars$TDR <- TDR/100
  g1pars$Date <- as.Date(Date)
  
  
  #get confidence interval of g1 parameters
  g1conf <- as.data.frame(do.call(rbind,lapply(fits.list,FUN=confint)))
  g1pars$g1_2.5 <- g1conf[,1]
  g1pars$g1_97.5 <- g1conf[,2]
  g1pars$names <- with(g1pars,paste(Species,Treat,sep=": "))
  g1pars$TDRjit <- jitter(g1pars$TDR,factor=0.5)
  
  #gap fill the data with uncertain g1 values
  g1pars$g1[which(g1pars$g1_2.5<=0)] <- 0
  g1pars$g1_2.5[which(is.na(g1pars$g1_2.5))] <- 0
  g1pars$g1_97.5[which(is.na(g1pars$g1_97.5))] <- 0
  
  return(g1pars)
}
#-----------------------------------------------------------------------------------------







#-----------------------------------------------------------------------------------------
#-- function to return d13C data from the multiple drought experiment at the ROS

#-- notes: As a reminder, the samples ending -21 were collected on 21/12/13, after the second dry down, 
#        while those ending 26 were collected at the end of the final dry down (26/3/14). If you recall, 
#        we tagged the twigs just before the first dry down began, so the first set of “leaves” sampled had formed 
#        at some point during - or between - the first two dry downs, while those collected later had developed 
#        under even more varied (and extreme) conditions.

# I found two possible mislabellings: One of the wet cacu in shelter 3 was labeled 318 on the first sampling and 378 on the second.
# One of the wet Eute in shelter 5 was labeled 505 on the first sampling and 503 on the second. These were confirmed with Sally
# and changed in the original dataset (the correct labels were 318 and 503)
#-----------------------------------------------------------------------------------------
get_d13C <- function(frac1=27,frac2=4.4,aird13C=-8){
  #frac1 is the kinetic fractionation of photosynthesis
  #frac2 is the kinetic fractionation of differential diffusion of 12CO2 and 13Co2
  #aird13C is the assumed isotopic composition of the source CO2 (atmosphere)
  
  #---- read in the data, do a little processing
  d1 <- read.csv("Data/Isotopes/ROS isotope data.csv")
  names(d1)[6] <- "Treat"
  d1$plant <- as.factor(substr(d1$ID,start=1,stop=3))
  d1$Treat <- relevel(d1$Treat,ref="wet")
  d1$Date <- as.Date(d1$Date,format="%d/%m/%Y")
  d1$bigDelta <- (aird13C/1000-d1$deltaC/1000)/(1+d1$deltaC/1000) #calculate discrimination, assume atmosphere is -8 permil
  d1$CiCa <- with(d1,(bigDelta-frac2/1000)/(frac1/1000-frac2/1000))
  d1$time <- as.factor(d1$time)
  d1$Shelter <- as.factor(d1$Shelter)
  
  return(d1)
}
#-----------------------------------------------------------------------------------------














#---------------------------------------------------------------------------------------------------------------------
#- Function to fit beta models to data (either g1 or NSL datasets)
#  Assumes data have been averaged across dates and species
#---------------------------------------------------------------------------------------------------------------------
fitBeta_TDR <- function(dat,type="g1",startlist = list(Xlow = 0.0, Xhigh=0.5, q = 0.6)){
  
  #- split into list of species
  dat.l <- split(dat,dat$Species)
  
  #- preallocate empty lists to catch output
  fit.sp <- newdat <-  list()
  
  #- loop over each element of the list, fit data
  for (i in 1:length(dat.l)){
    dat.temp <- dat.l[[i]]
    dat.temp$Species <- factor(dat.temp$Species)
    
    #- fit the model
    if (type=="g1") dat.temp$Yval <- dat.temp$g1/max(dat.temp$g1)
    if (type=="NSL") dat.temp$Yval <- dat.temp$NSL
    fit.sp[[i]] <- nls(Yval ~ ((TDR-Xlow)/(Xhigh-Xlow))^q,start=startlist,data=dat.temp,algorithm="port",
                       lower=c(0,0.01,0.01),upper=c(0.007,0.6,3))
    
    
    # get predicted values and 95% confidence intervals by bootstrapping
    if (type=="g1") newdat[[i]] <- expand.grid(Species=levels(dat.temp$Species), TDR=seq(from=coef(fit.sp[[i]])["Xlow"]+0.01,to=max(dat.temp$TDR),length.out=99),lower=NA,upper=NA)
    if (type=="NSL") newdat[[i]] <- expand.grid(Species=levels(dat.temp$Species), TDR=seq(from=coef(fit.sp[[i]])["Xlow"]+0.01,to=coef(fit.sp[[i]])["Xhigh"],length.out=99),lower=NA,upper=NA)
    newdat[[i]]$wpred <- predict(fit.sp[[i]],newdat[[i]],level=0,se.fit=T)
    
    
    rm(b)
    b <- bootCase(fit.sp[[i]],B=300)
    for(j in 1:nrow(newdat[[i]])){
      rm(b02)
      b02 <- ((newdat[[i]]$TDR[j]-b[,"Xlow"])/(b[,"Xhigh"]-b[,"Xlow"]))^b[,"q"]
      
      newdat[[i]]$lower[j] <- unname(quantile(b02,probs=c(0.025,0.975)))[1]
      newdat[[i]]$upper[j] <- unname(quantile(b02,probs=c(0.025,0.975)))[2]
      
    }
    
  }
  
  #- return a list of two lists
  return(list(fit.sp,newdat))
}
#---------------------------------------------------------------------------------------------------------------------





#---------------------------------------------------------------------------------------------------------------------
#- Function to write-out a CSV file of the table of parameter values
#---------------------------------------------------------------------------------------------------------------------
writeBetaParams <- function(){
  #- extract summary tables of beta models for g1 and VWC
  summaries <- lapply(c(g1_TDR_beta[[1]],g1_TDR_beta2[[1]]),summary)
  df1 <- data.frame(summaries[[1]]$coefficients[,1:2])
  df2 <- data.frame(summaries[[2]]$coefficients[,1:2])
  df3 <- data.frame(summaries[[3]]$coefficients[,1:2])
  df4 <- data.frame(summaries[[4]]$coefficients[,1:2])
  df5 <- data.frame(summaries[[5]]$coefficients[,1:2])
  df6 <- data.frame(summaries[[6]]$coefficients[,1:2])
  df7 <- data.frame(summaries[[7]]$coefficients[,1:2])
  df8 <- data.frame(summaries[[8]]$coefficients[,1:2])
  
  #- extract vectors of parameters and se's
  vals1 <- paste(sprintf("%.2f",round(df1$Estimate,2))," (",sprintf("%.2f",round(df1$Std..Error,2)),")",sep="")
  vals2 <- paste(sprintf("%.2f",round(df2$Estimate,2))," (",sprintf("%.2f",round(df2$Std..Error,2)),")",sep="")
  vals3 <- paste(sprintf("%.2f",round(df3$Estimate,2))," (",sprintf("%.2f",round(df3$Std..Error,2)),")",sep="")
  vals4 <- paste(sprintf("%.2f",round(df4$Estimate,2))," (",sprintf("%.2f",round(df4$Std..Error,2)),")",sep="")
  vals5 <- paste(sprintf("%.2f",round(df5$Estimate-c(11,11,0),2))," (",sprintf("%.2f",round(df5$Std..Error,2)),")",sep="")
  vals6 <- paste(sprintf("%.2f",round(df6$Estimate-c(11,11,0),2))," (",sprintf("%.2f",round(df6$Std..Error,2)),")",sep="")
  vals7 <- paste(sprintf("%.2f",round(df7$Estimate-c(11,11,0),2))," (",sprintf("%.2f",round(df7$Std..Error,2)),")",sep="")
  vals8 <- paste(sprintf("%.2f",round(df8$Estimate-c(11,11,0),2))," (",sprintf("%.2f",round(df8$Std..Error,2)),")",sep="")
  g1table <- data.frame(rbind(vals1,vals2,vals3,vals4,vals5,vals6,vals7,vals8),
                        row.names=c("Cacu","Eusi","Eute","Pira","Cacu2","Eusi2","Eute2","Pira2"))
  names(g1table) <- c("Xo","Xh","q")
  g1table$Xvar <- c(rep("VWC",4),rep("LWP",4))
  g1table$Yvar <- "g1"
  
  #---
  #- extract summary tables of beta models for NSL and VWC
  summaries2 <- lapply(c(NSL_TDR_beta[[1]],NSL_TDR_beta2[[1]]),summary)
  df1 <- data.frame(summaries2[[1]]$coefficients[,1:2])
  df2 <- data.frame(summaries2[[2]]$coefficients[,1:2])
  df3 <- data.frame(summaries2[[3]]$coefficients[,1:2])
  df4 <- data.frame(summaries2[[4]]$coefficients[,1:2])
  df5 <- data.frame(summaries2[[5]]$coefficients[,1:2])
  df6 <- data.frame(summaries2[[6]]$coefficients[,1:2])
  df7 <- data.frame(summaries2[[7]]$coefficients[,1:2])
  df8 <- data.frame(summaries2[[8]]$coefficients[,1:2])
  
  #- extract vectors of parameters and se's
  vals1 <- paste(sprintf("%.2f",round(df1$Estimate,2))," (",sprintf("%.2f",round(df1$Std..Error,2)),")",sep="")
  vals2 <- paste(sprintf("%.2f",round(df2$Estimate,2))," (",sprintf("%.2f",round(df2$Std..Error,2)),")",sep="")
  vals3 <- paste(sprintf("%.2f",round(df3$Estimate,2))," (",sprintf("%.2f",round(df3$Std..Error,2)),")",sep="")
  vals4 <- paste(sprintf("%.2f",round(df4$Estimate,2))," (",sprintf("%.2f",round(df4$Std..Error,2)),")",sep="")
  vals5 <- paste(sprintf("%.2f",round(df5$Estimate-c(11,11,0),2))," (",sprintf("%.2f",round(df5$Std..Error,2)),")",sep="")
  vals6 <- paste(sprintf("%.2f",round(df6$Estimate-c(11,11,0),2))," (",sprintf("%.2f",round(df6$Std..Error,2)),")",sep="")
  vals7 <- paste(sprintf("%.2f",round(df7$Estimate-c(11,11,0),2))," (",sprintf("%.2f",round(df7$Std..Error,2)),")",sep="")
  vals8 <- paste(sprintf("%.2f",round(df8$Estimate-c(11,11,0),2))," (",sprintf("%.2f",round(df8$Std..Error,2)),")",sep="")
  NSLtable <- data.frame(rbind(vals1,vals2,vals3,vals4,vals5,vals6,vals7,vals8),
                         row.names=c("Cacu","Eusi","Eute","Pira","Cacu2","Eusi2","Eute2","Pira2"))
  names(NSLtable) <- c("Xo","Xh","q")
  NSLtable$Xvar <- c(rep("VWC",4),rep("LWP",4))
  NSLtable$Yvar <- "A/Ae"
  
  table <- rbind(g1table,NSLtable)[,c(4:5,1:3)]
  write.csv(table,"Output/table2.csv",row.names=T)
}
#-------------------------------------------------------------------------------------------------------








#------------------------------------------------------------------------------------------------------
#-- get the growth data from ROS
#------------------------------------------------------------------------------------------------------
returnGrowthROS <- function(plotson=F){
  
  #- download the data
  #downloadHIEv(searchHIEv("ROS_MD_PM_GROWTH"),topath="Data")
  
  #- read it in, do some processing to harmonize the variables
  dat <- read.csv("Data/ROS_MD_PM_GROWTH_20111222-20130906_L1_edited.csv")
  dat$Date <- as.Date(dat$date,format="%Y-%m-%d")
  dat$date <- NULL
  dat$Species <- as.factor(ifelse(dat$sp == "casuari","cacu",
                                  ifelse(dat$sp=="radiata","pira",
                                         ifelse(dat$sp=="sideroxylon","eusi","eute"))))
  dat$sp <- NULL
  dat$Treat <- as.factor(ifelse(dat$trt=="water","wet","dry"))
  dat$trt <- NULL
  dat$diam <- rowMeans(dat[,5:6])
  #dat$ht <- as.numeric(as.character(dat$ht))
  dat <- dat[complete.cases(dat),]
  
  
  #- download the harvest files to build an allometry
  #files <- searchHIEv(filename="ROS_MD_PM_HARVEST_")
  #downloadHIEv(files,topath="Data")
  
  #- read in the harvest files. Each file has a different set of data, which is annoying.
  filestoread <- list.files(path="Data",pattern="ROS_MD_PM_HARVEST_",full.names=T)
  dat1 <- dat.subset <- list()
  for(i in 1:3){
    dat1[[i]] <- read.csv(filestoread[i])
    dat.subset[[i]] <- dat1[[i]][,c("date","sp","treenumber","rootsDM","stemDM","branchDM","leafDM","ht","diam.1","diam.2","nrleaves",
                                    "tenleafarea.1","tenleafarea.2")]
  }
  dat.allom <- do.call(rbind,dat.subset)
  
  #----------------------------------------------------------------------------------------
  #- the fourth harvest was different in how leaf area was estimated.
  dat4 <- read.csv(filestoread[4])
  
  dat4.measured <- subset(dat4,is.na(leafsubsampleDM)==T)
  dat4.est <- subset(dat4,is.na(leafsubsampleDM)==F)
  dat4.est$SLA <- with(dat4.est,(tenleafarea.1+tenleafarea.2)/2/leafsubsampleDM) # in cm/g
  dat4.est$leafArea <- with(dat4.est,(leafDM+leafsubsampleDM)*SLA/10000)
  
  dat4.measure2 <- subset(dat4.measured,select=c("date","sp","treenumber","rootsDM","stemDM","branchDM","leafDM","ht","diam.1","diam.2","nrleaves",
                                                 "tenleafarea.1","tenleafarea.2"))
  dat4.est2 <- subset(dat4.est,select=c("date","sp","treenumber","rootsDM","stemDM","branchDM","leafDM","ht","diam.1","diam.2","nrleaves",
                                        "tenleafarea.1","tenleafarea.2","leafArea"))
  #---------------------------------------------------------------------------------------- 
  
  dat.allom <- rbind(dat.allom,dat4.measure2)
  #- process allometry dataset
  dat.allom$leafArea <- rowMeans(dat.allom[,c("tenleafarea.1","tenleafarea.2")]/10000)
  
  dat.allom <- rbind(dat.allom,dat4.est2)
  
  dat.allom$diam <- rowMeans(dat.allom[,c("diam.1","diam.2")])
  dat.allom$d2h <- with(dat.allom,(diam/10)^2*ht)
  dat.allom$totDM <- rowSums(dat.allom[,c("rootsDM","stemDM","branchDM","leafDM")])
  dat.allom$logd2h <- log(dat.allom$d2h)
  dat.allom$logtotDM <- log(dat.allom$totDM)
  dat.allom$logLA <- log(dat.allom$leafArea)
  
  linkdf <- data.frame("allom.species" = levels(dat.allom$sp),"Species"= levels(dat$Species))
  dat.allom <- merge(dat.allom,linkdf,by.x="sp",by.y="allom.species")
  
  #------------------------------------------
  #-- ancova for mass
  #- plot 
  if(plotson==T) plotBy(logtotDM~logd2h|Species,data=dat.allom)
  
  #- ancova
  lm.full <- lm(logtotDM~logd2h*Species,data=dat.allom)
  lm2 <- lm(logtotDM~logd2h+Species,data=dat.allom)
  #anova(lm.full,lm2)
  if (plotson==1) anova(lm2) # common slope, different intercept
  
  #- predict total mass, add to dataframe with size
  dat$d2h <- with(dat,(diam/10)^2*ht)
  dat$logd2h <- log(dat$d2h)
  dat$totDM <- exp(predict(lm2,newdata=data.frame("logd2h" = dat$logd2h,"Species" = dat$Species)))
  
  #------------------------------------------
  
  
  
  
  #------------------------------------------
  #- ancova for leaf area
  if(plotson==T) windows();plotBy(logLA~logd2h|Species,data=dat.allom)
  
  #- ancova
  lm.full.la <- lm(logLA~logd2h*Species,data=dat.allom)
  lm2.la <- lm(logLA~logd2h+Species,data=dat.allom)
  #anova(lm.full.la,lm2.la) # use simpler model
  if (plotson==1) anova(lm2.la) # common slope, different intercept
  
  #- predict total mass, add to dataframe with size
  dat$d2h <- with(dat,(diam/10)^2*ht)
  dat$logd2h <- log(dat$d2h)
  dat$totLA <- exp(predict(lm2.la,newdata=data.frame("logd2h" = dat$logd2h,"Species" = dat$Species)))
  
  dat2 <- subset(dat,plottype=="shelter")
  
  
  #------------------------------------------
  
  #- get AGR and RGR of each timepoint for each plant
  growth <- dat2[with(dat2,order(treenumber,Date)),]
  growth.list <- split(growth,growth$treenumber)
  for (i in 1:length(growth.list)){
    growth.list[[i]]$AG <- c(0,diff(growth.list[[i]]$totDM))
    growth.list[[i]]$AGR <- c(0,diff(growth.list[[i]]$totDM)/as.numeric(diff(growth.list[[i]]$Date)))
    growth.list[[i]]$RGR <- c(0,diff(log(growth.list[[i]]$totDM))/as.numeric(diff(growth.list[[i]]$Date)))
  }
  growth2 <- do.call(rbind,growth.list)
  growth2$Treat <- relevel(growth2$Treat,ref="wet")
  return(growth2)
}
#------------------------------------------------------------------------------------------------------






#---------------------------------------------------------------------------------------------------------------
#--- read in the soil moisture release curve data
#--- note that the pressure is in bars, and the moisture is in a gravimetric percentage
getMoistCurve <- function(){
  curve1.1 <- read.csv("Data/VWC/ROS_soil_moisture_release_curve.csv")
  curve1 <- subset(curve1.1,Sample.No.=="ROS1")[,c("actual_pressure","moisture")]
  
  #- read in the 20bar data that Burhan sent on 15/05/2015.
  curve2 <- read.csv("Data/VWC/ROS_soil_moisture_release_curve_20bar_190515.csv")[1:4,c("actual_pressure","moisture")]
  
  #- put dataframes together
  curve <- rbind(curve1,curve2)
  
  curve$pressure_MPa <- curve$actual_pressure/10
  bd <- 1.32 #bulk density is 1.32 g cm-3. See bulk density excel file in /Data/VWC/
  curve$VWC <- curve$moisture*bd/100
  
  #plot(pressure_MPa~VWC,data=curve)
  return(curve)
}
#---------------------------------------------------------------------------------------------------------------








#---------------------------------------------------------------------------------------------------------------
#- functions to add lines to a plot
#---------------------------------------------------------------------------------------------------------------
addpoly <- function(x,y1,y2,col=alpha("lightgrey",0.7),...){
  ii <- order(x)
  y1 <- y1[ii]
  y2 <- y2[ii]
  x <- x[ii]
  polygon(c(x,rev(x)), c(y1, rev(y2)), col=col, border=NA,...)
}

predline <- function(fit, from=NULL, to=NULL, col=alpha("lightgrey",0.7), ...){
  
  if(is.null(from))from <- min(fit$model[,2], na.rm=TRUE)
  if(is.null(to))to <- max(fit$model[,2], na.rm=TRUE)
  
  newdat <- data.frame(X = seq(from,to, length=101))
  names(newdat)[1] <- names(coef(fit))[2]
  
  pred <- as.data.frame(predict(fit, newdat, se.fit=TRUE, interval="confidence")$fit)
  
  addpoly(newdat[[1]], pred$lwr, pred$upr, col=col)
  
  ablinepiece(fit, from=from, to=to, ...)
  
}

#'@title Add a line to a plot
#'@description As \code{abline}, but with \code{from} and \code{to} arguments. 
#'If a fitted linear regression model is used as asn argument, it uses the min and max values of the data used to fit the model.
#'@param a Intercept (optional)
#'@param b Slope (optional)
#'@param reg A fitted linear regression model (output of \code{\link{lm}}).
#'@param from Draw from this X value
#'@param to Draw to this x value
#'@param \dots Further parameters passed to \code{\link{segments}}
#'@export
ablinepiece <- function(a=NULL,b=NULL,reg=NULL,from=NULL,to=NULL,...){
  
  # Borrowed from abline
  if (!is.null(reg)) a <- reg
  
  if (!is.null(a) && is.list(a)) {
    temp <- as.vector(coefficients(a))
    from <- min(a$model[,2], na.rm=TRUE)
    to <- max(a$model[,2], na.rm=TRUE)
    
    if (length(temp) == 1) {
      a <- 0
      b <- temp
    }
    else {
      a <- temp[1]
      b <- temp[2]
    }
  }
  
  segments(x0=from,x1=to,
           y0=a+from*b,y1=a+to*b,...)
  
}
#---------------------------------------------------------------------------------------------------------------




#-------------------------------------------------------------------------------------------------------
#- function to return beta values of a thing (g1 or apparent Vcmax). Range from 0 to 1.
#  Accepts (1) a nls model based on a beta function, (2) a vector of values to predict from,
#   and (3) a type variable refering to the predictor (TDR or LWP)
#-------------------------------------------------------------------------------------------------------
predictBeta <- function(model,data,type="TDR"){
  params <- coef(model)
  
  #- find the data lower than Xlow (to return 0)
  lows <- which(data < unname(params[1]))
  
  #- find the data higher than Xhigh (to return 1.0)
  highs <- which(data > unname(params[2]))
  
  #- predict output value for all x-values
  if(type=="TDR") yvals <- predict(model,newdat=data.frame(TDR=data))
  if(type=="LWP") yvals <- predict(model,newdat=data.frame(LWP=data))
  
  #- override the low and high values in the predicted value
  yvals[lows] <- 0.1
  yvals[highs] <- 1
  
  return(yvals)
}
#-------------------------------------------------------------------------------------------------------

