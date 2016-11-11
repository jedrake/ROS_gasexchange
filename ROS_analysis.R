#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
#- This is the central analysis script for a leaf-level gas exchange analysis of potted trees experiencing
#   controlled droughts under large outdoor rainout shelters (ROS). Most of the actual data manipulation occurs
#   in functions in dataFunction.R and plotFunctions.R. These functions are called by this script.
#   The intention is to keep this script usefully simple, but still recreate all of the tables and figures.
#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------------------
#- load the libraries and scripts that do all of the actual work.
source("R/loadLibraries.R")
#------------------------------------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------------------------------------
#- Download the data (from HIEv). This creates a "Data" and "Output" 
#  directories in the working directory, if they do not exist.
source("R/downloadData.R")
#------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------
#- plot soil moisture over time (Figure 1).
plotVWC(ptsize=1.8,output=T,type="1panel") # can plot type="1panel" or "4panel"
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------  
#- Plot photosynthetic variables relative to soil water content (Figure 2)
plotGX_theta(output=T,xlims=c(0,0.33),colors= brewer.pal(5,"Accent")[c(1,2,3,5)])#c("yellow3","blue","forestgreen","red"))

#- make a big 16-panel plot of photosynthetic variables over time (Figure S2).
plotGX(output=T)
#-------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------
#- plot the C isotope composition of ROS leaves, relate to gas exchange (Figure 3)
plotd13C_gx(output=T)
#-------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------
#- plot the soil moisture release curve (Figure S1).
pars <- plotMoistCurve(output=T)
#-------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------
#- plot the leaf water potential data over time (Figure S3).
plotLWP(fillcol="lightgrey",size=1.75,output=T,labsize=1.8)
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
#- plot normalized g1 and non-stomatal limitation as a function of VWC with beta functions.
#  Figure 4 and S4.

#- get the data, fit the beta functions. For some reason this cannot be dropped into a function.
#    Keep it as a script!
source("R/fitBeta_g1_nsl.R")
plotBetasG1NSL(output=T,g1data=g1values,NSLdata=NSLpars,g1list=g1_TDR_beta,NSLlist=NSL_TDR_beta)

source("R/fitBeta_g1_nsl_LWP.R")
plotBetasG1NSL_LWP(output=T,g1data=g1values2,NSLdata=NSLpars2,g1list=g1_TDR_beta2,NSLlist=NSL_TDR_beta2)

#- fit based on transpirable soil water
source("R/fitBeta_g1_nsl_MTSW.R")
plotBetasG1NSL.TSW(output=T,g1data=g1values,NSLdata=NSLpars,g1list=g1_TDR_beta.TSW,NSLlist=NSL_TDR_beta.TSW)
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
#- make a table of the beta parameters and standard errors. This table is written out as a csv in "Output".
#  After some manual manipulation, this becomes Table 1.
writeBetaParams()
#-------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------
#- Fit and plot the Tuzet model.

#- NOTE I need to update the non-stomatal parameters that are hard-coded inside the wrapper function
#   in this script!
source("R/fitTuzets.R")
#-------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------
#- Model Asat's response to drought given changes in g1, photosynthetic capacity, or both (Figure 7).

source("R/model photo vs observed.R")

#-- this bit is the old code
# Assumes the "fitBeta" scripts have been run.
# modelAsatVWC(output=T,fit.spg1=fit.spg1,fit.spNSL=fit.spNSL)
#-------------------------------------------------------------------------------------------------------
  
  
