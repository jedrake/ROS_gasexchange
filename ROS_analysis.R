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

#- plot data and predicted values.
plotBetasG1NSL(output=F,g1data=g1values,NSLdata=NSLpars,g1list=g1_TDR_beta,NSLlist=NSL_TDR_beta)
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
#- plot the C isotope composition of ROS leaves
plotd13C(export=T)
#-------------------------------------------------------------------------------------------------------




      