#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
#- This is the central analysis script for a leaf-level gas exchange analysis of potted trees experiencing
#   controlled droughts under large outdoor rainout shelters (ROS). Most of the actual data manipulation occurs
#   in functions in dataFunction.R and plotFunctions.R. These functions are called by this script.
#   The intention is to keep this script usefully simple, but still recreate all of the tables and figures.
#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

#- load the libraries and scripts that do all of the actual work.
source("R/loadLibraries.R")


#-------------------------------------------------------------------------------------------------------
#- plot soil moisture over time (Figure 1).
plotVWC(ptsize=1.8,output=T,type="1panel") # can plot type="1panel" or "4panel"
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
#- plot the leaf water potential data over time (Figure 2).
plotLWP(fillcol="lightgrey",size=1.75,output=T,labsize=1.8)
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------  
#- make a big 16-panel plot of photosynthetic variables over time (Figure 3).
plotGX(output=T)

#- OR, make an alternative plot of photosynthetic variables relative to soil water content (Figure 3)
plotGX_theta(output=T)
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
#- plot the soil moisture release curve (Figure S2).
pars <- plotMoistCurve(output=T)
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
#- plot Photo and Cond vs. leaf water potential, to show iso vs. anisohydry (Figure S3).
plotHydry(output=T)
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
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
#- make a table of the beta parameters and standard errors. This table is written out as a csv in "Output".
#  After some manual manipulation, this becomes Table 1.
writeBetaParams()
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
#- Model Asat's response to drought given changes in g1, photosynthetic capacity, or both (Figure 5).
# Assumes the "fitBeta" scripts have been run.
modelAsatVWC(output=T,fit.spg1=fit.spg1,fit.spNSL=fit.spNSL)
#-------------------------------------------------------------------------------------------------------
  
  
#-------------------------------------------------------------------------------------------------------
#- plot the C isotope composition of ROS leaves, relate to gas exchange (Figure 6)
plotd13C_gx(output=T)
#-------------------------------------------------------------------------------------------------------