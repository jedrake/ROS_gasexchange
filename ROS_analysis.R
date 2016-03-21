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
plotGX(output=F)
  