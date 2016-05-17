#-----------------------------------------------------------------------------------------
#- Libraries and scripts required to run stuff in the ROS repo
#-----------------------------------------------------------------------------------------
library(doBy)
library(plotBy)
library(HIEv)
library(stringr)
library(reshape)
library(plyr)
library(dplyr)
library(plantecophys)
library(stringr)
library(magicaxis)
library(colorRamps)
library(car)
library(scales)
library(DEoptim)

source("R/dataFunctions.R")
source("R/plotFunctions.R")
setToken(tokenfile="HIEv_token.txt")
