setup<-function(){
  cat(" Checking package availability and installing required packages...",
      "\n     Note: Some packages require 'rtools' or external github calls to specific repositories \n")
  
  required_packages<-c("raster","plyr","dplyr","rgeos","rgdal","lidR","dbscan","grDevices","sp",
                      "approximator","mapview","ggplot2","ggpubr","grid","MeanShiftR","future")
  
  install.packages(setdiff(required_packages,rownames(installed.packages())))
  
  cat("\n Check complete, all required packages are installed.","\n Proceeding...")
  
  for(i in 1:length(required_packages)){
    library(required_packages[i],character.only = TRUE)
  }
}
