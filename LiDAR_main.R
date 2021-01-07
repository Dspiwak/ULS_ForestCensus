#LiDAR_main.R
#   "quickscript" for continuing thesis research
#   Reformatted LiDAR_TreeChulls.R 7/31/2020 into sub-scripts which will be easier to work with and tweak in future
# 
#This script will "call" all other scripts and check to ensure all files/ functions are present.

##### Set Working Directory and Local CRS's----------
wd<-"C:/Users/note2/Desktop/thesis2021/" #Should be located where all other files / functions are on system
setwd(wd)


#####
##### Initialize All Packages and Local CRS----------
library(raster)
library(plyr)
library(dplyr)
library(rgeos)
library(rgdal)
library(lidR)
library(dbscan)
library(grDevices)
library(sp)
library(approximator)
library(mapview)
library(ggplot2)
library(ggpubr)
library(grid)
library(MeanShiftR) #this package requires you to install devtools to download off of github, but is not needed for you colin
library(future)

AllFile_CRS<-6346 #NAD83_2011 UTM zone 17N epsg code
NAD83_2011<-crs(paste0("+init=epsg:",AllFile_CRS))#creates CRS object with the Allfile epsg code
SetCRS<-crs("+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs ") #this is the projection used by the Smithsonian for their datasets


#####
##### Check Folder Structures / Etc. ----------
source("Scripts/Parent_structure_check.R")

#####
##### Preform Forest Geo Processing ----------
source("Scripts/FGeo_Processing.R")


#####
##### LiDAR Processing ( Chunking / Slicing / Clustering ) ----------
source("Scripts/LiDAR_structure_check.R")
source("Scripts/LiDAR_initial_processing.R")
source("Scripts/LiDAR_clustering.R")

LiDAR_structure_check("Processed_Data")
preprocessed_las<-list.files(path=paste0(wd,"Raw_Data"),pattern=".laz")

if(length(preprocessed_las)>1){ #checks if only 1 .las file or many, currently only implemented for 1 scan file
  print("Too many .laz files present in raw data directory, consider revising (Only applicable following final code)")
}else(lasdf<-las_process(("Raw_Data"),1,2)) #generates a singular slice for this (1-2m) height

clustered_points<-cluster_lidar_dbscan(lasdf,20,0.0067)


#####
##### LiDAR Processing ( Convex Hulls / Fourier Transform / Matching) ----------
source("Scripts/Stem_reconstructions.R")
source("Scripts/Stem_matching.R")

stem_hulls<-conv_hulls(clustered_points)
matched_stems<-match_stems(stem_hulls,2.5)


#####
##### Stat. Analyses

rmse<-sqrt(mean((matched_stems$chull_dbh_cm-matched_stems$dbh_cm)^2))
lin_model<-lm(matched_stems$chull_dbh_cm~matched_stems$dbh_cm)
cor(matched_stems$chull_dbh_cm,matched_stems$dbh_cm,method="pearson")

summary(lin_model)
cat(paste0("With an RMSE of ",rmse," cm",
           "\n Searching in tree sizes >",Smallest.Tree.Size,
           "\n Total trees present:",nrow(Clipped.GTruth),
           "\n Total trees detected:",nrow(stem_hulls)))

p<-ggplot(data=as.data.frame(matched_stems),aes(x=dbh_cm,y=chull_dbh_cm,ymin=0,ymax=max(chull_dbh_cm),xmin=0,xmax=max(dbh_cm)))+
  geom_smooth(method="lm",formula=y~x)+
  geom_point()+
  coord_equal()+
  stat_cor(label.y=max(matched_stems$chull_dbh_cm))+
  stat_regline_equation(label.y=max(matched_stems$chull_dbh_cm)-5)
p
plot(lin_model)
