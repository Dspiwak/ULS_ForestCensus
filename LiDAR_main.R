#LiDAR_main.R
#   "quickscript" for continuing thesis research
#   Reformatted LiDAR_TreeChulls.R 7/31/2020 into sub-scripts which will be easier to work with and tweak in future
# 
#This script will "call" all other scripts and check to ensure all files/ functions are present.

##### Set Working Directory and Initial Setup of Packages and CRS Projections----------
wd<-dirname(rstudioapi::getSourceEditorContext()$path) #Should be located where all other files / functions are on system
setwd(wd)

source("Scripts/Parent_structure_check.R")
source("Scripts/Setup.R")

setup()

AllFile_CRS<-6346 #NAD83_2011 UTM zone 17N epsg code
NAD83_2011<-crs(paste0("+init=epsg:",AllFile_CRS))#creates CRS object with the Allfile epsg code
SetCRS<-crs("+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs ") #this is the projection used by the Smithsonian for their datasets


#####
##### Preform Forest Geo Processing ----------
source("Scripts/FGeo_Processing.R")


#####
##### LiDAR Processing ( Chunking / Slicing / Clustering ) ----------
source("Scripts/LiDAR_structure_check.R")
source("Scripts/LiDAR_initial_processing.R")
source("Scripts/LiDAR_clustering.R")
source("Scripts/Stem_reconstructions.R")

LiDAR_structure_check("/Processed_Data")
preprocessed_las<-list.files(path=paste0(wd,"/Raw_Data"),pattern=".laz")

if(length(preprocessed_las)>1){ #checks if only 1 .las file or many, currently only implemented for 1 scan file
  print("Too many .laz files present in raw data directory, consider revising (Only applicable following final code)")
}else(lasdf<-las_process(("Raw_Data"),1,2)) #generates a singular slice for this (1-2m) height

#To adjust clipping / analysis region, set it in the forest geo processing section (for now)


points_df<-lasdf
coordinates(points_df)<-c("X","Y") #convert dataframe into spatial points data type
chunks<-generate_chunks(Fullscan)

offchunk<-generate_offset_chunks(Fullscan)





loop_hulls<-list()
for(i in 1:length(chunks)){
  chunk<-chunks[i]
  clipped_lasdf<-intersect(points_df,chunks)
  clipped_lasdf<-data.frame(clipped_lasdf)
  clustered_points<-cluster_lidar_dbscan(clipped_lasdf,20,0.0067)
  loop_hulls[i]<-conv_hulls(clustered_points)
}

temp<-loop_hulls[[1]]
for(i in 2:length(loop_hulls)){
temp<-rbind(temp,loop_hulls[[i]],makeUniqueIDs=TRUE)
}
writeOGR(temp,dsn=paste0(wd,'/Processed_Data'),layer=paste0(deparse(substitute(temp)),"_AutomatedOuput"),driver="ESRI Shapefile",overwrite_layer=TRUE)

#####
##### LiDAR Processing ( Convex Hulls / Fourier Transform / Matching) ----------
source("Scripts/Stem_matching.R")
matched_stems<-match_stems(temp,2.5)


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
