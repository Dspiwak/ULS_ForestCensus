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

LiDAR_structure_check("/Processed_Data")
preprocessed_las<-list.files(path=paste0(wd,"/Raw_Data"),pattern=".laz")

if(length(preprocessed_las)>1){ #checks if only 1 .las file or many, currently only implemented for 1 scan file
  print("Too many .laz files present in raw data directory, consider revising (Only applicable following final code)")
}else(lasdf<-las_process(("Raw_Data"),1,2)) #generates a singular slice for this (1-2m) height


temp<-split_BufferedPointCloud(pc.dt=lasdf,plot.width=myproj@chunk_options$size,buffer.width=2) #chunks into 1ha grids with 2m buffer

chunk_size=sqrt(5000) #current chunk size is .5ha -or- 1 acre (~size of initial testing area)
grdpnts<-makegrid(testarea,cellsize=chunk_size) #makes a grid over the given shapefile but only contains "point" centers of the grids
spgrd<-SpatialPoints(grdpnts,proj4string = NAD83_2011) #converts into spatial points
spgrd<-SpatialPixels(spgrd) 
spgrd<-as(spgrd,"SpatialPolygons") #converts into "raster" chunks but each chunk is its own polygon which can be manipulated / focused on
Chull_Points<-list()
for(i in 1:length(spgrd)){
  grid_cell_points<-data.frame(spgrd[i]@polygons[[1]]@Polygons[[1]]@coords) #extracts border points about each chunk as coordinate dataframe
  coordinates(grid_cell_points)<-c('x','y') 
  tempmat<-cbind(grid_cell_points$x,grid_cell_points$y) #creates a temporary matrix for each chunk and contains grid coordinates
  chunkpnts<-chull(tempmat[,1],tempmat[,2]) #detects the convex hull points about the matrix of points and allows for 'reverse' engineering (may not need?)
  
  temp_points<-grid_cell_points[c(chunkpnts,chunkpnts[1]),] #remaps the points as spatial points
  
  Chull_Points[[i]]<-(temp_points) #could clean potentially and remap to temp_points?
  Chull_Points[[i]]<-as(Chull_Points[[i]],"SpatialLines") # creates a series of spatial lines which could be buffered along  to detect closeby clusters to merge
}




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
