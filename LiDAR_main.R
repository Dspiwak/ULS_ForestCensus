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
clipregion<-testarea
cliptruth<-intersect(Ground.Truth,clipregion)

if(length(preprocessed_las)>1){ #checks if only 1 .las file or many, currently only implemented for 1 scan file
  print("Too many .laz files present in raw data directory, consider revising (Only applicable following final code)")
}else(lasdf<-las_process(("Raw_Data"),1,4)) #generates a singular slice for this (1-2m) height

#To adjust clipping / analysis region, set it in the forest geo processing section (for now)


points_df<-lasdf
coordinates(points_df)<-c("X","Y","Z") #convert dataframe into spatial points data type
chunks<-generate_chunks(clipregion,chunk_size=sqrt(5000))


chunks_with_points<-intersect(chunks,points_df) #extract only the chunks that have points in them...primarily useful for testing

#Preform initial runthrough of clustering and convex hull generation of points.
minPnts<-20
KnnThreshold<-0.0067

for(i in 1:length(chunks_with_points)){
  chunk<-chunks_with_points[i]
  clipped_lasdf<-intersect(points_df,chunk) %>%
    data.frame()
  
  clustered_points<-cluster_lidar_dbscan(clipped_lasdf,minPnts,KnnThreshold)
  
  current_chunk_hulls<-conv_hulls(clustered_points,plot=FALSE)
  current_chunk_hulls$chunk_hull_ID<-paste(i,'-',current_chunk_hulls$ClusterID)
  
  if(i==1){
    chunked_hulls<-current_chunk_hulls
  }
  else(chunked_hulls<-rbind(chunked_hulls,current_chunk_hulls,makeUniqueIDs=TRUE))
}

stems<-resolve_stems(chunked_hulls) #will "resolve" chunking of convex hulls...but rewrites / recalcs the dbh data

#Fullslice Circle Fitting
Fitted_CirclePratt<-fit_circle(stems,method='Pratt',plot=TRUE) 

Fitted_CircleLM<-fit_circle(stems,method='LM',plot=TRUE) 


#Multislice Fitting (Needs "vertical profile ggplot)
stems$ConvH<-multislice_process(stems,fit_convhull) #Convex Hull approach

stems$Circle_Pratt<-multislice_process(stems,fit_circle,"Pratt")#Pratt Circle Fit

stems$Circle_LM<-multislice_process(stems,fit_circle,"LM")#Reduced Levenberg-Marquardt Method


#RANSAC Fitting
Fitted_Ransac_Circles<-RANSAC_circle(points_df,stems,.95,.5,3)

#Fourier Curve Fitting
Adjusted_CoordData<-conv2_polar(points_df,Fitted_CircleData,stems)
for(i in 1:8){
  dat<-Adjusted_CoordData[[1]][[i]]
  res<-fit_fourier(dat,n=8,up=10,plot=TRUE) #maximum value of harmonics is nrow(dat)-2 (wang et all set this to 8 / 3 depending on AOI)
}



#Testing fourier fitting 1
Test_FourierDat<-data.frame(dist=c(1,2,2,2,3,4,5,2,1),azimuth=(seq(-pi,pi,length.out=length(Test_FourierDat$dist))))
res<-fit_fourier(Test_FourierDat,n=3,up=10,plot=TRUE)

plot(Test_FourierDat$azimuth,Test_FourierDat$dist)
lines(res$time,res$y,col='red')
lines(Test_FourierDat$azimuth,blek,col='blue')



#Testing fourier fitting 2 (tree model)
#Test_FourierDatApp<-read.csv(paste0(getwd(),"/Processed_Data/TestFourierSeries")) #Challenge Tree with all cases
Test_FourierDatApp<-read.csv(paste0(getwd(),"/Processed_Data/TestTree_OptimalTree"))
temp<-chull(Test_FourierDatApp)
Test_FourierDatApp_Hull<-Test_FourierDatApp[c(temp,temp[1]),]
testchull<-SpatialPolygons(list(Polygons(list(Polygon(Test_FourierDatApp_Hull)),ID=1)))
testchull$id<-1
coordinates(Test_FourierDatApp)<-c("x","y")
mockdat<-data.frame(seq(1:length(Test_FourierDatApp)))
spdf<-SpatialPointsDataFrame(coords=Test_FourierDatApp,data=mockdat)
spdf$X<-spdf$x
spdf$Y<-spdf$y
#testcircdat<-data.frame("cx"=14.0,"cy"=10.06554,"radius"=7.205898) #original cx=5.434675 / cy=13.06554 (ONLY FOR CHALLANGE TREE)
testcircdat<-fit_circle(spdf,hulls = testchull)
#coordinates(testcircdat)<-c("cx","cy")
testpolar<-conv2_polar(spdf,testcircdat,testchull,plot=TRUE)
res2<-fit_fourier(testpolar[[1]][[1]],n=2,up=10,plot=TRUE)

testfit_in<-testpolar[[1]][[1]]$dist
testfit_out<-fft(testfit_in)
barplot(Mod(testfit_out[2:(length(testfit_in)/2+1)]),main='Supposed "Energy" of Harmonics')


#testmod<-lm(testpolar[[1]][[1]]$azimuth+testpolar[[1]][[1]]$dist~res2$time+res2$y)
testpolarcleaned<-remove_fourieroutliers(res2,testpolar[[1]][[1]])

calc_arclength(res2,testpolar[[1]][[1]])

temppolarplot<-ggplot()+
  geom_point(data=testpolar[[1]][[1]],aes(azimuth,dist))+
  coord_polar(theta="x",direction=-1,start=pi/2)+
  scale_y_continuous(limits=c(0,max(testpolar[[1]][[1]]$dist)+.2))+
  geom_line(data=res2,aes(time,y))+
  geom_point(data=testpolarcleaned,aes(azimuth,dist),color='red')
tempcartplot<-ggplot()+
  geom_point(data=testpolar[[1]][[1]],aes(azimuth,dist))+
  geom_line(data=res2,aes(time,y))+
  geom_point(data=testpolarcleaned,aes(azimuth,dist),color='red')+
  scale_y_continuous(limits=c(0,max(testpolar[[1]][[1]]$dist)+.05))
grid.arrange(tempcartplot,temppolarplot,ncol=2)

#write.csv(res2,file='C:/Users/note2/Documents/GitHub/ULS_ForestCensus/Processed_Data/Example_FourierCurve_Data.csv')
#writeOGR(hulls,dsn=getwd(),layer=paste0(deparse(substitute(hulls)),"__TestCheck__AutomatedOuput"),driver="ESRI Shapefile",overwrite_layer=TRUE)


plot(res$time,res$y)

# ggplot(data=AdjustedCoords,aes(azimuth,dist))+
#   geom_point()+
#   geom_line(data=res,aes(time,y))


#####
##### LiDAR Processing ( Convex Hulls / Fourier Transform / Matching) ----------
source("Scripts/Stem_matching.R")
#A temporary fix (will need to re-write the ID somehow)
hulls$clusid<-c(1:length(hulls))
matched_stems<-match_stems(stems,cliptruth,2.5)


#####
##### Stat. Analyses
source("Scripts/Stats_Outputs.R")

measures<-list("ConvH","Circle_Pratt","Circle_LM")
stats<-basic_stats(cliptruth,stems,matched_stems,measures)

plot_regression(matched_stems$dbh_cm,matched_stems$ConvH)
plot_regression(matched_stems$dbh_cm,matched_stems$Circle_Pratt)
plot_regression(matched_stems$dbh_cm,matched_stems$Circle_LM)

lin_model<-lm(matched_stems$ConvH~matched_stems$dbh_cm)
plot(lin_model)













temphistdat<-as.data.frame(matched_stems)
bucket<-list(chull_dbh=temphistdat$chull_dbh_cm,truth_dbh=temphistdat$dbh_cm)
p2<-ggplot(melt(bucket),aes(value,fill=L1))+
  geom_histogram(position="stack",binwidth=10)
p2

tab<-table(temphistdat$dbh_cm,temphistdat$chull_dbh_cm)
barplot(tab)
