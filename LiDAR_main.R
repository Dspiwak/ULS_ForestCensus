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
#clipregion<-testarea #This is a small .5ha 'test area' which the bulk of this processing was built apon and later expanded to the fullsize
clipregion<-Fullscan
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

#####
##### Stem Reconstructions ( Convexhull Fitting and Circle Fitting ) ----------

#Single Stem Fullslice Fitting for demonstration
Fitted_CirclePratt<-fit_circle(stems[1,],method='Pratt') 
Fitted_CircleLM<-fit_circle(stems[1,],method='LM') 
Fitted_Ransac_Circles<-fit_RANSAC(stems[90,],thresh = .5,.95,.5)

#Multislice Fitting (Needs "vertical profile ggplot)
stems$ConvH<-multislice_process(stems,fit_convhull) #Convex Hull approach

stems$Circle_Pratt<-multislice_process(stems,fit_circle,"Pratt")#Pratt Circle Fit

stems$Circle_LM<-multislice_process(stems,fit_circle,"LM")#Reduced Levenberg-Marquardt Method

#RANSAC Fitting
stems$Circle_RANSAC<-multislice_process(stems,fit_RANSAC,thresh=.5,prob=.95,w=.5)


#####
##### Matching & Producing Product to Assess Accuracy ----------
source("Scripts/Stem_matching.R")

#Automated Matching
matched_stems<-match_stems(stems,cliptruth,2.5)#double check this is using the correct ground truth (without dead stems)

#Convert to traditional SPDF and Export for Manual Matching
stems_spdf<-conv2_spdf(stems)
if(!file.exists("Processed_Data/Exported_Stem_Data.shp")){
  writeOGR(stems_spdf,"Processed_Data","Exported_Stem_Data","ESRI Shapefile",)
}
if(!file.exists("Processed_Data/Exported_GTruth_Data.shp")){
  Buffered_truth<-gBuffer(cliptruth,byid=TRUE,width=cliptruth$rbh_m)
  writeOGR(Buffered_truth,"Processed_Data","Exported_GTruth_Data","ESRI Shapefile")
}

#####
##### Stat. Analyses
source("Scripts/Stats_Outputs.R")

Matched_Stems<-read.csv("Processed_Data/Manually_Matched_Trees_Dataset_AllStems.csv")%>% #If manually matched
  filter(Confidence>1)
Matched_Stems<-data.frame(matched_stems)#If automatically matched

#Methods to check
measures<-list("ConvH","Circle_Pratt","Circle_LM","Circle_RANSAC")

#General Statistics (MAE, RMSE, Pearson's Correlation, Total Detection, etc. )
stats<-basic_stats(cliptruth,stems_spdf,Matched_Stems,measures)

#General Regressions
plot_regression(Matched_Stems$Gtruth_dbh_cm,Matched_Stems$ConvH)
plot_regression(Matched_Stems$Gtruth_dbh_cm,Matched_Stems$Circle_Pratt)
plot_regression(Matched_Stems$Gtruth_dbh_cm,Matched_Stems$Circle_LM)
plot_regression(Matched_Stems$Gtruth_dbh_cm,Matched_Stems$Circle_RANSAC)

#General Statistics

statdf<-data.frame(DBH=Matched_Stems$Gtruth_dbh_cm,
                   ConvH=Matched_Stems$ConvH,
                   Pratt=Matched_Stems$Circle_Pratt,
                   LM=Matched_Stems$Circle_LM,
                   RANSAC=Matched_Stems$Circle_RANSAC)%>%
  mutate(bin=cut_width(DBH,width=10,boundary=0))

#BIAS Boxplots

Biases<-calc_bias(statdf,10,plot=TRUE)

#MAE Plots
MAEs<-calc_MAE(statdf,plot=TRUE)




#extra binned count histograms (not working yet)
temphistdat<-as.data.frame(matched_stems)
bucket<-list(chull_dbh=temphistdat$chull_dbh_cm,truth_dbh=temphistdat$dbh_cm)
p2<-ggplot(melt(bucket),aes(value,fill=L1))+
  geom_histogram(position="stack",binwidth=10)
p2

tab<-table(temphistdat$dbh_cm,temphistdat$chull_dbh_cm)
barplot(tab)
