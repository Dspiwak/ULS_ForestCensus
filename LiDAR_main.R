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

#Display "knee in curve" graph example for a chunk
templasdf<-intersect(points_df,chunks_with_points[4])%>%
  data.frame()%>%
  auto_EPS(20,0.0067,plot = TRUE)


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
Fitted_ConvHull<-fit_convhull(stems[886,],plot=TRUE)
Fitted_CirclePratt<-fit_circle(stems[886,],method='Pratt',plot=TRUE) 
Fitted_CircleLM<-fit_circle(stems[886,],method='LM',plot=TRUE) 
Fitted_Ransac_Circles<-fit_RANSAC(stems[886,],thresh = .5,.95,.5,plot=TRUE)

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

Matched_Stems<-read.csv("Processed_Data/Manually_Matched_Trees_Dataset_AllStems_run2.csv")%>% #If manually matched
  filter(Confidence>1)
#Matched_Stems<-data.frame(matched_stems)#If automatically matched

#Methods to check
measures<-list("ConvH","Circle_Pratt","Circle_LM","Circle_RANSAC")

#General Statistics (MAE, RMSE, Pearson's Correlation, Total Detection, etc. )
stats<-basic_stats(Clipped.GTruth,stems_spdf,Matched_Stems,measures)

#Remove Excesive noise >1.5x the max known gtruth
Matched_Stems<-Matched_Stems%>%
  filter_at(vars(!!!syms(measures)),~.<max(Clipped.GTruth$dbh_cm)*1.5)#filters at all columns in the 'measures' list

#General Regressions
p1<-plot_regression(Matched_Stems,"Gtruth_dbh_cm","ConvH")
p2<-plot_regression(Matched_Stems,"Gtruth_dbh_cm","Circle_Pratt")
p3<-plot_regression(Matched_Stems,"Gtruth_dbh_cm","Circle_LM")
p4<-plot_regression(Matched_Stems,"Gtruth_dbh_cm","Circle_RANSAC")
wrap_plots(p1,p2,p3,p4)
#General Statistics

statdf<-data.frame(DBH=Matched_Stems$Gtruth_dbh_cm,
                   ConvH=Matched_Stems$ConvH,
                   Pratt=Matched_Stems$Circle_Pratt,
                   LM=Matched_Stems$Circle_LM,
                   RANSAC=Matched_Stems$Circle_RANSAC)%>%
  mutate(bin=cut_width(DBH,width=10,boundary=0))

#Detection Plots
plot_detected(Clipped.GTruth,Matched_Stems,width=10,plot_percentages=TRUE)

#BIAS Boxplots
Biases<-calc_bias(statdf,10,RemoveBinsUnder=10,plot=TRUE,plot_stacked=TRUE)

#MAE Plot
MAEs<-calc_MAE(statdf,plot=TRUE)

#MAPE Plot
MAPEs<-calc_MAPE(statdf,plot=TRUE)


#TestArea Results (Recalc / reclip the matched stems to only include those in the test area)
TestArea_Truthdf<-cliptruth%>%
  intersect(testarea)%>%
  data.frame()
TestArea_MatchedStems<-Matched_Stems%>%
  dplyr::filter(Gtruth_tag%in%TestArea_Truthdf$tag)%>%
  dplyr::filter(Confidence>1)

TestArea_stats<-basic_stats(cliptruth%>%intersect(testarea),stems_spdf%>%intersect(testarea),TestArea_MatchedStems,measures)

pTA1<-plot_regression(TestArea_MatchedStems,"Gtruth_dbh_cm","ConvH")
pTA2<-plot_regression(TestArea_MatchedStems,"Gtruth_dbh_cm","Circle_Pratt")
pTA3<-plot_regression(TestArea_MatchedStems,"Gtruth_dbh_cm","Circle_LM")
pTA4<-plot_regression(TestArea_MatchedStems,"Gtruth_dbh_cm","Circle_RANSAC")
wrap_plots(pTA1,pTA2,pTA3,pTA4)

statdf2<-data.frame(DBH=TestArea_MatchedStems$Gtruth_dbh_cm,
                    ConvH=TestArea_MatchedStems$ConvH,
                    Pratt=TestArea_MatchedStems$Circle_Pratt,
                    LM=TestArea_MatchedStems$Circle_LM,
                    RANSAC=TestArea_MatchedStems$Circle_RANSAC)%>%
  mutate(bin=cut_width(DBH,width=10,boundary=0))

plot_detected(cliptruth%>%intersect(testarea),TestArea_MatchedStems,width=10,plot_percentages=TRUE)
calc_bias(statdf2,10,plot=TRUE,plot_stacked=TRUE)
calc_MAE(statdf2,plot=TRUE)
calc_MAPE(statdf2,plot=TRUE)
