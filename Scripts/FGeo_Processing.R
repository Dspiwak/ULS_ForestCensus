#FGeo_Processing.R
#   Reformatted version from LiDAR_TreeChulls.R 7/31/2020 which will preform initial FGeo processing
#Will convert local x/y coords to NAD83 lat/lon coords and remove non-dominante trees from census data

file<-as.data.frame(read.csv(paste0(getwd(),"/Raw_Data/scbi.stem3.csv"))) #Raw FGEO Census data

#SIGEO Local coords conversion (Adapted from SI script)
NAD83.SW<-c(747385.521,4308506.438)
NAD83.NW<-c(747370.676,4309146.156)
offset<-atan2(NAD83.NW[1]-NAD83.SW[1],NAD83.NW[2]-NAD83.SW[2])
grid2nad83<-function(x,y){
  NAD83.X<-NAD83.SW[1]+(x*cos(offset)+y*sin(offset))
  NAD83.Y<-NAD83.SW[2]+(-x*sin(offset)+y*cos(offset))
  nad83<-list(NAD83.X,NAD83.Y)
  names(nad83)<-c("NADX","NADY")
  nad83
}

## DF processing and conversion to spatial points
Smallest.Tree.Size<-18 #cm DBH
loclist<-grid2nad83(file$gx,file$gy)
Ground.Truth<-file
Ground.Truth$dbh[Ground.Truth$dbh==0]<-NA
Ground.Truth$dbh[Ground.Truth$dbh=="NULL"]<-NA#removes trees with missing / 0 dbh
Ground.Truth<-mutate(Ground.Truth,quadrat_x=as.numeric(as.character(file$quadrat))%/% 100) %>%
  mutate(quadrat_y=as.numeric(as.character(file$quadrat))%%100) %>%
  mutate(grid_x=loclist$NADX,grid_y=loclist$NADY,NADX=loclist$NADX,NADY=loclist$NADY)%>%
  mutate(dbh_cm=as.numeric(dbh)/10)%>%
  mutate(rbh_m=(dbh_cm/2)/100)%>%#create a radius field which is in m (for buffering points later)
  dplyr::select(-dbh,-ExactDate,-date,-MeasureID,-CensusID)%>%
  filter(dbh_cm>=Smallest.Tree.Size,!is.na(tag),!is.na(dbh_cm))%>% #filters for trees greater than 12cm dbh and ensures it has tag info
  filter(!grepl('B|P|X|DC|DT|DN',codes))%>%#ensure tree is not broken,prostrate(parallel to ground),dead fallen/missing
  filter(grepl('caco|cagl|caovl|cato|fagr|fram|litu|qual|qupr|quru',sp)) #ensure it is is a dominant species
coordinates(Ground.Truth)<-Ground.Truth[,c("NADX","NADY")] #converts from dataframe to spatial points data frame
crs(Ground.Truth)<-SetCRS #sets crs of the ground truth to that of the known deer exclusion zone/ shapefile
Ground.Truth<-spTransform(Ground.Truth,NAD83_2011)
Clipped.GTruth<-Ground.Truth
#Clipped.GTruth<-intersect(Ground.Truth,deerexclusion) #clips the full dataset down to the deer exclusion zone
Fullscan<-readOGR(paste0(getwd(),"/Raw_Data/FullScanExtent.shp")) #full extents of scan(~12ha.)
deerexclusion<-readOGR(paste0(getwd(),"/Raw_Data/deer_exclosure_2011.shp")) %>% #deer exclusion (~5ha.)
  spTransform(NAD83_2011)
testarea<-readOGR(paste0(getwd(),"/Raw_Data/testarea.shp")) #test region where I "tuned" values and collected manual digitization (.5ha)
clipregion<-Fullscan
Clipped.GTruth<-intersect(Clipped.GTruth,clipregion)
Stem_Buff<-gBuffer(Clipped.GTruth,byid=TRUE,width=Clipped.GTruth$rbh_m) # data is in mm but conv to cm for radius buff