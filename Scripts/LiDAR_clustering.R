generate_chunks<-function(lasdf,extent_shp){
  chunk_size=sqrt(5000)
  grdpnts<-makegrid(extent_shp,chunk_size)
  spgrd<-SpatialPoints(grdpnts,proj4string = NAD83_2011)
  spgrdWithin<-SpatialPixels(spgrd[extent_shp])
  spgrdWithin<-as(spgrdWithin,"SpatialGrid")
  #spgrdWithin<-as(spgrdWithin,"SpatialPolygons")
  
  return(spgrdWithin)
}

get.elbow.points.indices<-function(x,y,threshold){
  d1<-diff(y)/diff(x) #first derivative
  d2<-diff(d1)/diff(x[-1]) #second derivative
  indices<-which(abs(d2)>threshold)               #### Maybe change to < ? to detect the max curve from bottom? ####
  return(indices)
}

auto_EPS<-function(df,minpnts,threshold){
  temp<-kNNdist(df,k=minpnts,all=TRUE) #distances of the knn of any selected point within dataset
  index<-as.numeric(rownames(data.frame(temp)))
  val<-rowMeans(temp)
  temp<-data.frame("x"=index,"y"=val)%>% #creates a temporary table containing the average knn distances for each point
    arrange(desc(y))%>%
    mutate(index=row_number(y))
  indices<-get.elbow.points.indices(temp$index,temp$y,threshold) #detects the "knee" of the data based on 2nd derivative
  crit_points<<-data.frame("x"=temp$index[indices],"y"=temp$y[indices])
  p<<-ggplot(data=temp,aes(x=index,y=y))+ #allows the data to be plotable
    geom_point()+
    geom_point(data=crit_points,aes(x=x,y=y),color="red")+
    geom_hline(aes(yintercept=min(crit_points$y)),color="blue",linetype="dashed")+
    geom_text(aes(0,min(crit_points$y),label=min(crit_points$y),vjust=-.6),hjust=-2)
  return(min(crit_points$y))
}

cluster_lidar_dbscan<-function(las_dataframe,minpnts,KnnThreshold){
  #DBSCAN clustering of LiDAR points
  #minpnts=minimum # of points per "tree"
  #KnnThreshold is the sensitivity of the automatic calculation of epsilon (may need refinement)
  
  clusters_db<-dbscan(las_dataframe,eps=auto_EPS(las_dataframe,minpnts,threshold=KnnThreshold),minPts=minpnts) #clusters based on params
  las_dataframe$id<-clusters_db$cluster #creates a new column under the input dataframe to record the cluster #id
  
  points.df<-las_dataframe[,c("id","X","Y","Z")] %>% 
    filter(id>0) #removes the "noise" (denoted with id=0) from the dataset
  coordinates(points.df)<-c("X","Y") #converts point data back to spatial such that convex hulls can be created later
  crs(points.df)<-NAD83_2011
  
  return(points.df)
}