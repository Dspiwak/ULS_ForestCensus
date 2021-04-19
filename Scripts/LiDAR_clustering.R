generate_chunks<-function(shp,chunk_size){
 #chunk size of sqrt(5000) = .5ha or 1 acre

  shp_buffer<-gBuffer(shp,width=chunk_size) #creates buffer which is needed to ensure a cell is created along edge when offset
  
  grdpnts<-makegrid(shp_buffer,cellsize=chunk_size) #makes a grid over the given shapefile but only contains "point" centers of the grids
  
  spgrd<-SpatialPoints(grdpnts,proj4string = NAD83_2011) %>% #converts into spatial points
    SpatialPixels()%>% #generate chunks &clip
    as("SpatialPolygons")
  offset_spgrd<-elide(spgrd,shift=c(chunk_size/2,chunk_size/2))
  
  crs(offset_spgrd)<-crs(spgrd)
  
  chunks<-rbind(spgrd,offset_spgrd,makeUniqueIDs=TRUE) %>%
    gIntersection(shp,byid=TRUE)
  
  plot(shp)
  plot(chunks[1:(length(chunks)/2)],col=NA,add=TRUE)
  plot(chunks[((length(chunks)/2)+1):length(chunks)],col=NA,border='red',add=TRUE)
  
  return(chunks)
}

resolve_stems<-function(chunked_polys){
  #Takes chunks (in form of a SP input) and convex hulls (or other reconstruction in form of SPDF) which have been 'chunked'
  #Will return a set of sliced convex hulls / SP which have been resolved due to their "splitting" between chunks.
  #This process recalculates the convex hull after UnaryUnion of the overlapping polys.
  
  
  mergedHulls<-gUnaryUnion(chunked_polys)%>% #merges and dissolves boundaries across all overlapping polygons. 
    disaggregate()
  
  SingleStems<-data.frame('PointsInStem'=matrix(NA,ncol=1,nrow=length(mergedHulls)))
  for(i in 1:length(mergedHulls)){
    SingleTree_Points<-intersect(points_df,mergedHulls[i])
    
    SingleTree_Coords<-data.frame(SingleTree_Points)
    
    SingleStems$PointsInStem[i]=list(SingleTree_Coords)
    SingleStems$StemID[[i]]=i
  }
  
  return(SingleStems)
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
  coordinates(points.df)<-c("X","Y","Z") #converts point data back to spatial such that convex hulls can be created later
  crs(points.df)<-NAD83_2011
  
  return(points.df)
}