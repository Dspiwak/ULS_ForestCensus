#Useful random functions (Testing area)

#####
#Minimum distanct to chunk boundary partial function
#Will 'chunk' a region, and allow for seperation of individual lines of the chunk...
#Cluster points will be checked/ intersected with buffered region and if falls inside... 
#append to list and recalculate convex hulls based on 'matched' / closest convex hulls that have been fragmented
chunk_size=sqrt(5000) #current chunk size is .5ha -or- 1 acre (~size of initial testing area)
grdpnts<-makegrid(Fullscan,cellsize=chunk_size) #makes a grid over the given shapefile but only contains "point" centers of the grids
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

plot(Fullscan)
plot(testarea,add=TRUE,col='red')
for(i in Chull_Points){
  plot(i,add=TRUE,col='blue')
}
#####


#clean the spatial polygons weirdly
tempy<-data.frame(ID=character(),stringsAsFactors = FALSE)
for(i in spgrd@polygons){
  tempy<-rbind(tempy,data.frame(ID=i@ID,stringsAsFactors = FALSE))
}
row.names(tempy)<-tempy$ID
spgrd<-SpatialPolygonsDataFrame(spgrd,tempy)

#####


clustered_points<-cluster_lidar_dbscan(lasdf,20,0.0067)
#####

