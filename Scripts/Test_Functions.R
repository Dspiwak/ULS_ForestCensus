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


offset_grdpnts<-makegrid(shp_buffer,cellsize=chunk_size,offset=c(1,1))

library(maptools)

scanregion_buff<-gBuffer(Fullscan,width=sqrt(5000))

grdpnts<-makegrid(scanregion_buff,cellsize=sqrt(5000)) #makes a grid over the given shapefile but only contains "point" centers of the grids

spgrd<-SpatialPoints(grdpnts,proj4string = NAD83_2011) %>% #converts into spatial points
  SpatialPixels()%>% #generate chunks &clip
  as("SpatialPolygons")%>%
  gIntersection(Fullscan,byid=TRUE)

offset_grdpnts<-makegrid(scanregion_buff,cellsize=sqrt(5000))
offset_spgrd<-SpatialPoints(offset_grdpnts,proj4string = NAD83_2011) %>% #converts into spatial points
  SpatialPixels()%>%
  as("SpatialPolygons")
offset_spgrd<-elide(offset_spgrd,shift=c(sqrt(5000)/2,sqrt(5000)/2))%>%
  gIntersection(Fullscan,byid=TRUE)
  

chunks<-rbind(spgrd,offset_spgrd,makeUniqueIDs=TRUE)


#####

all_hulls<-chunked_hulls
all_hulls@data$intersects<-NA
all_hulls@data$mergedIDs<-NA

for(i in 1:length(all_hulls)){
  for(k in 1:length(all_hulls)){
    if(gIntersects(all_hulls[i,],all_hulls[k,])==TRUE && i!=k){
      all_hulls@data[i,]$intersects<-TRUE
      all_hulls@data[k,]$intersects<-TRUE
      temp_poly<-raster::bind(all_hulls[i,],all_hulls[k,])
      temp_poly@data$mergedIDs<-paste0(all_hulls@data[i,]$clusid,'-',all_hulls@data[k,]$clusid)
      #plot(temp_poly,col='grey')
      #plot(all_hulls[i,],add=TRUE,col=NA,border='red')
      #plot(all_hulls[k,],add=TRUE,col=NA,border='blue')
      
      
      for(i in 1:length(temp_poly)){
        print(length(temp_poly))
        print(nrow(temp_poly@polygons[[i]]@Polygons[[1]]@coords))
        ooga<-rbind(temp_poly@polygons[[i]]@Polygons[[1]]@coords)
      }
      
      
      
      temp_poly_chull<-chull(ooga)
      
      
      
    }
    else(next)
  }
}

writeOGR(merged_hulls,dsn=paste0(wd,'/Processed_Data'),layer=paste0(deparse(substitute(temp)),"_AutomatedOuput"),driver="ESRI Shapefile",overwrite_layer=TRUE)



## Working "test" version of the circle fitting to initial fourier curve series fitting. 

calc_circbounds<-function(circle_dat,npoints=100){ #from stack overflow plot ggplot circle https://stackoverflow.com/questions/6862742/draw-a-circle-with-ggplot2
  r=circle_dat$radius
  tt<-seq(0,2*pi,length.out=npoints)
  xx<-circle_dat$cx+r*cos(tt)
  yy<-circle_dat$cy+r*sin(tt)
  return(data.frame(x=xx,y=yy))
}

holderframe<-data.frame('Stem'=NA)
#test extraction of points per chull to calc fourier trans
for(i in 1:(length(hulls)-1)){
  SingleTree_dat<-intersect(points_df,hulls[i,]) #get set of points per hull
  Fitted_Circ<-data.frame(fitCircle(SingleTree_dat$X,SingleTree_dat$Y)) #fit circle (least squares method...from package...should switch to RANSAC)
  coordinates(Fitted_Circ)<-c("cx","cy")
  
  AdjustedCoords<-data.frame(original_x=NA,original_y=NA,adjusted_x=NA,adjusted_y=NA,dist=NA,azimuth=NA)
  
  for(k in 1:length(SingleTree_dat)){#recalculate the individual tree's point's coords w/ respect to the fitted circle cent. as origin
    dist<-sqrt(((SingleTree_dat[k,]$X-Fitted_Circ$cx)^2)+((SingleTree_dat[k,]$Y-Fitted_Circ$cy)^2))
    azimuth<-atan2((SingleTree_dat[k,]$Y-Fitted_Circ$cy),(SingleTree_dat[k,]$X-Fitted_Circ$cx))
    
    newx<- (dist*sin(azimuth))
    newy<- (dist*cos(azimuth))
    
    
    AdjustedCoords[k,]<-c(SingleTree_dat[k,]$X,SingleTree_dat[k,]$Y,newx,newy,dist,azimuth)
  }
  
  
  circ<-calc_circbounds(Fitted_Circ)
  
  polarplot<-ggplot(data=AdjustedCoords, aes(azimuth,dist))+
    geom_point()
  #coord_polar("y")
  cartplot<-ggplot(data=data.frame(SingleTree_dat),aes(X,Y))+
    geom_point()+
    geom_point(data=data.frame(Fitted_Circ),aes(cx,cy),col='red')+
    geom_path(data=circ,aes(x,y),col='red')
  grid.arrange(cartplot,polarplot,ncol=2,top=paste('Tree ',hulls[i,]$id))
  
}

#plot(AdjustedCoords$adjusted_x,AdjustedCoords$adjusted_y)

#write.csv(AdjustedCoords,file='C:/Users/note2/Documents/GitHub/ULS_ForestCensus/Processed_Data/Example_Polar_StemData.csv')
#writeOGR(hulls,dsn=getwd(),layer=paste0(deparse(substitute(hulls)),"__TestCheck__AutomatedOuput"),driver="ESRI Shapefile",overwrite_layer=TRUE)


dat<-AdjustedCoords$dist
#t<-1:max(AdjustedCoords$azimuth)
#rg<-diff(range(dat))

nff<-function(x=NULL,n=NULL,up=10L,plot=TRUE,add=FALSE,main=NULL,...){
  dff<<-fft(x)
  t=seq(from=-pi,to=pi)#may need to adjust to match radians measurment?
  nt=seq(from=-pi,to=pi-1/up,by=1/up)
  ndff<<-array(data=0,dim=c(length(nt),1L))
  ndff[1]<<-dff[1]
  if(n!=0){
    ndff[2:(n+1)]<<-dff[2:(n+1)]
    ndff[length(ndff):(length(ndff)-n+1)]<<-dff[length(x):(length(x)-n+1)]
  }
  indff<<-fft(ndff/max(x),inverse=TRUE)
  idff<<-fft(dff/max(x),inverse=TRUE)
  if(plot){
    if(!add){
      plot(x=t,y=x,pch=16L,xlab="Time",ylab="Measurment",
           main=ifelse(is.null(main),paste(n," harmonics"),main))
      lines(y=Mod(idff),x=t,col=adjustcolor(1L,alpha=0.5))
    }
    lines(y=Mod(indff),x=nt,...)
  }
  
  ret=data.frame(time=nt,y=Mod(indff))
  ret$y=((ret$y-min(ret$y)) / (max(ret$y)-min(ret$y)))*(max(x)-min(x))+min(x)
  return(ret)
}

res=nff(x=dat,n=18L,up=10L,col=2L,add=TRUE)

plot(res$time,res$y)

ggplot(data=AdjustedCoords,aes(azimuth,dist))+
  geom_point()+
  geom_line(data=res,aes(time,y))

