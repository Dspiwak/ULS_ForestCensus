conv_hulls<-function(clusters_df){
  #create a bunch of empty lists which will be needed later
  polylist<-list()
  PPConvex<-list()
  dbh1<-list()
  dbh2<-list()
  dbh3<-list()
  maxdist<-list()
  
  #loop over individual clusters (singular trees in theory) and generate convex hulls
  for(i in 1:length(unique(clusters_df$id))){
    SingleTree_Points<-data.frame(clusters_df) %>%
      filter(clusters_df$id==i)
    SingleTree_Matrix<-cbind(SingleTree_Points$X,SingleTree_Points$Y)
    temp<-chull(SingleTree_Matrix)
    Chull_Points<-SingleTree_Matrix[c(temp,temp[1]),]
    temp_df<-data.frame(Chull_Points)
    coordinates(temp_df)<-c("X1","X2")
    crs(temp_df)<-SetCRS
    polylist[[i]]<-SpatialPolygons(list(Polygons(list(Polygon(Chull_Points)),ID=i)))
    cent1<-gCentroid(polylist[[i]])
    crs(cent1)<-SetCRS
    dist2cent<-list()
    dist2cent<-spDistsN1(Chull_Points,cent1,longlat=FALSE)#a method of dist calculation ##MEASURMENT WILL BE OFF BECAUSE REPEAT LAST POLY POINT IN CHULLPOINTS //BIAS IN MEAN CALC
    dist3cent<-pointDistance(cent1,Chull_Points,lonlat=FALSE)#another method of distance calculation
    dist4cent<-0
    for(k in 1:(nrow(Chull_Points)-1)){
      dist4cent[[k]]<-(sqrt(((Chull_Points[k,1]-cent1@coords[1,1])^2)+((Chull_Points[k,2]-cent1@coords[1,2])^2))) #calc euclidean distance manually
    }
    PPConvex[[i]]<-nrow(SingleTree_Points)
    dbh1[[i]]<-(mean(dist2cent)*2)*100
    dbh2[[i]]<-(mean(dist3cent)*2)*100
    dbh3[[i]]<-(mean(dist4cent)*2)*100
    #calc max distances in polys -- to remove fallen trees which have a point distance greater than max dbh in plot + sum
    tempmax<-0
    for(j in 1:nrow(Chull_Points)){
      tempmax[[j]]<-(sqrt(((Chull_Points[j,1]-Chull_Points[1,1])^2)+((Chull_Points[j,2]-Chull_Points[1,2])^2)))
    }
    maxdist[[i]]<-(max(tempmax)*2)*100
    
    
    #manually plot the graphs with measurments in it
    #title<-paste0("DBH of tree ID: ",i,"\n func: ", dbh1[[i]], " cm or euc_dist: ", dbh3[[i]]," cm")
    #plot(SingleTree_Points$X,SingleTree_Points$Y,main=title)
    #points(Chull_Points,col="green")
    #lines(Chull_Points, col="red")
    #plot(cent1,add=TRUE,col="blue")
    #startpnt<-as.numeric(cent1@coords)
    #for(p in 2:nrow(Chull_Points)){
    #measurement<-paste0(dbh3[[i]]/2," cm")#measurment from convexhull point to center ~cm
    #endpnt<-as.numeric(Chull_Points[p,])
    #segments(startpnt[1],startpnt[2],endpnt[1],endpnt[2],lty="dashed",col="darkgrey")#draws lines between the borderpoint and center of poly
    #text(endpnt[1],endpnt[2],labels=measurement)#adds text
    #}
    
    
  }
  #data conversion from dataframes and lists to spatialpolygondataframes for export to GIS prgm
  PPConvex<-do.call(rbind,PPConvex)
  dbh1<-do.call(rbind,dbh1)
  dbh2<-do.call(rbind,dbh2)
  dbh3<-do.call(rbind,dbh3)
  maxdist<-do.call(rbind,maxdist)
  dbhframe<-data.frame("PointsPerHull"=PPConvex,"diam"=dbh1,"diam2"=dbh2,"diam3_euc"=dbh3,"maxdist"=maxdist)%>%
    dplyr::mutate(clusid=row_number())%>%
    dplyr::select(PointsPerHull,diam,diam2,diam3_euc,maxdist,clusid)
  
  #Create spatial polygons and remove noise (fallen debris, branches, & shrubs)
  polys<-do.call(bind,polylist)
  polys<-SpatialPolygonsDataFrame(Sr=polys,data=dbhframe,FALSE)
  polys<-polys[polys$diam3_euc<=(max(Ground.Truth$dbh_cm)+20),] #removes polygons that have ANY SECTIONAL DIAMETER (DBH) larger than the maximum ground truth DBH size +20cm
  polys<-polys[polys$maxdist<=(max(Ground.Truth$dbh_cm)+20),] #removes polygons that have ANY MEASURMENT (convex-hull side lengths) larger than the maximum ground truth DBH size +20cm
  crs(polys)<-NAD83_2011
  
  return(polys)
}