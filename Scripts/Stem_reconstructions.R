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
    temp<-chull(SingleTree_Matrix) #could use gConvexHull rather than chull (would make code more readable later)
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

calc_circbounds<-function(circle_dat,npoints=100){ #from stack overflow plot ggplot circle https://stackoverflow.com/questions/6862742/draw-a-circle-with-ggplot2
  r=circle_dat$radius
  tt<-seq(0,2*pi,length.out=npoints)
  xx<-circle_dat$cx+r*cos(tt)
  yy<-circle_dat$cy+r*sin(tt)
  return(data.frame(x=xx,y=yy))
} #creates "points" along curve of ellipse...primarily for plotting

fit_circle<-function(points,hulls){ #fits a circle to each of the hulls in a dataset and returns a spatial points dataframe based on the centroids and contains attribute for radius
  dat<-data.frame("cx"=NA,"cy"=NA,"radius"=NA,"clusid"=NA)
  for(i in 1:length(hulls)){
    SingleTree_points<-intersect(points,hulls[i,]) #extracts points from an individual hull when overlapping
    temp<-data.frame(fitCircle(SingleTree_points$X,SingleTree_points$Y))#fit circle (least squares method...from package...should switch to RANSAC)
    dat[i,]<-c(temp$cx,temp$cy,temp$radius,hulls$id[[i]])
    }
  coordinates(dat)<-c("cx","cy")#makes spatial points data to act as origin for polar coords...radius field can be used as a DBH method technically (although somewhat poor accuracy currently)
  return(dat)
}

conv2_polar<-function(points_df,Fitted_CircleData,hulls,plot=FALSE){
  polar_points_dfls<-data.frame("dfs"=matrix(NA,ncol=1,nrow=length(hulls)))
  for(i in 1:length(hulls)){
    minibuff_hull<-gBuffer(hulls[i,],width=0.1) #must make a small buffer around the chull because points along the boundary will not be "detected" with the intersection
    SingleTree_points<-intersect(points_df,minibuff_hull)
    SingleTree_CircleDat<-Fitted_CircleData[i,]
    SingleTree_AdjustedCoords<-data.frame(original_x=NA,original_y=NA,dist=NA,azimuth=NA)
    
    dist<-NA
    azimuth<-NA
    
    for(k in 1:length(SingleTree_points)){#recalculate the individual tree's point's coords w/ respect to the fitted circle cent. as origin
      dist<-sqrt(((SingleTree_points[k,]$X-SingleTree_CircleDat$cx)^2)+((SingleTree_points[k,]$Y-SingleTree_CircleDat$cy)^2))
      azimuth<-atan2((SingleTree_points[k,]$Y-SingleTree_CircleDat$cy),(SingleTree_points[k,]$X-SingleTree_CircleDat$cx))
      
      #newx<- (dist*sin(azimuth))
      #newy<- (dist*cos(azimuth))
      
      SingleTree_AdjustedCoords[k,]<-c(SingleTree_points[k,]$X,SingleTree_points[k,]$Y,dist,azimuth)
    }
     print(k)
     polar_points_dfls$dfs[i]<-list(SingleTree_AdjustedCoords)
     polar_points_dfls$hullid[i]<-hulls[i,]$id
     
     
     
    if(plot){
      circ<-calc_circbounds(SingleTree_CircleDat)
      polarplot<-ggplot(data=SingleTree_AdjustedCoords, aes(azimuth,dist))+
        geom_point()
      cartplot<-ggplot(data=data.frame(SingleTree_points),aes(X,Y))+
        geom_point()+
        geom_point(data=data.frame(SingleTree_CircleDat),aes(cx,cy),col='red')+
        geom_path(data=circ,aes(x,y),col='red')
      grid.arrange(cartplot,polarplot,ncol=2,top=paste('Tree ',hulls[i,]$id))
    }
  }
  return(polar_points_dfls)
}
 
fit_fourier<-function(AdjustedCoords_dists_df,n,up=10L,plot=FALSE){
  dat<-AdjustedCoords_dists_df$dist
  
  dff=fft(dat)#discrete fourier transform of the distances
  t=seq(from=-pi,to=pi,length.out=length(dat))#domain range ("time") [same length as the data]
  nt=seq(from=-pi,to=pi,length.out=length(dat))#upsampled domain range ... this may be affecting the overall domain /shift?
  ndff=array(data=0,dim=c(length(nt),1L))
  ndff[1]=dff[1]
  if(n!=0){
    ndff[2:(n+1)]=dff[2:(n+1)]
    ndff[length(ndff):(length(ndff)-n+1)]=dff[length(dat):(length(dat)-n+1)]
  }
  indff=fft(ndff/length(dat),inverse=TRUE)#calculates the 
  idff=fft(dff/length(dat),inverse=TRUE)

  ret=data.frame(time=nt,y=Mod(indff))
  ret$y=((ret$y-min(ret$y)) / (max(ret$y)-min(ret$y)))*(max(dat)-min(dat))+min(dat) #normalize data to same period as data (-pi -pi)
  blek<<-Mod(idff)
  
  print(paste0("t= ",length(t)))
  print(paste0("nt= ",length(nt)))
  print(paste0("ndff length= ",length(ndff)))
  print(paste0("dff length= ",length(dff)))
  print(paste0("indff length= ",length(indff)))
  print(paste0("idff length= ",length(idff)))
  
  
  if(plot){
    #plot(ret$time,ret$y)
    # if(!add){
    #   plot(x=t,y=dat,pch=16L,xlab="Time",ylab="Measurment",
    #        main=ifelse(is.null(main),paste(n," harmonics"),main))
    #   lines(y=Mod(idff),x=t,col=adjustcolor(1L,alpha=0.5))
    # }
    graph<-ggplot(data=AdjustedCoords_dists_df,aes(azimuth,dist))+
      geom_point()+
      geom_line(data=ret,aes(time,y),color='red')
    print(graph)
    
  }
  
  
  L=  2#arc length in polar coordinates
  ds=2*pi #effective domain range (set to 2pi for now because ULS "assumes" equal coverage about face of tree)
  dbh=(2*L)/ds
  
  
  #using fourier curve and fitted circle to "fill in" fractions of missing areas of fourier curve (where not sufficient data coverage) [currently not implemented]
  #dbh2=(L/pi)+(1-(ds/2*pi))*pc (where pc is diam of fitted circle)
  return(ret)
}

remove_fourieroutliers<-function(fouriercurve_df,polarcord_data){
  newpolardat=polarcord_data
  newpolardat$dist2curve=NA
  
  for(i in 1:nrow(newpolardat)){
    newpolardat$dist2curve[i]=sqrt((fouriercurve_df$y[i]^2+newpolardat$dist[i]^2)+
                                     (2*fouriercurve_df$y[i]*newpolardat$dist[i]*cos(fouriercurve_df$time[i]-newpolardat$azimuth[i])))
  }
  
  newpolardat=filter(newpolardat,dist2curve<3)
  return(newpolardat)
}

calc_arclength<-function(fouriercurve_df,polarcord_data){
  arclength=list()
  for(i in 1:length(fouriercurve_df$y)-1){#calculate arc length given a dataframe with datapoints along a fourier curve
    arclength[i]=sqrt((fouriercurve_df$y[i]^2+fouriercurve_df$y[i+1]^2)+
                        (2*fouriercurve_df$y[i]*fouriercurve_df$y[i+1]*cos(fouriercurve_df$time[i]-fouriercurve_df$time[i+1])))
  }
  arclength=do.call(sum,arclength)
  
  effectivedomain=sum(abs(range(polarcord_data$azimuth))) #calculates the effective domain range as the sum of min & max theta values (essentially 2pi right now)
  
  fourier_dbh=(2*arclength)/effectivedomain
  
  return(fourier_dbh)
}