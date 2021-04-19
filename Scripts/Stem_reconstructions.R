conv_hulls<-function(clusters_df,plot=FALSE){
  #create a bunch of empty lists which will be needed later
  polylist<-list()
  PPConvex<-list()
  dbh1<-list()
  maxdist<-list()
  clusids<-list()
  
  #loop over individual clusters (singular trees in theory) and generate convex hulls
  for(i in 1:length(unique(clusters_df$id))){
    tempclusid<-i
    SingleTree_Points<-data.frame(clusters_df) %>%
      filter(clusters_df$id==tempclusid)
    SingleTree_Matrix<-cbind(SingleTree_Points$X,SingleTree_Points$Y)
    temp<-chull(SingleTree_Matrix) #could use gConvexHull rather than chull (would make code more readable later)
    Chull_Points<-SingleTree_Matrix[c(temp,temp[1]),]
    temp_df<-data.frame(Chull_Points)
    coordinates(temp_df)<-c("X1","X2")
    crs(temp_df)<-SetCRS
    polylist[[i]]<-SpatialPolygons(list(Polygons(list(Polygon(Chull_Points)),ID=tempclusid)))
    cent1<-gCentroid(polylist[[i]])
    crs(cent1)<-SetCRS
    dist2cent<-list()
    dist2cent<-spDistsN1(Chull_Points,cent1,longlat=FALSE)#a method of dist calculation ##MEASURMENT WILL BE OFF BECAUSE REPEAT LAST POLY POINT IN CHULLPOINTS //BIAS IN MEAN CALC
    PPConvex[[i]]<-nrow(SingleTree_Points)
    dbh1[[i]]<-(mean(dist2cent)*2)*100
    #calc max distances in polys -- to remove fallen trees which have a point distance greater than max dbh in plot + sum
    tempmax<-0
    for(j in 1:nrow(Chull_Points)){
      tempmax[[j]]<-(sqrt(((Chull_Points[j,1]-Chull_Points[1,1])^2)+((Chull_Points[j,2]-Chull_Points[1,2])^2)))
    }
    maxdist[[i]]<-(max(tempmax)*2)*100
    clusids[[i]]<-tempclusid
    
    
    #manually plot the graphs with measurments in it
    if(plot){
      title<-paste0("DBH of tree ID: ",i,"\n func: ", dbh1[[i]], " cm ")
      plot(SingleTree_Points$X,SingleTree_Points$Y,main=title,asp=1)
      points(Chull_Points,col="green")
      lines(Chull_Points, col="red")
      plot(cent1,add=TRUE,col="blue",asp=1)
      startpnt<-as.numeric(cent1@coords)
      for(p in 2:nrow(Chull_Points)){
        measurement<-paste0(dbh1[[i]]/2," cm")#measurment from convexhull point to center ~cm
        endpnt<-as.numeric(Chull_Points[p,])
        segments(startpnt[1],startpnt[2],endpnt[1],endpnt[2],lty="dashed",col="darkgrey")#draws lines between the borderpoint and center of poly
        #text(endpnt[1],endpnt[2],labels=measurement)#adds text
      }
    }
    
  }
  #data conversion from dataframes and lists to spatialpolygondataframes for export to GIS prgm
  PPConvex<-do.call(rbind,PPConvex)
  dbh1<-do.call(rbind,dbh1)
  maxdist<-do.call(rbind,maxdist)
  ClusterID<-do.call(rbind,clusids)
  dbhframe<-data.frame("ClusterID"=ClusterID,"PointsInHull"=PPConvex,"diam"=dbh1,"maxdist"=maxdist)%>%
    dplyr::select(ClusterID,PointsInHull,diam,maxdist)
  
  #Create spatial polygons and remove noise (fallen debris, branches, & shrubs)
  polys<-do.call(bind,polylist)
  polys<-SpatialPolygonsDataFrame(Sr=polys,data=dbhframe,FALSE)
  polys<-polys[polys$diam<=(max(Ground.Truth$dbh_cm)+20),] #removes polygons that have ANY SECTIONAL DIAMETER (DBH) larger than the maximum ground truth DBH size +20cm
  polys<-polys[polys$maxdist<=(max(Ground.Truth$dbh_cm)+20),] #removes polygons that have ANY MEASURMENT (convex-hull side lengths) larger than the maximum ground truth DBH size +20cm
  crs(polys)<-NAD83_2011
  
  return(polys)
}

fit_convhull<-function(stem,quantile,plot=FALSE){
  #Fits a convex hull to a set of data and will calculate average distance from gravity center to convex hull points
  
  
  stem_info=stem$PointsInStem[[1]][,1:2]
  coordinates(stem_info)=c("X","Y")
  
  chull=gConvexHull(stem_info)
  chull_coords=data.frame(chull@polygons[[1]]@Polygons[[1]]@coords)
  coordinates(chull_coords)=c("x","y")
  
  cent<-gCentroid(stem_info)
  rDists<-spDistsN1(chull_coords,cent,longlat=FALSE)
  
  quantile_dists=spDistsN1(stem_info,cent,longlat=FALSE)
  quantile_dists=quantile_dists*quantile
  
  dat=c( cent$x, cent$y , mean(rDists)*2*100, stem$StemID) #write to DF (*100 * 2 converts to cm & diam)
  
  if(plot){
    print('function not yet implemented')
  }
  return(dat)
}

calc_circbounds<-function(circle_dat,npoints=100){ #from stack overflow plot ggplot circle https://stackoverflow.com/questions/6862742/draw-a-circle-with-ggplot2
  if(class(circle_dat)=="data.frame"){
    r=circle_dat$radius
    tt=seq(0,2*pi,length.out=npoints)
    xx=circle_dat$cx+r*cos(tt)
    yy=circle_dat$cy+r*sin(tt)
  }
  else if(class(circle_dat)=="list" | is.vector(circle_dat)){
    r=circle_dat[[3]]
    tt=seq(0,2*pi,length.out=npoints)
    xx=circle_dat[[1]]+r*cos(tt)
    yy=circle_dat[[2]]+r*sin(tt)
  }
  return(data.frame(x=xx,y=yy))
} #creates "points" along curve of ellipse...primarily for plotting

fit_circle<-function(stem,method,plot=FALSE){ #fits a circle to each of the hulls in a dataset and returns a spatial points dataframe based on the centroids and contains attribute for radius
  
  dat<-data.frame("cx"=NA,"cy"=NA,"radius"=NA,"StemID"=NA) #create DF for function
  
  if(method=="Pratt"){#Pratt Algebraic fit
    circ_info=CircleFitByPratt(stem$PointsInStem[[1]][,1:2])
    dat=c( circ_info[1], circ_info[2] , circ_info[3]*2*100, stem$StemID) #write to DF (*100 * 2 converts to cm & diam)
  }
  
  else if (method=="LM"){ #Levenberg-Marquardt Reduced Method
    tempcoords=stem$PointsInStem[[1]]
    coordinates(tempcoords)=c("X","Y")
    
    cent=gCentroid(tempcoords) #get centroid of points for initial LM-Reduced (a,b) guess
    
    circ_info=LMreducedCircleFit(stem$PointsInStem[[1]][,1:2], c(data.frame(cent@coords)$x,data.frame(cent@coords)$y) ,1,1,20)
    dat=c( circ_info[1], circ_info[2], circ_info[3]*2*100, stem$StemID)
  }
  
  if(plot){
    tempdat=dat
    tempdat[[3]]=tempdat[[3]]/2/100
    
    circ=calc_circbounds(tempdat)#get a dataframe of points along circle curve
    
    p=ggplot(data=stem$PointsInStem[[1]],aes(X,Y))+
      geom_point()+
      coord_fixed()+ #retain aspect ratio
      labs(title=paste0("DBH of tree ID: ",dat[4],"\n", dat[3], " cm ","\n Method:",method))+
      geom_path(data=circ,aes(x,y),col='red') #plot circle
    print(p)
  }
  #coordinates(dat)<-c("cx","cy")
  return(dat)
}

RANSAC_circle<-function(points,hulls,prob,w,numsample=3,plot=FALSE){
  '%!in%'= Negate('%in%')  #temporary function / value which will be used to check if values are not in a list
  
  random_sample<-function(x,numsample){#select 3 points randomly from the data
    sample_dat=data.frame()
    sample_indices=sample(nrow(x),size=numsample) #randomly select a number from the length of the point dataset to extract those x/y corods (numsample typically set to 3)
    for(i in 1: length(sample_indices)){
      sample_dat=rbind(sample_dat,c(x[sample_indices[i],])) #append the selected data to a list of paired x/y coords
    }
    return(sample_dat)
  }
  
  make_model<-function(sample_dat){
    pt1=sample_dat[1,]
    pt2=sample_dat[2,]
    pt3=sample_dat[3,]
    
    #perform matrix algebra to calculate A,B,C matrices to determine centroid of circle
    A=cbind( c( pt2[[1]]-pt1[[1]] , pt2[[2]]-pt1[[2]]) , c(pt3[[1]]-pt2[[1]] , pt3[[2]]-pt2[[2]]) )
    B=rbind( c(pt2[[1]]^2-pt1[[1]]^2 + pt2[[2]]^2-pt1[[2]]^2) , c(pt3[[1]]^2-pt2[[1]]^2+pt3[[2]]^2-pt2[[2]]^2) )
    inv_A=solve(A) #calcualtes inverse of matrix 'A'
    
    c_x=(inv_A %*% B)/2 #dot product
    c_y=(inv_A %*% B)/2
    c_x=c_x[1,] #centroid x-coord
    c_y=c_y[2,] #centroid y-coord
    
    r=sqrt((c_x-pt1[[1]])^2 + (c_y-pt1[[2]])^2) #radius
    
    return(data.frame('c_x'=c_x , 'c_y'=c_y , 'r'=r))
  }
  
  eval_model<-function(sample_dat,model_dat){
    d=0
    c_x=model_dat$c_x
    c_y=model_dat$c_y
    r=model_dat$r
    
    for(i in 1:nrow(sample_dat)){ #may need to adjust to fix loop length
      dis=sqrt( (sample_dat$x_data[i]-c_x)^2 + (sample_dat$y_data[i]-c_y)^2 ) #problem here? warning message
      
      if(dis >= r){
        d= d+(dis-r)
      }
      else{ d=d+(r-dis) }
    }
    
    return(d)
  }
  
  fitted_RANSAC=data.frame('c_x'=NA,'c_y'=NA,'r'=NA,'clusid'=NA)
  for( i in 1:length(hulls)){
    SingleTree_points = intersect(points,hulls[i,]) #isolate a singular tree's set of points
    SingleTree_points<-data.frame('x_data'=SingleTree_points$X,'y_data'=SingleTree_points$Y)
    
    if(plot){
      plot(SingleTree_points$x_data,SingleTree_points$y_data,asp=1)
    }
    
    d_min=99999 #initialize a default values which will be iteratively improved upon for each tree
    best_model=NA
    
    Ek<-w^-numsample #expected value of k (which should be exceeded by 2-3 SD of k...thus probabilistic k is determined below)
    
    k= (log(1-prob)) / (log(1-(1-w)^numsample)) #Number of iterations, constrained by probability that the generated model is accurate
    k<-k*Ek #k ~ = probablistic k * expected value of k 
    
    
    for(j in 1:k){ #n should be calculated optimally for a .95 percentile...this is how many times to try to calc the optimal circle
      sample_dat=random_sample(SingleTree_points,numsample)
      model_dat=make_model(sample_dat)
      d_temp=eval_model(sample_dat,model_dat)
      
      if(d_min>d_temp){ #if current model is better than previous model, then update the model
        best_model=model_dat
        d_min=d_temp
        if(plot){
          draw.circle(best_model$c_x,best_model$c_y,best_model$r)
        }
      }
    }
    fitted_RANSAC[i,]<-c(best_model$c_x,best_model$c_y,best_model$r,hulls[i,]$id)
  }
  return(fitted_RANSAC)
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
  
  m=mean(dat)
  sigma=sd(dat)
  vec=vector()
  for(i in 1:length(dat)){
    if(abs(dat[i]-m)<sigma){
      vec=c(vec,dat[i])
    }
  }
  dat=vec
  
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
  
  if(plot){ #allows plotting of data in polar and cartesian coordinate systems
    cartgraph=ggplot(data=AdjustedCoords_dists_df,aes(azimuth,dist))+
      geom_point()+
      geom_line(data=ret,aes(time,y),color='red')
    polargraph=ggplot(data=AdjustedCoords_dists_df,aes(azimuth,dist))+
      geom_point()+
      #geom_point(data=data.frame(dat),aes(time,dist),color='red')+
      geom_line(data=ret,aes(time,y),col='red')+
      coord_polar(theta="x",direction=-1,start=pi/2)+
      #scale_y_continuous(limits=c(0,max(AdjustedCoords_data[[1]][[1]]$dist)+.2))+
      labs(title=paste0("Harmonics: ",n))
    grid.arrange(cartgraph,polargraph,ncol=2)
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

multislice_process<-function(stems,method,additional_args=FALSE,plot_args=FALSE,numslice=7){
  #A function which will allow the user to define the number of slices with which the main slice should be divided into.
  #This method will average the multiple slice diams together for each respective stem
  ave_diam=NA
  
  sliceHbins=seq( round(min(stems$PointsInStem[[1]]$Z)) , round(max(stems$PointsInStem[[1]]$Z)) , length.out=numslice) #sequence of sub slice heights
  
  
  for( i in 1:nrow(stems)){
    stem_df<-data.frame('PointsInStem'=matrix(NA,ncol=1,nrow=1)) #temporary dataframe which holds subslice datapoints
    
    diam_list=NA#temporary list of diameters to calculate average for the singular stem across multiple sub slices
    
    for(k in 1:(numslice-1)){
      #set a upper and lower bounds to filter out data to fulfill each subslice req. Z
      lowerBounds=sliceHbins[[k]]
      upperBounds=sliceHbins[[k+1]]
      
      
      stem_df$PointsInStem=list(filter(stems$PointsInStem[[i]],Z>=lowerBounds & Z<=upperBounds)) #get subslice of points for each tree
      
      if(!is.null(stem_df$PointsInStem[[1]]$X) & nrow(stem_df$PointsInStem[[1]])>=4){ #Ensures that whatever method is being run will have at least 4 points per subslice to define a model
        diam_list[[k]]=method(stem_df,additional_args,plot_args)[3]
      }else(diam_list[[k]]=NA)
      
    }
    
    ave_diam[[i]]=mean(diam_list,na.rm=TRUE) #append the average of the subslices to the list
  }
  return(ave_diam)
}


ooga<-function(dat,somearg=NA){
  print(sum(dat))
  if(!is.na(somearg)){
    print("somearg works")
  }
}

testfunc<-function(dat,method,add_args){
  print("loaded data")
  method(dat,add_args)
  print("ran func")
}



