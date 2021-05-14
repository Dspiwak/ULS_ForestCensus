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

rodrigues_rot<-function(P,n0,n1){
  #Rotate a set of points between two normal vectors using Rodrigues' formula
  #P= dataset of points (x,y,z coords)
  #n0= origin vector , n1= destination vector
  #Returns dataset but rotated
  
  P=as.matrix(P)
  if(size(P)[2]==1){ #reformat matrix if it is the '1D'...if is coordinates of a center
    P=t(P)
  }
  
  #get rotation k and angle theta
  n0=n0/norm(n0,type="2")
  n1=n1/norm(n1,type="2")
  k_rot=xprod(n0,n1)
  P_rot<-zeros(nrow(P),m = 3)
  
  if(norm(k_rot,type="2")!=0){
    k_rot<-k_rot/norm(k_rot,type="2")
    theta<-acos(dot(n0,n1))
    for(i in 1:nrow(P)){
      P_rot[i,]=(P[i,]*cos(theta))+(as.numeric(xprod(k_rot,P[i,]))*sin(theta))+(k_rot*dot(k_rot,P[i,])*(1-cos(theta)))
    }
  }else(P_rot=P)
  
  return(P_rot)
}

fit_RANSAC<-function(stem,thresh=.5,prob=.95,w=.5,plot=FALSE){
  #stem: A data.frame which contains a list called 'PointsInStem' which has the 'X' , 'Y' ,& 'Z' coordinates for each point in the stem
  #thresh: distance from cylinder hull which is considered an inlier
  #prob: the probability at least one random selection is error free of set n points
  #w: the probability selected data is within error tolerance
  
  temp_cent=stem%>%fit_circle(method='LM') #get a roughly fitted center by fitting a circle using 'LM' method
  
  pts=stem$PointsInStem[[1]][,1:3]%>% 
    transmute(x_data=X,y_data=Y,z_data=Z)%>% #reshape data into necessary headings
    mutate(x_data=(x_data-temp_cent[[1]])*100,y_data=(y_data-temp_cent[[2]])*100)#normalize the point coordinates about the circle fitted center
 
  Ek<-w^-3 #expected value of k (which should be exceeded by 2-3 SD of k...thus probabilistic k is determined below)
  
  k= (log(1-prob)) / (log(1-(1-w)^3)) #Number of iterations, constrained by probability that the generated model is accurate '3' is number of samples
  k<-k*Ek #k ~ = probablistic k * expected value of k 
  
  n_points=nrow(pts) #returns dimension of the pts dataset being loaded (proba dataframe which means this fucntion should be changes)
  best_inliers=NA
  outputdat=data.frame("cx"=NA,"cy"=NA,"radius"=NA,"cz"=NA)
  
  #Plot initial point data being fitted
  if(plot){
    plot(pts$x_data,pts$y_data,asp=1)
  }
  
  #Iteratively solve for Circle Params
  for(it in 1:k){
    
    sample_indices=sample(nrow(pts),size=3) #randomly select a number from the length of the point dataset to extract those x/y coords
    sample_dat=pts[sample_indices,]
    ptsmpl=sample_dat
    
    #Determine plane which describes the 3 points sampled
    #A=pt2-pt1
    #B=pt3-pt1
    
    vecA=as.numeric(ptsmpl[2,]-ptsmpl[1,]) #Generates Vector A between points 2 & 1
    vecA_norm=vecA/norm(vecA,type="2")
    vecB=as.numeric(ptsmpl[3,]-ptsmpl[1,])
    vecB_norm=vecB/norm(vecB,type="2")
    
    #Cross product of vecA and vecB results in vecC which is normal to plane
    vecC=xprod(vecA_norm,vecB_norm) #Uses a 'cross product' function which is used for physics based modelling of 3D vectors
    vecC=vecC/norm(vecC,type="2")
    
    kplane=-sum(vecC*ptsmpl[2,])
    plane_eq=c(vecC[1],vecC[2],vecC[3],kplane)
    
    #Calculate rotaion of points with rodrigues rotation eq.
    n1_rot=as.numeric(c(0,0,1))
    P_rot=rodrigues_rot(ptsmpl,vecC,n1_rot)
    
    #Find center from 3 points & intersecting lines with these points and 'center'
    ma=0
    mb=0
    while(ma==0){
      ma=(P_rot[2,2]-P_rot[1,2])/(P_rot[2,1]-P_rot[1,1]) # (y2-y1) / (x2-x1)
      mb=(P_rot[3,2]-P_rot[2,2])/(P_rot[3,1]-P_rot[2,1]) # (y3-y2) / (x3-x2)
      if(ma==0){
        P_rot=c(tail(P_rot,-1),head(P_rot,1)) #equivalent to np.roll -1 (This rearranges point order if an two points are considered a vertical line)
      }else(break)
    }
    
    #Calculate center by verifying intersection of orthogonal lines
    p_center_x=(ma*mb*(P_rot[1,2]-P_rot[3,2])+mb*(P_rot[1,1]+P_rot[2,1])-ma*(P_rot[2,1]+P_rot[3,1]))/(2*(mb-ma))
    p_center_y=-1/(ma)*(p_center_x-(P_rot[1,1]+P_rot[2,1])/2)+(P_rot[1,2]+P_rot[2,2])/2
    p_center=c('x_data'=p_center_x,'y_data'=p_center_y,'z_data'=0)
    radius=norm(p_center-P_rot[1,],type="2") #Radius is the distance from the center to the first rotation point
    
    #Recalc Rodrigues rotation
    center=rodrigues_rot(p_center,n1_rot,vecC)
    
    
    #Distance from a point to the circle's plane
    dist_pt_plane=(plane_eq[1]*pts[,1]+plane_eq[2]*pts[,2]+plane_eq[3]*pts[,3]+plane_eq[4])/sqrt(plane_eq[1]^2+plane_eq[2]^2+plane_eq[3]^2)
    #vecC_stakado=matrix(vecC, ncol=3, nrow=n_points, byrow=TRUE)
    vecC_stakado=vecC
    
    #distance from a point to the circle hull if as its perpendicular to plane
    dist_pt_inf_circle=list()
    
    tempdat=center-pts
    tempdf=data.frame(xprod(vecC,tempdat))
    func<-function(x){
      norm(x,'2')-radius
    }
    dist_pt_inf_circle=apply(tempdf,MARGIN=1,FUN=func)
    
    #distance from each point to the circle
    dist_pt=matrix(sqrt((dist_pt_inf_circle^2)+(dist_pt_plane^2)))
    
    
    
    #select indices where distance is greater than threshold
    #Distance from a point to a line
    pt_id_inliers=c(which(dist_pt<=thresh)) #Search each column and if any row has a value under the threshold, it is an inlier and get index
    
    
    if(is.null(nrow(best_inliers))){ #If initial run, then set params to default to initial calc
      best_inliers=pts[pt_id_inliers,]
      outputdat=mutate(outputdat,cx=center[1],cy=center[2],radius=radius,
                       cz=center[3],
                       Axis=list(vecC),
                       inliers=list(best_inliers))
      if(plot){
        draw.circle(outputdat$cx,outputdat$cy,radius=outputdat$radius,col=NA,border='red')
      }
    }
    else if(length(pt_id_inliers) > nrow(best_inliers) & !is.null(best_inliers)){ #If subsequent runs, 'recalc' / store new values if there are more inliers than previous random samples
      best_inliers=pts[pt_id_inliers,]
      outputdat[1,1:2]=center[1:2]
      outputdat$radius=radius
      outputdat[1,4]=center[3]
      outputdat$axis=list(vecC)
      outputdat$inliers=list(best_inliers)
      
      if(plot){
        draw.circle(outputdat$cx,outputdat$cy,radius=outputdat$radius,col=NA,border='red')
      }
    }
  }
  
  #Plot 'Final' solution
  if(plot){
    draw.circle(outputdat$cx,outputdat$cy,radius=outputdat$radius,col=NA,border='green')
  }
  #Output final solution
  return(outputdat)
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

multislice_process<-function(stems,method,...,numslice=7){
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
        diam_list[[k]]=method(stem_df,...,)[[3]]
      }else(diam_list[[k]]=NA)
      
    }
    
    ave_diam[[i]]=mean(diam_list,na.rm=TRUE) #append the average of the subslices to the list
  }
  return(ave_diam)
}

tempdat<-c(1,2,3)
ooga<-function(dat,...){
  print(sum(dat))  
  }

testfunc<-function(dat,method,add_args){
  print("loaded data")
  method(dat,add_args)
  print("ran func")
}



