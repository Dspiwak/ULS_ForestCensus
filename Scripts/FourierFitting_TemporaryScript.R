#Temporary Test Script for Trying Fourier Curve Fitting Methodologies
#This script serves as a testing area for various fourier curve fitting methods, and may contain errors.
#These errors are most likely to arise in file loading / directory setup
#Will be added to git ignore and a formal commit to the Stem_reconstructions script will be added once a methodology is decided.

#Load required packages----
required_libraries<-list("rgeos","dplyr","sp","rgdal","lidR","raster","grDevices","ggplot2","gridExtra","smacof")

if(menu(c("Yes","No"),title="Do you wish to install the required packages?")==1){
  install.packages(setdiff(required_libraries,rownames(installed.packages()))) #will install required packages
}else{print("The script may not run properly without the proper packages")}

for(package in required_libraries){
  library(package,character.only = TRUE)
}

#Set working directory---- (Currently is pointing to where this script is located...and assumes the test data is located here as well)
wd<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

#Define Required Functions----
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

calc_circbounds<-function(circle_dat,npoints=100){ #Allows plotting of fitted circle and calculates its datapoints (S.O)
  r=circle_dat$radius
  tt<-seq(0,2*pi,length.out=npoints)
  xx<-circle_dat$cx+r*cos(tt)
  yy<-circle_dat$cy+r*sin(tt)
  return(data.frame(x=xx,y=yy))
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
    print(paste0("Points in dataset:",k))
    polar_points_dfls$dfs[i]<-list(SingleTree_AdjustedCoords)
    polar_points_dfls$hullid[i]<-hulls[i,]$id

  }
  return(polar_points_dfls)
}

fit_fourier<-function(AdjustedCoords_dists_df,n,up=10L,plot=FALSE){
  dat<-AdjustedCoords_dists_df$dist
  
  dff=fft(dat)#discrete fourier transform of the distances
  t=seq(from=-pi,to=pi,length.out=length(dat))#domain range ("time") which is -pi:pi (should be only extents of the dataset?)
  ndff=array(data=0,dim=c(length(t),1L))
  ndff[1]=dff[1]
  if(n!=0){
    ndff[2:(n+1)]=dff[2:(n+1)]
    ndff[length(ndff):(length(ndff)-n+1)]=dff[length(dat):(length(dat)-n+1)]
  }
  indff=fft(ndff/length(dat),inverse=TRUE)#calculates the inverse FFT
  idff=fft(dff/length(dat),inverse=TRUE)
  
  ret=data.frame(time=t,y=Mod(indff))
  ret$y=((ret$y-min(ret$y)) / (max(ret$y)-min(ret$y)))*(max(dat)-min(dat))+min(dat) #normalize data to same period as data (-pi -pi)

  
  if(plot){ #allows plotting of data in polar and cartesian coordinate systems
    cartgraph=ggplot(data=AdjustedCoords_dists_df,aes(azimuth,dist))+
      geom_point()+
      geom_line(data=ret,aes(time,y),color='red')
    polargraph=ggplot(data=AdjustedCoords_dists_df,aes(azimuth,dist))+
      geom_point()+
      geom_line(data=ret,aes(time,y),col='red')+
      coord_polar(theta="x",direction=-1,start=pi/2)+
      scale_y_continuous(limits=c(0,max(testpolar[[1]][[1]]$dist)+.2))+
      labs(title=paste0("Harmonics: ",n))
    grid.arrange(cartgraph,polargraph,ncol=2)
  }
  return(ret)
}

#Load test datasets----
OptimalTree<-read.csv(paste0(getwd(),"/TestTree_OptimalTree"))
SubOptimalTree<-read.csv(paste0(getwd(),"/TestTree_SubOptimalTree"))
OccludedTree<-read.csv(paste0(getwd(),"/TestTree_UnevenOcclusion"))

#Convert a test dataset to a 'spatial' object and calculate initial convex hull----

dataset<-OptimalTree #Specific selected test dataset

ConvHull_points<-chull(dataset)#Calculate the convex hull (required for when looping over actual entire dataset)
dataset_Hull<-dataset[c(ConvHull_points,ConvHull_points[1]),]
dataset_Hull<-SpatialPolygons(list(Polygons(list(Polygon(dataset_Hull)),ID=1)))#SpatialPolygon for this dataset
dataset_Hull$id<-1 #set a dummy id

coordinates(dataset)<-c("x","y")#Convert the original dataset points into spatial objects
dummy_ids<-data.frame(seq(1:length(dataset)))
spdf<-SpatialPointsDataFrame(coords=dataset,data=dummy_ids)
spdf$X<-spdf$x
spdf$Y<-spdf$y


#Fit a circle to the dataset to "center" the polar coordinates about----
CircleData<-fit_circle(spdf,hulls = dataset_Hull)#If the dataset selected is the "challenge tree" then the circle fitting will not be adiquate and should use the commented code below
#CircleData<-data.frame("cx"=14.0,"cy"=10.06554,"radius"=7.205898) #original cx=5.434675 / cy=13.06554 (ONLY FOR CHALLANGE TREE)
#coordinates(CircleData)<-c("cx","cy")

#Convert to polar coordinates and fit fourier curve----
testpolar<-conv2_polar(spdf,CircleData,dataset_Hull,plot=TRUE)
res2<-fit_fourier(testpolar[[1]][[1]],n=3,up=10,plot=TRUE)
