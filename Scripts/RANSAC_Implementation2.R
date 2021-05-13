#R Implementation of RANSAC Circle Fitting
#Based / transcribed from leomariga's python implementation: https://github.com/leomariga/pyRANSAC-3D/blob/master/pyransac3d/circle.py#L68
#Great resource for this material! : http://paulbourke.net/geometry/circlesphere/

#probablistic number of iterations to run (double check this to ensure it is based off of the number of points in the dataset?)
n=3 #number of sampled points per iteration
w=.5 #probability selected data is within error tolerance

Ek<-w^-n #expected value of k (which should be exceeded by 2-3 SD of k...thus probabalistic k is determined below)
prob=.95 #probability at least one random selection is error free of set n points

k<-(log(1-prob)) / (log(1-(1-w)^3)) #probablistic approach for determining k
k<-k*Ek #k ~ = probablistic k * expected value of k 


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


Circle_RANSAC<-function(pts,thresh=.2,iterations,plot=FALSE){
  #thresh = distance from cylinder hull which is considered an inlier
  
  n_points=nrow(pts) #returns dimension of the pts dataset being loaded (proba dataframe which means this fucntion should be changes)
  best_inliers=NA
  outputdat=data.frame("centX"=NA,"centY"=NA,"centZ"=NA,"Radius"=NA)
  
  #Plot initial point data being fitted
  if(plot){
    plot(pts$x_data,pts$y_data,asp=1)
  }
  
  #Iteratively solve for Circle Params
  for(it in 1:iterations){
    
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
      outputdat=mutate(outputdat,centX=center[1],centY=center[2],centZ=center[3],
                       Radius=radius,
                       Axis=list(vecC),
                       Inliers=list(best_inliers))
      if(plot){
        draw.circle(outputdat$centX,outputdat$centY,radius=outputdat$Radius,col=NA,border='red')
      }
    }
    else if(length(pt_id_inliers) > nrow(best_inliers) & !is.null(best_inliers)){ #If subsequent runs, 'recalc' / store new values if there are more inliers than previous random samples
      best_inliers=pts[pt_id_inliers,]
      outputdat[1,1:3]=center
      outputdat$Radius=radius
      outputdat$Axis=list(vecC)
      outputdat$Inliers=list(best_inliers)
      
      if(plot){
        draw.circle(outputdat$centX,outputdat$centY,radius=outputdat$Radius,col=NA,border='red')
      }
    }
  }
  
  #Plot 'Final' solution
  if(plot){
    draw.circle(outputdat$centX,outputdat$centY,radius=outputdat$Radius,col=NA,border='green')
  }
  #Output final solution
  return(outputdat)
}



ooga<-Circle_RANSAC(pts = testdat, thresh = .1, iterations = 100,plot=TRUE)

