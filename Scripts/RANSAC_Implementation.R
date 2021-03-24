#R Implementation of RANSAC Circle Fitting
#Based / transcribed from SeongHyunBae's python implementation: https://github.com/SeongHyunBae/RANSAC-circle-python

#Define a circle
circle_func<-function(a,b,r,x){
  return(c(sqrt(r^2-(x-a)^2)+b, -sqrt(r^2-(x-a)^2)+b))
}

#gen test circle data
init_circ<-function(){
  x_data=c()
  y_data=c()
  
  a=2
  b=3
  r=100
  
  #for(i in (a-r, a+r, 1)){ #what does this do??
  for(i in (a-r):(a+r)){
    x=i
    y=0
    y1=circle_func(a,b,r,x)[1] #could compress these down to a dataframe of "circles" containing each respective info
    y2=circle_func(a,b,r,x)[2]
    
    if(rnorm(1)>0){
      y=y1+rnorm(1)*5
    }
    else{
      y=y2+rnorm(1)*5
    }
    x_data=append(x_data,x)
    y_data=append(y_data,y)
  }
  return(data.frame('x_data'=x_data,'y_data'=y_data))
}

#create a function to make R equivalent of 'not in'
'%!in%'<-Negate('%in%')


#select 3 points randomly from the data
random_sample<-function(x){
  sample_dat=data.frame()
  sample_indices=sample(nrow(x),size=3) #randomly select a number from the length of the point dataset to extract those x/y corods
  for(i in 1: length(sample_indices)){
    sample_dat=rbind(sample_dat,c(x[sample_indices[i],])) #append the selected data to a list of paired x/y coords
  }
  return(sample_dat)
}

#Calculate A,B,C values form the sampled points

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

execute_ransac<-function(x){
  d_min=99999 #initialize a default values which will be iteratively improved
  best_model=NA
  
  for(i in 1:k){ #n should be calculated optimally for a .95 percentile...this is how many times to try to calc the optimal circle
    sample_dat=random_sample(x)
    model_dat=make_model(sample_dat)
    d_temp=eval_model(sample_dat,model_dat)
    
    if(d_min>d_temp){ #if current model is better than previous model, then update the model
      best_model=model_dat
      d_min=d_temp
      draw.circle(best_model$c_x,best_model$c_y,best_model$r)
    }
    
  }
  return(best_model)
}


testdat<-init_circ() #generate test data
plot(testdat,asp=1) #plot test data while preserving aspect ratio

k<-50 #iterate and sample 50 times (should update to the probabilistic approach for determining k)

optimal_circ<-execute_ransac(testdat)

