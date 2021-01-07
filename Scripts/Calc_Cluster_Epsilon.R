get.elbow.points.indices<-function(x,y,threshold){
  d1<-diff(y)/diff(x) #first derivative
  d2<-diff(d1)/diff(x[-1]) #second derivative
  indices<-which(abs(d2)>threshold)
  return(indices)
}

auto_EPS<-function(df,minpnts,threshold){
  temp<-kNNdist(df,k=minpnts,all=TRUE)
  index<-as.numeric(rownames(data.frame(temp)))
  val<-rowMeans(temp)
  temp<-data.frame("x"=index,"y"=val)%>%
    arrange(desc(y))%>%
    mutate(index=row_number(y))
  indices<-get.elbow.points.indices(temp$index,temp$y,threshold)
  crit_points<<-data.frame("x"=temp$index[indices],"y"=temp$y[indices])
  p<<-ggplot(data=temp,aes(x=index,y=y))+
    geom_point()+
    geom_point(data=crit_points,aes(x=x,y=y),color="red")+
    geom_hline(aes(yintercept=min(crit_points$y)),color="blue",linetype="dashed")+
    geom_text(aes(0,min(crit_points$y),label=min(crit_points$y),vjust=-.6),hjust=-2)
  p
  return(min(crit_points$y))
}