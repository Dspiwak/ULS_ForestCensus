basic_stats<-function(GTruth,Stems,Matched_stems,Measurement_ls){
  #Will return basic statistics regarding accuracy of each reconstruction type
  #Returns rmse,mape,etc
  #Returns a printout of the num stems detected vs total present
  #Measurement_ls is a list of strings of the reconstruction techniques to assess in the Matched_stems df
  
  stat_dat=data.frame('mae'=matrix(NA,ncol=1,nrow=length(Measurement_ls)))%>%
    mutate(mse=NA,rmse=NA,mape=NA,nmse=NA,rstd=NA)
  
  for(i in 1:length(Measurement_ls)){
    stat_dat[i,]=c(data.frame(rmserr(Matched_stems$dbh_cm, Matched_stems[,as.character(Measurement_ls[[i]])]@data[[1]])))
  }
  rownames(stat_dat)<-Measurement_ls
  stat_dat=round(stat_dat,2)
  
  #Stat Readout
  for(i in 1:50){cat("-")}
  
  cat("\n")
  print(stat_dat)
  cat("\n")
  
  for(i in 1:25){cat("-")}
  
  for(i in 1:length(Measurement_ls)){
    cat("\n",Measurement_ls[[i]]," Pearsons correlation: ",
      round( cor(Matched_stems$dbh_cm,Matched_stems[,as.character(Measurement_ls[[i]])]@data[[1]],method="pearson"), 4))
  }
  
  cat("\n")
  
  for(i in 1:25){cat("-")}
  cat("\n Searching in tree sizes >",Smallest.Tree.Size,
             "\n Total trees present:",nrow(cliptruth),
             "\n Total trees detected:",nrow(stems))
  
  return(stat_dat)
}

plot_regression=function(x,y){
  
  lab_str=sub(".*matched_stems","",deparse(substitute(y)))
  lab_str=substr(lab_str ,2, (nchar(lab_str)) )
  
  p=ggplot(data=as.data.frame(matched_stems),aes_string(y,x,ymin=0,ymax=max(y),xmin=0,xmax=max(x)))+
    geom_smooth(method="lm",formula=y~x)+
    geom_point()+
    coord_equal()+
    stat_cor(label.y=max(y))+
    stat_regline_equation(label.y=max(y)-5)+
    labs(title=paste0("Reconstruction Method: ",lab_str))+
    xlab("True DBH (cm)")+
    ylab("Predicted DBH (cm)")
  print(p)
}