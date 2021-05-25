basic_stats<-function(GTruth,Stems,Matched_stems,Measurement_ls){
  #Will return basic statistics regarding accuracy of each reconstruction type
  #Returns rmse,mape,etc
  #Returns a printout of the num stems detected vs total present
  #Measurement_ls is a list of strings of the reconstruction techniques to assess in the Matched_stems df
  
  stat_dat=data.frame('mae'=matrix(NA,ncol=1,nrow=length(Measurement_ls)))%>%
    mutate(mse=NA,rmse=NA,mape=NA,nmse=NA,rstd=NA)
  
  for(i in 1:length(Measurement_ls)){
    stat_dat[i,]=c(data.frame(rmserr(Matched_stems$Gtruth_dbh_cm, Matched_stems[,as.character(Measurement_ls[[i]])])))
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
      round( cor(Matched_stems$Gtruth_dbh_cm,Matched_stems[,as.character(Measurement_ls[[i]])],method="pearson"), 4))
  }
  
  cat("\n")
  
  for(i in 1:25){cat("-")}
  cat("\n Searching in tree sizes >",Smallest.Tree.Size,
             "\n Total trees present:",nrow(cliptruth),
             "\n Total trees detected:",nrow(stems),
             "\n Total trees matched:",nrow(Matched_stems))
  
  return(stat_dat)
}

plot_regression=function(x,y){
  
  lab_str=sub(".*Matched_Stems","",deparse(substitute(y)))
  lab_str=substr(lab_str ,2, (nchar(lab_str)) )
  
  p=ggplot(data=as.data.frame(Matched_Stems),aes_string(y,x,ymin=0,ymax=max(y),xmin=0,xmax=max(x)))+
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

calc_bias<-function(statdf,binwidth,plot=FALSE){
  #A function which will calculate the binned biases across a dataframe containing stem measurements
  #binwidth in cm 
  
  #Create a df to store the bias values per bin per measurement type
  biasdf<-data.frame(DBH=statdf$DBH,
                     ConvH_Bias=statdf$DBH-statdf$ConvH,
                     Pratt_Bias=statdf$DBH-statdf$Pratt,
                     LM_Bias=statdf$DBH-statdf$LM,
                     RANSAC_Bias=statdf$DBH-statdf$RANSAC)%>%
    mutate(bin=cut_width(DBH,width=binwidth,boundary=0))#Create a field with teh respective bin index value for each row
  
  #loop through the columns of the dataframe and return a bias boxplot plot 
  if(plot){
    for(i in 1:(ncol(biasdf)-2)){
      biasplot<-biasdf%>%
        ggplot(aes(x=bin,y=biasdf[,i+1]))+
        geom_boxplot()+
        ggtitle(names(biasdf)[i+1])+
        xlab("DBH (cm)")+
        ylab("Bias (cm)")
      print(biasplot)
    }
  }
  return(biasdf)
}

calc_MAE<-function(statdf,plot=FALSE){
  #A function which will calculate the binned MAEs across a dataframe containing stem measruments
  #Assumes the dataframe already contains a bin column
  
  #Create a df to store the MAE values per bin per measurment type
  MAE_df<-data.frame(ConvH_MAE=zeros(length(unique(statdf$bin)),1),
                     Pratt_MAE=zeros(length(unique(statdf$bin)),1),
                     LM_MAE=zeros(length(unique(statdf$bin)),1),
                     RANSAC_MAE=zeros(length(unique(statdf$bin)),1))
  
  MAEbins=unique(statdf$bin) #A list which stores the bin values to be looped through
  
  for(i in 2:(ncol(statdf)-1)){ #loop through the statdf, but exclude the first column (DBH) & last column (bin)
    MAE<-list() #A list of the MAE values for a specific measrument type to be populated by loop
    for(k in 1:length(MAEbins)){ #loop through all bins
      tempdf<-filter(statdf,bin==MAEbins[k])#Filter the input df to only include values of a specific bin at one time
      MAE[k]<-rmserr(tempdf$DBH,tempdf[,i])$mae #calc MAE
    }
    MAE_df[[i-1]]<-MAE
  }
  MAE_df<-mutate(MAE_df,bin=unique(statdf$bin))
  
  #Plot the MAE using GGPLOT
  if(plot){
    for(i in 1:(ncol(MAE_df)-1)){
      MAE_plot<-ggplot(MAE_df,aes(bin,MAE_df[,i]))+
        geom_bar(stat="identity")+
        ggtitle(names(MAE_df)[i])+
        xlab("DBH (cm)")+
        ylab("MAE (cm)")
      print(MAE_plot)
    }
  }
  return(MAE_df)
}
