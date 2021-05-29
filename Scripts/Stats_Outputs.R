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

plot_regression=function(df,x,y){
  
  lab_str=sub(".*Matched_Stems","",deparse(substitute(y)))
  lab_str=substr(lab_str ,2, (nchar(lab_str)) )
  
  p=ggplot(data=df)+
    geom_smooth(aes_string(y,x),method="gam",formula=y~x)+
    geom_pointdensity(aes_string(y,x))+
    scale_color_viridis_c()+
    coord_cartesian(xlim=c(0,quantile(df$x,.9)),
                    ylim=c(0,quantile(df$y,.9)))+
    #stat_cor(label.y=max(y))+
    #stat_regline_equation(label.y=max(y)-5)+
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
                     ConvH_Bias=statdf$ConvH-statdf$DBH,
                     Pratt_Bias=statdf$Pratt-statdf$DBH,
                     LM_Bias=statdf$LM-statdf$DBH,
                     RANSAC_Bias=statdf$RANSAC-statdf$DBH)%>%
    mutate(bin=cut_width(DBH,width=binwidth,boundary=0))#Create a field with teh respective bin index value for each row
  
  #loop through the columns of the dataframe and return a bias boxplot plot 
  if(plot){
    for(i in 1:(ncol(biasdf)-2)){
      biasplot<-biasdf%>%
        ggplot(aes(x=bin,y=biasdf[,i+1]))+
        geom_boxplot(outlier.shape = NA)+
        ggtitle(names(biasdf)[i+1])+
        xlab("DBH (cm)")+
        ylab("Bias (cm)")
      
      ylim1=boxplot.stats(biasdf[,i+1],coef=3)$stats[c(1,5)]
      
      biasplot=biasplot+coord_cartesian(ylim=c(-110,110))
      
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
    MAE_df[[i-1]]=MAE
  }
  MAE_df=mutate(MAE_df,bin=unique(statdf$bin))%>%
    mutate(stepval=ceiling(seq(from=10,to=max(as.numeric(statdf[,1])),length.out=length(unique(statdf$bin)))))%>%
    mutate_if(is.list,as.numeric)%>%
    dplyr::select(-bin)
  
  #Plot the MAE using GGPLOT
  if(plot){
    Molten=reshape2::melt(MAE_df,id.vars="stepval")
    MAE_plot=ggplot(Molten,aes(x=stepval,y=value,colour=variable))+
      geom_smooth(method="loess",se=FALSE,formula=y~x)+
      coord_equal()+
      xlab("DBH (cm)")+
      ylab("MAE (cm)")
    print(MAE_plot)
  }
  
  return(MAE_df)
}
