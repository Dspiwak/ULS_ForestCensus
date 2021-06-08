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

plot_detected<-function(GTruth,Matched_Stems,plot_totals=FALSE,plot_percentages=FALSE){
  
  '%!in%'=Negate('%in%')#create a temporary 'not in' function
  
  presencedf=data.frame(GTruth)%>% #temp dataframe which assesses presence of a stem in both the matched stems and ground truth
    mutate(detected=case_when(treeID%in%Matched_Stems$Gtruth_ID~1, #1 if in both datasets (ie present / detected)
                              treeID%!in%Matched_Stems$Gtruth_ID~0))#0 if not in both and 'not detected'
  if(plot_percentages){
  p<-ggplot(presencedf,aes(dbh_cm,fill=factor(detected)))+
    geom_histogram(binwidth=15)+
    stat_bin(
      aes(label=ifelse(..count..>30,
                       paste(round((..count../tapply(..count..,..x..,sum)[as.character(..x..)]),3)*100,'%')
                       ,''),group=factor(detected),fill=factor(detected)),
      geom="label_repel",direction='y',force=0,alpha=.8,binwidth=15, show.legend=FALSE)+
    labs(title="Number of Detected Stems",x="DBH",y="Number of Stems")+
    guides(fill=guide_legend(title="Detection Status"))+
    scale_fill_discrete(labels=c("Missing","Present"))
  }
  
  if(plot_totals){
    p<-ggplot(presencedf,aes(dbh_cm))+
      geom_histogram(aes(fill=factor(detected)),binwidth=15)+
      stat_bin(aes(y=min(..count..),label=..count..),
               geom="label",binwidth=15,position='stack',colour='black')+
      labs(title="Number of Detected Stems",x="DBH",y="Number of Stems")+
      guides(fill=guide_legend(title="Detection Status"))+
      scale_fill_discrete(labels=c("Missing","Present"))
  }
  
  return(print(p))
}

plot_regression<-function(df,x,y){
  
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

calc_bias<-function(statdf,binwidth,plot=FALSE,plot_stacked=FALSE){
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
  if(plot|plot_stacked){
    
    #Creates a named list of ggplot objects based on the available columns in the biasdf
    plotlist=sapply(names(biasdf)[-grep("bin",names(biasdf))], function(col){
      ggplot(biasdf,aes_string(x="bin",y=col))+
        geom_boxplot(outlier.shape = NA)+
        coord_cartesian(ylim=c(-110,110))+
        labs(title=col,
             x="DBH (cm)",
             y="Bias (cm)")+
        theme(axis.text.x=element_text(angle=90))
    },simplify=FALSE)
    
    if(plot){ #plot columns individually
      print(plotlist) 
    }
    
    if(plot_stacked){ #plot columns stacked together in one figure
      print(wrap_plots(plotlist[2:length(plotlist)]))
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

calc_MAPE<-function(statdf,plot=FALSE){
  #A function which will calculate the binned MAEs across a dataframe containing stem measruments
  #Assumes the dataframe already contains a bin column
  
  #Create a df to store the MAE values per bin per measurment type
  MAPE_df<-data.frame(ConvH_MAPE=zeros(length(unique(statdf$bin)),1),
                     Pratt_MAPE=zeros(length(unique(statdf$bin)),1),
                     LM_MAPE=zeros(length(unique(statdf$bin)),1),
                     RANSAC_MAPE=zeros(length(unique(statdf$bin)),1))
  
  MAPEbins=unique(statdf$bin) #A list which stores the bin values to be looped through
  
  for(i in 2:(ncol(statdf)-1)){ #loop through the statdf, but exclude the first column (DBH) & last column (bin)
    MAPE<-list() #A list of the MAE values for a specific measrument type to be populated by loop
    for(k in 1:length(MAPEbins)){ #loop through all bins
      tempdf<-filter(statdf,bin==MAPEbins[k])#Filter the input df to only include values of a specific bin at one time
      MAPE[k]<-rmserr(tempdf$DBH,tempdf[,i])$mape #calc MAE
    }
    MAPE_df[[i-1]]=MAPE
  }
  MAPE_df=mutate(MAPE_df,bin=unique(statdf$bin))%>%
    mutate(stepval=ceiling(seq(from=10,to=max(as.numeric(statdf[,1])),length.out=length(unique(statdf$bin)))))%>%
    mutate_if(is.list,as.numeric)%>%
    dplyr::select(-bin)
  
  #Plot the MAE using GGPLOT
  if(plot){
    Molten=reshape2::melt(MAPE_df,id.vars="stepval")
    MAPE_plot=ggplot(Molten,aes(x=stepval,y=value,colour=variable))+
      geom_smooth(method="loess",se=FALSE,formula=y~x)+
      coord_cartesian(ylim=c(0,1.2))+
      xlab("DBH (cm)")+
      ylab("MAPE (%)")
    print(MAPE_plot)
  }
  
  return(MAPE_df)
}
