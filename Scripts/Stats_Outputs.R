'%!in%'=Negate('%in%')#create a temporary 'not in' function

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
      "\n Total stems present:",nrow(GTruth),
      "\n Total stem IDs not detected:",nrow(GTruth[GTruth$stemID%!in%Matched_stems$Gtruth_stemID,]),
      "\n Total stems detected (With False Positives):",nrow(Stems),
      "\n Total stems matched:",nrow(Matched_stems),
      "\n Total stems matched & detected (absolute):",nrow(GTruth[GTruth$tag%in%Matched_stems$Gtruth_tag,]),
      "\n Number of excessively noisey stems >1.5x the max dbh in the gtruth: ",sum(rowSums(Matched_stems[,7:10]>max(GTruth$dbh_cm)*1.5)),"\n",
      
      "\n ** Note, the number of stems may not add up due to the inclusion of duplicate stemIDs...",
      "\n ** However, it may be undesirable to remove these due to detection of shared treeID values where the stem tag may differ",
      "\n ** Furthermore, assess possibility for 'double matching' to occur ")
  
  cat("\n")
  
  for(i in 1:25){cat("-")}
  cat("\n Average MAPE across all bins")
  for(i in 2:(ncol(statdf)-1)){
    cat("\n")
    cat(paste0(names(statdf)[i],": ",round((rmserr(statdf$DBH,statdf[,i])$mape*100),2),"%"))
    }
  
  return(stat_dat)
}

plot_detected<-function(GTruth,Matched_Stems,width=10,plot_totals=FALSE,plot_percentages=FALSE){
   presencedf=data.frame(GTruth)%>% #temp dataframe which assesses presence of a stem in both the matched stems and ground truth
    mutate(detected=case_when(tag%in%Matched_Stems$Gtruth_tag~1, #1 if in both datasets (ie present / detected)
                              tag%!in%Matched_Stems$Gtruth_tag~0))#0 if not in both and 'not detected'
  if(plot_percentages){
    p<-ggplot(presencedf,aes(dbh_cm,fill=factor(detected)))+
      geom_histogram(binwidth=width)+
      stat_bin(
        aes(label=ifelse(..count..>width*.1,
                         paste(round((..count../tapply(..count..,..x..,sum)[as.character(..x..)]),3)*100,'%')
                         ,''),group=factor(detected),fill=factor(detected)),
        geom="label_repel",direction='y',force=0,alpha=.8,binwidth=width, show.legend=FALSE)+ #Percentage per bin per class factor
      #stat_bin(data=presencedf%>%filter(detected==1),aes(y=min(..count..)-5,label=..count..),
               #geom="label",position='identity',binwidth=width,fill='white')+ #Total number detected & matched stems per column
      theme(text=element_text(size=16))+
      labs(title="Stem Detection Rates",x="DBH",y="Number of Stems Accuratly Matched")+
           #caption=paste('Total Accuratley Matched: ',nrow(Matched_Stems),'/',nrow(GTruth)))+
      guides(fill=guide_legend(title="Detection Status"))+
      scale_fill_discrete(labels=c("Missing","Present & Matched"))+
      scale_x_continuous(breaks=seq(0,max(presencedf$dbh_cm)+width,width))
  }

  if(plot_totals){
    p<-ggplot(presencedf,aes(dbh_cm))+
      geom_histogram(aes(fill=factor(detected)),binwidth=15)+
      stat_bin(aes(y=min(..count..),label=..count..),
               geom="label",binwidth=15,position='stack',colour='black')+
      labs(title="Number of Detected Stems",x="DBH",y="Number of Stems Accuratly Matched")+
      guides(fill=guide_legend(title="Detection Status"))+
      scale_fill_discrete(labels=c("Missing","Present & Matched"))
  }
  
  return(print(p))
}

plot_regression<-function(df,x,y,GTruth=NULL){
  
  lab_str=sub(".*Matched_Stems","",deparse(substitute(y)))#append the args string to 'point' to that column in the df
  lab_str=substr(lab_str ,2, (nchar(lab_str)-1) ) #rextract the name as a string so that it can be used to label title (may not be needed?)
  
  p=ggplot(data=df,aes_string(x,y))+
    geom_pointdensity()+ #density plot of scatter points
    scale_color_viridis_c()+ #set color scale
    geom_segment(aes(x=0,y=0,xend=100,yend=100),linetype='dashed',size=1.2)+ #plot 1:1 reference line
    geom_smooth(method="gam",formula=y~x)+ #plot linear regression line
    stat_cor(method='pearson',label.x.npc="left",label.y.npc="top")+ #add pearson's correlation coef. and p-val
    stat_poly_eq(aes(label=paste(..rr.label..,..eq.label..,sep="~~~")),
                 label.x.npc="left",label.y.npc=.9,formula=y~x,parse=TRUE)+ #add r2 and linreg line eq vals #https://stackoverflow.com/questions/37494969
    theme(text=element_text(size=16))+
    labs(title=paste0("Reconstruction Method: ",lab_str))+
    xlab("True DBH (cm)")+
    ylab("Predicted DBH (cm)")
  
  #adjust y limits due to outliers (zoom plot)
  if(!is.null(GTruth)){
    p=p+ylim(c(NA,max(GTruth$dbh_cm)*1.5))
  }
  else if(is.null(GTruth)){
    p=p+coord_cartesian(ylim=c(0,quantile(df$y,.9)))
  }
  
  
  return(p)
}

calc_bias<-function(statframe,binwidth,RemoveBinsUnder=10,plot=FALSE,plot_stacked=FALSE){
  #A function which will calculate the binned biases across a dataframe containing stem measurements
  #binwidth in cm 
  
  #remove bins with few returns
  selectedbins<-statframe%>%
    dplyr::group_by(bin)%>%
    summarise(n=n())%>%
    filter(n>RemoveBinsUnder)
  statdf<-statframe%>%
    filter(bin%in%selectedbins$bin)
  
  
  
  
  
  
  
  #Create a df to store the bias values per bin per measurement type
  biasdf<-data.frame(DBH=statdf$DBH,
                     ConvH_Bias=statdf$ConvH-statdf$DBH,
                     Pratt_Bias=statdf$Pratt-statdf$DBH,
                     LM_Bias=statdf$LM-statdf$DBH,
                     RANSAC_Bias=statdf$RANSAC-statdf$DBH)%>%
    mutate(bin=cut_width(DBH,width=binwidth,boundary=0))#Create a field with the respective bin index value for each row
  
  #Create a df to store the respective median values for each bin step for each column
  summstats<-biasdf%>%
    dplyr::group_by(bin)%>%
    summarize(DBH=mean(DBH),ConvH_Bias=median(ConvH_Bias),Pratt_Bias=median(Pratt_Bias),LM_Bias=median(LM_Bias),RANSAC_Bias=median(RANSAC_Bias))
  
  #loop through the columns of the dataframe and return a bias boxplot plot 
  if(plot|plot_stacked){
    
    #Creates a named list of ggplot objects based on the available columns in the biasdf
    plotlist=sapply(names(biasdf)[-grep("bin",names(biasdf))], function(col){ #loops through available columns in above biasdf
      ggplot(biasdf,aes_string(x="bin",y=col))+ #plot the current columns data
        geom_boxplot(outlier.shape = NA)+
        geom_label(data=summstats,aes(x=bin,label=round(..y..,1)),nudge_y = 30,size=3)+ #add a median value label to each bin...note '..y..' accesses the values in the current 'col' column
        coord_cartesian(ylim=c(-50,50))+
        labs(title=col,
             x="DBH (cm)",
             y="Bias (cm)")+
        theme(text=element_text(size=16),axis.text.x=element_text(angle=90))
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
  MAE_df=mutate(MAE_df,bin=as.character(unique(statdf$bin)))%>%
    mutate(bin=gsub('\\(','',bin))%>%
    mutate(bin=gsub('\\[|\\]','',bin))%>%
    mutate(startbin=strsplit(bin,','), endbin=strsplit(bin,','))
  for(i in 1:nrow(MAE_df)){
    MAE_df$startbin[[i]]=as.numeric(MAE_df$startbin[[i]][1])
    MAE_df$endbin[[i]]=as.numeric(MAE_df$endbin[[i]][2])
  }
  MAE_df=MAE_df%>%
    mutate_if(is.list,as.numeric)
  for(i in 1:nrow(MAE_df)){
    MAE_df$bin[[i]]=mean(c(MAE_df$startbin[[i]],MAE_df$endbin[[i]]))
  }
  MAE_df=mutate(MAE_df,stepval=as.numeric(bin))%>%
    arrange(startbin)%>%
    dplyr::select(-bin,-startbin,-endbin)
  
  #Plot the MAE using GGPLOT
  if(plot){
    Molten=reshape2::melt(MAE_df,id.vars="stepval")
    MAE_plot=ggplot(Molten,aes(x=stepval,y=value,colour=variable))+
      geom_smooth(method="loess",se=FALSE,formula=y~x)+
      scale_x_continuous(breaks=seq(0,max(Molten$stepval),by=15),limits=c(0,max(Molten$stepval)))+
      scale_y_continuous(breaks=seq(0,max(Molten$stepval),by=15))+
      coord_equal(ratio=1)+
      theme(text=element_text(size=16))+
      guides(color=guide_legend(title="Method"))+
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
      MAPE[k]<-rmserr(tempdf$DBH,tempdf[,i])$mape*100 #calc MAE
    }
    MAPE_df[[i-1]]=MAPE
  }
   MAPE_df=mutate(MAPE_df,bin=as.character(unique(statdf$bin)))%>%
    mutate(bin=gsub('\\(','',bin))%>%
    mutate(bin=gsub('\\[|\\]','',bin))%>%
    mutate(startbin=strsplit(bin,','), endbin=strsplit(bin,','))
  for(i in 1:nrow(MAPE_df)){
    MAPE_df$startbin[[i]]=as.numeric(MAPE_df$startbin[[i]][1])
    MAPE_df$endbin[[i]]=as.numeric(MAPE_df$endbin[[i]][2])
  }
  MAPE_df=MAPE_df%>%
    mutate_if(is.list,as.numeric)
  for(i in 1:nrow(MAPE_df)){
    MAPE_df$bin[[i]]=mean(c(MAPE_df$startbin[[i]],MAPE_df$endbin[[i]]))
  }
  MAPE_df=mutate(MAPE_df,stepval=as.numeric(bin))%>%
    arrange(startbin)%>%
    dplyr::select(-bin,-startbin,-endbin)
  
  
  #Plot the MAE using GGPLOT
  if(plot){
    Molten=reshape2::melt(MAPE_df,id.vars="stepval")
    MAPE_plot=ggplot(Molten,aes(x=stepval,y=value,colour=variable))+
      geom_smooth(method="loess",se=FALSE,formula=y~x)+
      coord_cartesian(ylim=c(0,100))+
      theme(text=element_text(size=16))+
      guides(color=guide_legend(title="Method"))+
      xlab("DBH (cm)")+
      ylab("MAPE (%)")
    print(MAPE_plot)
  }
  
  return(MAPE_df)
}
