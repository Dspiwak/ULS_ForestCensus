match_stems<-function(stems,GTruth,DetectionRange){
  #A function which will match clustered LiDAR data to GroundTruth points based on distance to centroid 
  #A "faux" 2-knn approach which relies on a search area "DetectionRange" that will be focused on to attempt to match stems
  
  
  ##Really cludgy way of converting nested dataframe back to 1 dataframe of convex hulls
  tempstems=stems
  ids=data.frame("StemID"=matrix(NA,ncol=1,nrow=nrow(tempstems)))#ID Variable (Used for filtering cleanly)
  
  chulllist=list()
  
  for(i in 1:nrow(tempstems)){
    coordinates(tempstems$PointsInStem[[i]])=c("X","Y")
    chulllist[[i]]=gConvexHull(tempstems$PointsInStem[[i]])
    ids$StemID[[i]]=tempstems$StemID[[i]]
  }
  chulllist=do.call(bind,chulllist)
  loop_polys=SpatialPolygonsDataFrame(Sr=chulllist,data=ids,FALSE)
  stem_dat=stems%>%
    dplyr::select(-PointsInStem)%>%
    mutate(StemID=as.integer(StemID))
  polys=loop_polys
  ids=loop_polys$StemID
  ##
  
  
  Stem_Buff=gBuffer(GTruth,byid=TRUE,width=GTruth$rbh_m)
  matched_trees=data.frame()
  
  i=0
  while(1<length(loop_polys)){
    print(i)
    i<-i+1
    Chull_Centroid<-gCentroid(loop_polys[loop_polys$StemID==ids[i],])
    val<-ids[i]
    SearchArea<-gBuffer(Chull_Centroid,width=DetectionRange)
    detect_gtruth<-(intersect(SearchArea,Stem_Buff)) #matches by stem_buffers, could use "over" and match by gtruth "cents" aswell...but..this retains data and is easier
    detect_multipoly<-intersect(SearchArea,loop_polys)
    if(!is.null(detect_gtruth) && nrow(detect_gtruth)==1){#if there is only 1 tree in the search area then calc the min dist and ensure another poly isnt closer
      offset2_gtruth<-pointDistance(Chull_Centroid,gCentroid(detect_gtruth),lonlat=FALSE)
      azimuth<-"temp_val"
      if(nrow(detect_multipoly)>1){#if there is another chull poly in search dist, then calc its dist to this point, then compare and assign
        for(k in 1:length(detect_multipoly)){
          detect_multipoly$dist2gtruth[k]<-pointDistance(gCentroid(detect_multipoly[k,]),gCentroid(detect_gtruth),lonlat=FALSE)
        }
        detect_multipoly_df<-as.data.frame(detect_multipoly)%>%
          arrange(dist2gtruth)
        closest_chull<-slice(detect_multipoly_df,1)
        if(closest_chull$dist2gtruth<offset2_gtruth){
          matched_trees<-bind_rows(matched_trees,data.frame(detect_gtruth))
          matched_trees$StemID[i]<-val
          #matched_trees$chull_dbh_cm[i]<-closest_chull$diam
          matched_trees$offset[i]<-closest_chull$dist2gtruth
          matched_trees$azimuth[i]<-azimuth
          matched_trees$testcase[i]<-1
          Stem_Buff<-Stem_Buff[-which(gIntersects(Stem_Buff,detect_gtruth,byid=TRUE)),]
          loop_polys<-loop_polys[loop_polys$StemID!=val,]
        }
        else if(closest_chull>offset2_gtruth){
          matched_trees<-bind_rows(matched_trees,data.frame(detect_gtruth))
          matched_trees$StemID[i]<-val
          #matched_trees=bind_rows(stem_dat[stem_dat$StemID==ids[i],])
          matched_trees$offset[i]<-offset2_gtruth
          matched_trees$azimuth[i]<-azimuth
          matched_trees$testcase[i]<-1.2
          Stem_Buff<-Stem_Buff[-which(gIntersects(Stem_Buff,detect_gtruth,byid=TRUE)),]
          loop_polys<-loop_polys[loop_polys$StemID!=val,]
        }
      }else{
        matched_trees<-bind_rows(matched_trees,data.frame(detect_gtruth))
        matched_trees$StemID[i]<-val
        #matched_trees$chull_dbh_cm[i]<-loop_polys[loop_polys$StemID==ids[i],]$diam
        matched_trees$offset[i]<-offset2_gtruth
        matched_trees$azimuth[i]<-azimuth
        matched_trees$testcase[i]<-1.3
        Stem_Buff<-Stem_Buff[-which(gIntersects(Stem_Buff,detect_gtruth,byid=TRUE)),]
        loop_polys<-loop_polys[loop_polys$StemID!=val,]
      }
    }#detect if there is only 1 ground truth in region
    else if(!is.null(detect_gtruth) && nrow(detect_gtruth)>1){ #multiple gtruth trees in search area...then select the closest if there isnt another chull nearby
      for(j in 1:length(detect_gtruth)){
        offset2_gtruth<-pointDistance(Chull_Centroid,gCentroid(detect_gtruth[j,]),lonlat=FALSE)
        detect_gtruth$dist2closestgtruth[j]<-offset2_gtruth
      }
      detect_gtruth_df<-as.data.frame(detect_gtruth)%>%
        arrange(dist2closestgtruth)
      knn2<-slice(detect_gtruth_df,2)
      detect_gtruth_df<-slice(detect_gtruth_df,1)
      coordinates(detect_gtruth_df)<-detect_gtruth_df[,c("NADX","NADY")]
      crs(detect_gtruth_df)<-NAD83_2011
      coordinates(knn2)<-knn2[,c("NADX","NADY")]
      crs(knn2)<-NAD83_2011
      detect_gtruth<-detect_gtruth_df
      offset2_gtruth<-as.numeric(pointDistance(Chull_Centroid,detect_gtruth,lonlat=FALSE))
      offset2_knn2<-as.numeric(pointDistance(Chull_Centroid,knn2,lonlat=FALSE))
      if(nrow(detect_multipoly)>1){ #also if there is another chull poly, then calc their mindists to this closest point
        for(k in 1:length(detect_multipoly)){
          detect_multipoly$dist2gtruth[k]<-pointDistance(gCentroid(detect_multipoly[k,]),detect_gtruth,lonlat=FALSE)
        }
        detect_multipoly<-as.data.frame(detect_multipoly)%>%
          arrange(dist2gtruth)
        closest_chull<-slice(detect_multipoly,1)
        if(closest_chull$dist2gtruth<offset2_gtruth){
          matched_trees<-bind_rows(matched_trees,data.frame(knn2))
          matched_trees$StemID[i]<-val
          #matched_trees$chull_dbh_cm[i]<-loop_polys[loop_polys$StemID==ids[i],]$diam
          matched_trees$offset[i]<-offset2_knn2
          matched_trees$azimuth[i]<-azimuth
          matched_trees$testcase[i]<-2
          Stem_Buff<-Stem_Buff[-which(gIntersects(Stem_Buff,detect_gtruth,byid=TRUE)),]
          loop_polys<-loop_polys[loop_polys$StemID!=val,]
          #might need to remove second detected poly too? the knn2 poly/ closest chull?
        }else{
          matched_trees<-bind_rows(matched_trees,data.frame(detect_gtruth))
          matched_trees$StemID[i]<-val
          #matched_trees$chull_dbh_cm[i]<-loop_polys[loop_polys$StemID==ids[i],]$diam
          matched_trees$testcase[i]<-2.2
          Stem_Buff<-Stem_Buff[-which(gIntersects(Stem_Buff,detect_gtruth,byid=TRUE)),]
          loop_polys<-loop_polys[loop_polys$StemID!=val,]
        }
      }else{
        matched_trees<-bind_rows(matched_trees,data.frame(detect_gtruth))
        matched_trees$StemID[i]<-val
        #matched_trees$chull_dbh_cm[i]<-loop_polys[loop_polys$StemID==ids[i],]$diam
        matched_trees$testcase[i]<-2.3
        Stem_Buff<-Stem_Buff[-which(gIntersects(Stem_Buff,detect_gtruth,byid=TRUE)),]
        loop_polys<-loop_polys[loop_polys$StemID!=val,]
      }
    }else{ #otherwise, if there is no ground truth in the search area, plot the results, and continue on for later review
      matched_trees[i,]<-NA
      matched_trees$StemID[i]<-val
      #matched_trees$chull_dbh_cm[i]<-loop_polys[loop_polys$StemID==ids[i],]$diam
      matched_trees$testcase[i]<-3
      loop_polys<-loop_polys[loop_polys$StemID!=val,]
    }
  }
  
  matched_trees=left_join(matched_trees,stem_dat,by="StemID")%>%
    filter(!is.na(dbh_cm))
  coordinates(matched_trees)<-matched_trees[,c("NADX","NADY")]
  crs(matched_trees)<-NAD83_2011
  
  return(matched_trees)
}

conv2_spdf<-function(stems){
  polylist<-list() #create a temporary list which will be populated by polygons
  for(i in 1:nrow(stems)){
    tempcoords<-stems$PointsInStem[[i]] #get 1 row of data from 'stems' nested dataframe
    coordinates(tempcoords)<-c("X","Y")
    crs(tempcoords)<-SetCRS
    polylist[[i]]<-gConvexHull(tempcoords) #assign the convex hull as the shape of the polygon
  }
  polys<-do.call(bind,polylist)#bind the list to a 'SpatialPolygons' which is digestible to the next line
  tempspdf<-SpatialPolygonsDataFrame(polys,data=stems%>%select(-PointsInStem)%>%mutate(StemID=as.numeric(StemID))) #create the SPDF and bring in the additional columns
  
  return(tempspdf)
}