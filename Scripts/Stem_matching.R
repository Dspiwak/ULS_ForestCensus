match_stems<-function(polys,DetectionRange,is.multislice){
  loop_polys<-polys
  ids<-loop_polys$clusid
  Stem_Buff<-gBuffer(Clipped.GTruth,byid=TRUE,width=Clipped.GTruth$rbh_m)
  matched_trees<-data.frame()
  bufwidth<-DetectionRange #originally 2.5 but after discussion with shawn we increased to 4 to test
  i<-0
  while(1<length(loop_polys)){
    i<-i+1
    Chull_Centroid<-gCentroid(loop_polys[loop_polys$clusid==ids[i],])
    val<-ids[i]
    SearchArea<-gBuffer(Chull_Centroid,width=bufwidth)
    detect_gtruth<-(intersect(SearchArea,Stem_Buff)) #matches by stem_buffers, could use "over" and match by gtruth "cents" aswell...but..this retains data and is easier
    detect_multipoly<-intersect(SearchArea,loop_polys)
    if(!is.null(detect_gtruth) && nrow(detect_gtruth)==1){#if there is only 1 tree in the search area then calc the min dist and ensure another poly isnt closer
      offset2_gtruth<-pointDistance(Chull_Centroid,gCentroid(detect_gtruth))
      azimuth<-"temp_val"
      if(nrow(detect_multipoly)>1){#if there is another chull poly in search dist, then calc its dist to this point, then compare and assign
        for(k in 1:length(detect_multipoly)){
          detect_multipoly$dist2gtruth[k]<-pointDistance(gCentroid(detect_multipoly[k,]),gCentroid(detect_gtruth))
        }
        detect_multipoly_df<-as.data.frame(detect_multipoly)%>%
          arrange(dist2gtruth)
        closest_chull<-slice(detect_multipoly_df,1)
        if(closest_chull$dist2gtruth<offset2_gtruth){
          matched_trees<-bind_rows(matched_trees,data.frame(detect_gtruth))
          matched_trees$clusid[i]<-val
          matched_trees$chull_dbh_cm[i]<-closest_chull$diam3_euc
          matched_trees$offset[i]<-closest_chull$dist2gtruth
          matched_trees$azimuth[i]<-azimuth
          matched_trees$testcase[i]<-1
          Stem_Buff<-Stem_Buff[-which(gIntersects(Stem_Buff,detect_gtruth,byid=TRUE)),]
          loop_polys<-loop_polys[loop_polys$clusid!=val,]
        }
        else if(closest_chull>offset2_gtruth){
          matched_trees<-bind_rows(matched_trees,data.frame(detect_gtruth))
          matched_trees$clusid[i]<-val
          matched_trees$chull_dbh_cm[i]<-loop_polys[loop_polys$clusid==ids[i],]$diam3_euc
          matched_trees$offset[i]<-offset2_gtruth
          matched_trees$azimuth[i]<-azimuth
          matched_trees$testcase[i]<-1.2
          Stem_Buff<-Stem_Buff[-which(gIntersects(Stem_Buff,detect_gtruth,byid=TRUE)),]
          loop_polys<-loop_polys[loop_polys$clusid!=val,]
        }
      }else{
        matched_trees<-bind_rows(matched_trees,data.frame(detect_gtruth))
        matched_trees$clusid[i]<-val
        matched_trees$chull_dbh_cm[i]<-loop_polys[loop_polys$clusid==ids[i],]$diam3_euc
        matched_trees$offset[i]<-offset2_gtruth
        matched_trees$azimuth[i]<-azimuth
        matched_trees$testcase[i]<-1.3
        Stem_Buff<-Stem_Buff[-which(gIntersects(Stem_Buff,detect_gtruth,byid=TRUE)),]
        loop_polys<-loop_polys[loop_polys$clusid!=val,]
      }
    }#detect if there is only 1 ground truth in region
    else if(!is.null(detect_gtruth) && nrow(detect_gtruth)>1){ #multiple gtruth trees in search area...then select the closest if there isnt another chull nearby
      for(j in 1:length(detect_gtruth)){
        offset2_gtruth<-pointDistance(Chull_Centroid,gCentroid(detect_gtruth[j,]))
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
      offset2_gtruth<-as.numeric(pointDistance(Chull_Centroid,detect_gtruth))
      offset2_knn2<-as.numeric(pointDistance(Chull_Centroid,knn2))
      if(nrow(detect_multipoly)>1){ #also if there is another chull poly, then calc their mindists to this closest point
        for(k in 1:length(detect_multipoly)){
          detect_multipoly$dist2gtruth[k]<-pointDistance(gCentroid(detect_multipoly[k,]),detect_gtruth)
        }
        detect_multipoly<-as.data.frame(detect_multipoly)%>%
          arrange(dist2gtruth)
        closest_chull<-slice(detect_multipoly,1)
        if(closest_chull$dist2gtruth<offset2_gtruth){
          # matched_trees<-bind_rows(matched_trees,data.frame(detect_gtruth))
          # matched_trees$chull_dbh_cm[i]<-closest_chull$diam3_euc
          # matched_trees$offset[i]<-closest_chull$dist2gtruth
          # matched_trees$azimuth[i]<-azimuth
          # matched_trees$testcase[i]<-2
          matched_trees<-bind_rows(matched_trees,data.frame(knn2))
          matched_trees$clusid[i]<-val
          matched_trees$chull_dbh_cm[i]<-loop_polys[loop_polys$clusid==ids[i],]$diam3_euc
          matched_trees$offset[i]<-offset2_knn2
          matched_trees$azimuth[i]<-azimuth
          matched_trees$testcase[i]<-2
          Stem_Buff<-Stem_Buff[-which(gIntersects(Stem_Buff,detect_gtruth,byid=TRUE)),]
          loop_polys<-loop_polys[loop_polys$clusid!=val,]
          #might need to remove second detected poly too? the knn2 poly/ closest chull?
        }else{
          matched_trees<-bind_rows(matched_trees,data.frame(detect_gtruth))
          matched_trees$clusid[i]<-val
          matched_trees$chull_dbh_cm[i]<-loop_polys[loop_polys$clusid==ids[i],]$diam3_euc
          matched_trees$testcase[i]<-2.2
          Stem_Buff<-Stem_Buff[-which(gIntersects(Stem_Buff,detect_gtruth,byid=TRUE)),]
          loop_polys<-loop_polys[loop_polys$clusid!=val,]
        }
      }else{
        matched_trees<-bind_rows(matched_trees,data.frame(detect_gtruth))
        matched_trees$clusid[i]<-val
        matched_trees$chull_dbh_cm[i]<-loop_polys[loop_polys$clusid==ids[i],]$diam3_euc
        matched_trees$testcase[i]<-2.3
        Stem_Buff<-Stem_Buff[-which(gIntersects(Stem_Buff,detect_gtruth,byid=TRUE)),]
        loop_polys<-loop_polys[loop_polys$clusid!=val,]
      }
    }else{ #otherwise, if there is no ground truth in the search area, plot the results, and continue on for later review
      matched_trees[i,]<-NA
      matched_trees$clusid[i]<-val
      matched_trees$chull_dbh_cm[i]<-loop_polys[loop_polys$clusid==ids[i],]$diam3_euc
      matched_trees$testcase[i]<-3
      loop_polys<-loop_polys[loop_polys$clusid!=val,]
    }
  }
  
  matched_trees<-filter(matched_trees,!is.na(dbh_cm))
  coordinates(matched_trees)<-matched_trees[,c("NADX","NADY")]
  crs(matched_trees)<-NAD83_2011
  
  return(matched_trees)
}