las_process<-function(rawLASdirectory_name,lowersliceH,uppersliceH){
 
  #lidR Processing Engine Structure (allows for parallel processing later and automated "chunking" of the pointcloud)
  proc_las=function(las,lowersliceH,uppersliceH){
    las<-readLAS(las)
    if(is.empty(las))return(NULL)
    las<-classify_ground(las,csf()) %>%
      normalize_height(knnidw())%>%
      filter_poi(Z>lowersliceH)%>%
      filter_poi(Z<uppersliceH)
    las<-filter_poi(las,buffer==0)
    if(is.empty(las))return(NULL)
    return(las)
  }
  proc_las.LAScatalog=function(las,lowersliceH,uppersliceH)
  {
    opt_select(las)<-"*"
    options<-list(need_output_file=TRUE,need_buffer=TRUE,drop_null=TRUE)
    output<-catalog_sapply(las,proc_las,lowersliceH=lowersliceH,uppersliceH=uppersliceH, .options = options)
    return(output)
  }
  
  rawlas<-catalog(paste0(getwd(),"/",rawLASdirectory_name),select="icrna",pattern=".laz")
  opt_output_files(rawlas)<-paste0(getwd(),"/Processed_Data","/",file_structure[[2]],"/las_clipped_{ID}")
  crs(rawlas)<-NAD83_2011
  clip_roi(rawlas,clipregion)
  myproj<-catalog(paste0(getwd(),"/Processed_Data/",file_structure[[2]]))
  opt_output_files(myproj)<-paste0(getwd(),"/Processed_Data/",file_structure[[3]],"/processed_chunk_{ID}")
  opt_chunk_size(myproj)<-(myproj$Max.X-myproj$Min.X)/4
  opt_progress(myproj)<-TRUE
  crs(myproj)<-NAD83_2011
  myproj<<-myproj
  
  #preform point cloud normalization and slicing in parallel 
  plan(multisession,workers=availableCores()/2)
  set_lidr_threads(2)
  output<-proc_las.LAScatalog(myproj,lowersliceH,uppersliceH) 
  
  #create a list of the processed chunks and write las files to a directory
  las.list<-list.files(paste0(getwd(),"/Processed_Data/",file_structure[[3]],"/"),pattern="*.las",full.names=TRUE)
  import.list<-llply(las.list,readLAS)
  for(i in 1:length(import.list)){
    import.list[[i]]<-as.data.frame(import.list[[i]]@data)
  }
  lasdf<<-bind_rows(import.list)%>%
    dplyr::select(X,Y,Z)
  
  return(lasdf)
}