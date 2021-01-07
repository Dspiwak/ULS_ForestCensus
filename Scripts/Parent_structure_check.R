#Check Dirs of the working directory and flag if a specific folder is not present

present_dirs<-list.dirs(path=wd,full.names = FALSE)
required_dirs<-list("Raw_Data","Scripts","Processed_Data") #Minimum list of required directories to process data
required_scripts<-list("FGeo_Processing.R","LiDAR_structure_check.R","Parent_structure_check.R","LiDAR_clustering.R","LiDAR_initial_processing.R",
                       "Stem_reconstructions.R","Stem_matching.R") #list of scripts required to fully preform process

for(i in required_dirs){
  cat("\n")
  
  if(i %in% present_dirs==TRUE){ #Required folder is present
    if(i=='Scripts'){ #Checks to ensure all necessary scripts are present
      present_scripts=list.files(path=paste0(wd,'/Scripts'))
      
      for(k in required_scripts){
        if(k %in% present_scripts==FALSE){
          cat(paste0("Missing script: ",k))
          cat("Please ensure all scripts are present before continuing to process")
          }
        else if(k %in% present_scripts==TRUE){
          next
          }
        else(cat(paste0("Non-recognized script, or other file present in scripts directory: ",k)))
      }
    }
  }
  
  else if(i %in% present_dirs==FALSE){ #Required folder is missing
    if(i==""){next} #Skips the first "parent" directory listing (".")
    cat(paste0(i," is not present..Creating ",wd,"/",i))
    dir.create(i) #creates missing directory, but will not put necessary data in here (user must do)
    cat(" Ensure all neccisary files are placed in this newly created directory, or the rest of the script will fail")
  }
  
  else(cat(paste0("Non-required directory present '", i,"' but will continue running")))
}

