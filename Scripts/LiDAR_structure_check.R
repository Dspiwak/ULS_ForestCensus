#Structure_Check.R
#   Checks file structure and uses parts of script from LiDAR_TreeChulls.R 7/31/2020
#Creates a file structure for the LiDAR processing to help organize items into subdirectories...
#Relies upon presence of a "Raw_Data", "Scripts" directories in the working directory
LiDAR_structure_check<-function(slice_dir){
  file_structure<<-list(slice_dir,"las_clipped","processed_chunk") #Temporary string used to name extensions in the parent "slices" name
  for(i in 1:length(file_structure)){
    if(file_structure[[i]]==slice_dir && !dir.exists(file.path(wd,"/",file_structure[[i]]))){
      cat(paste0("Declared folder does not exist currently,","\n Do you want to create the folder and choose a LAS/LAZ file now?"))
      choice<-readline(prompt="y/n:  ")
      if(choice=="y"){
        dir.create(file.path(wd,"/",slice_dir))
        cat(paste0("Folder::","//",slice_dir," created \n","\nPlease now select a LAS/LAZ file to process"))
        file.copy(file.choose(),paste0(wd,"/",slice_dir))
      }
      else if(choice=="n"){
        print("Ok, please rerun function and select a file which exists and contains a LAS/LAZ")
        stop(cat(paste0("File and file structure does not exist.","\n Terminating script...")))
      }else(stop(cat(paste0("Please enter a valid response (y/n)"," \n Terminating Script..."))))
    }
    else if(file_structure[[i]]!=slice_dir&& !dir.exists(file.path(wd,"/",slice_dir,"/",file_structure[[i]]))){
      print(cat(paste0("Creating '",file_structure[[i]],"' folder which is needed to proceed with lidR processing")))
      dir.create(file.path(wd,"/",slice_dir,"/",file_structure[[i]]))
    }
    else if(file_structure[[i]]!=slice_dir && dir.exists(file.path(wd,"/",slice_dir,"/",file_structure[[i]]))){
      files_present<-list.files(path=paste0(wd,"/",slice_dir,"/",file_structure[[i]]),full.names = TRUE)
      if(length(files_present>0)){
        for(k in 1:length(files_present)){
          print(paste0("deleting files present:: ",files_present[[k]]))
          file.remove(files_present[[k]])
        }
      }else(print(paste0("no old files present in:: ",file_structure[[i]]," ...good!")))
    }
  }
  print("File structure check OK for lidR processing, proceeding...")
}
