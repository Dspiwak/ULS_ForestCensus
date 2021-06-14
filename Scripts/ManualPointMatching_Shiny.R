#to create a interactive graph which will generate a CSV with the points placed in the graph
library(dplyr)
library(raster)
library(ggplot2)
library(shiny)
library(shinyjs)
library(shinythemes)
library(rgdal)
library(rgeos)
library(sf)

####-------------UI-----------------------------------------------------------------------
ui<-navbarPage(theme=shinytheme("darkly"),
               fluid=TRUE,
  #Title for app--
  "Match The Convex Hull With The Respective Ground Truth 'Game'",
  
  #create second tab in UI for data initialization
  tabPanel("Load Data",
           
           #chull UI element
           fileInput('LiDAR_Data','Choose LiDAR Convex Hull Shapefile',multiple=TRUE,accept=c('.shp','.dbf','.sbn','sbx','shx','.prj')),
           fileInput('Ground_Truths_Data','Choose Ground Truth Data Shapefile',multiple=TRUE,accept=c('.shp','.dbf','.sbn','sbx','shx','.prj'))
  ),
  
  #create tab for matching stems game
  tabPanel("Match Stems",
           wellPanel(
             fluidRow(
               #Display of zoom buttons--
               column(1,
                      fluidRow(actionButton('zoomout','-'),actionButton('zoomin','+'),align='left'),offset=6)
               ),
             fluidRow(
               column(12,
                      #Display of map--
                      plotOutput("map",click="click")
               ),
               #Display of iteration controls and confidence value radio buttons
               fluidRow(
                 column(1,actionButton('inv_iterfeature','<<---- Previous Convex Hull Feature'),offset=2),
                 column(2,radioButtons("integer","Confidence: ( 1 Low - High 5)",
                                       c('1'=1,
                                         '2'=2,
                                         '3'=3,
                                         '4'=4,
                                         '5'=5),
                                       inline=TRUE),offset=2),
                 column(1,actionButton('iterfeature','Next Convex Hull Feature ---->>'),position='left',offset=1)),
               fluidRow(
                 column(1,actionButton('nomatch','No Match Present (Skip)'),offset=5)),
               fluidRow(
                 column(2,verbatimTextOutput("info"),offset=5) #small 'quick' information table depicting the selected data
               )
             )
           ),
           #first sidebar panel--
           fluidRow(
             column(10,
                    wellPanel(
                      #create widget in first panel--
                      tableOutput("MatchedTreesTable"),
                      align='center',
                      downloadButton('downloadData','Export Dataframe (.csv)'),
                      useShinyjs()
                    ),offset=1
             )
           )
  )
)

####--------------SERVER-----------------------------------------------------------------
server<-function(input,output,session){

  #Data Importation Controls----
  data<-reactiveValues(importfile=NULL,importfile2=NULL) #Creates reactive value which "holds" the two imported shapefiles
  
  #Chull data importation
  observeEvent(input$LiDAR_Data,{
    shpdf<-input$LiDAR_Data
    tempdir<-dirname(shpdf$datapath[1])
    for(i in 1:nrow(shpdf)){
      file.rename(
        shpdf$datapath[i],
        paste0(tempdir,'/',shpdf$name[i]))
    }
    data[['importfile']]<-readOGR(paste(tempdir,shpdf$name[grep(pattern='*.shp',shpdf$name)],sep='/'))%>%
      spTransform("+proj=utm +zone=17 +ellps=GRS80 +units=m +no_defs")
    data[['importfile']]$tempid<-1:nrow(data[['importfile']])
    data[['importfile']]<-st_as_sf(data[['importfile']])
  })
  #Ground truth data importation
  observeEvent(input$Ground_Truths_Data,{
    shpdf<-input$Ground_Truths_Data
    tempdir<-dirname(shpdf$datapath[1])
    for(i in 1:nrow(shpdf)){
      file.rename(
        shpdf$datapath[i],
        paste0(tempdir,'/',shpdf$name[i]))
    }
    data[['importfile2']]<-readOGR(paste(tempdir,shpdf$name[grep(pattern='*.shp',shpdf$name)],sep='/'))%>%
      spTransform("+proj=utm +zone=17 +ellps=GRS80 +units=m +no_defs")%>%
      st_as_sf()
  })
  
  
  #Feature Iteration Controls----
  featnum<-reactiveVal(0) #sets up a reactive variable to allow for iteratation through rows in DFs
  observeEvent(input$iterfeature,{
    featnum(featnum()+1)
    })
 
  observeEvent(input$inv_iterfeature,{
    if(featnum()-1>0){
      featnum(featnum()-1)#iterates backwards
      tryCatch(silent=TRUE,{
        if(featnum()%in%(df$dt$iteration)){
          if(featnum()-1==0){
            df$dt[1,]<-rbind(NA)
          }
          else(df$dt<-df$dt[-nrow(df$dt[featnum()]),])
        }
        },error=function(e){return()})
    }
    else(return())
    })
  
  observeEvent(input$nomatch,featnum(featnum()+1)) #skips feature if clicks no match present
  
  #Create 'confidence' table to allow user to enter in confidence value
  sliderValues<-reactive({
    data.frame(Name=c("Confidence"),
               Value=as.character(c(input$integer)),
               stringsAsFactors=FALSE)
  })
  output$values<-renderTable({
    sliderValues() #confidence values
  })
  
  
  #Plot & Plot Data Controls----
  #Define reactive values used throughout code
  vals<-reactiveValues(
    chulls=NULL,#initializes these to NULL, but will be updated once running
    truths=NULL,
    buffs=NULL
  )
  
  #Initialize vector of T/F for buffers present in data (only runs when importing grtuth data)
   observeEvent(input$Ground_Truths_Data,
     vals$buffs<-rep(TRUE,nrow(data$importfile2))#creates vector of TRUE for the #rows in the imported ground truth
   )
   
   #Zoom controls
   zoomval<-reactiveVal(4)
   observeEvent(input$zoomin,{
     if(zoomval()>2){zoomval(zoomval()-2)}
     else(zoomval())
   })
   observeEvent(input$zoomout,{
     zoomval(zoomval()+2)
   })
  
  #Create a ggplot 'map' that plots the currently selected feature and ground truth w/ centroids
    output$map<-renderPlot({
      if(featnum()==0) return()#will not display a "map" on first iteration where no features exist
      input$map
      
      allbuffs<-data$importfile2[vals$buffs,,drop=FALSE]#Determines all available gtruth buffers based off of all "TRUE" buffers in vals$buffs vector
      selectedbuffs<-data$importfile2[!vals$buffs,,drop=FALSE]#Determines all "FALSE" buffers based on those not in vals$buffs
      
      map<-ggplot()+
        geom_sf(data=data$importfile[featnum(),],fill='red',alpha=.8)+#add layer for current convex hull
        geom_sf(data=data$importfile%>%st_crop(
          y=c((st_bbox(data$importfile[featnum(),])$xmin)-zoomval()-1, #xmin crop
              (st_bbox(data$importfile[featnum(),])$ymin)-zoomval()-1,#ymin crop
              (st_bbox(data$importfile[featnum(),])$xmax)+zoomval()+1,#xmax crop
              (st_bbox(data$importfile[featnum(),])$ymax)+zoomval()+1))#ymax crop
          ,fill='red',alpha=.2)+ #add layer for 'context' convex hulls (ie not current convex hulls) & crop to speed up loading of map
        geom_sf(data=allbuffs,fill='blue')+#add layer for gtruth buffers
        geom_sf(data=selectedbuffs,fill='green')+#add layer for highlighted (selected) gtruth buffers
        coord_sf(datum=st_crs(data$importfile), #sets projection of the ggplot graph
          xlim=c((st_bbox(data$importfile[featnum(),])$xmin)-zoomval(),(st_bbox(data$importfile[featnum(),])$xmax)+zoomval()), #edits the extents of the graph
          ylim=c((st_bbox(data$importfile[featnum(),])$ymin)-zoomval(),(st_bbox(data$importfile[featnum(),])$ymax)+zoomval()))
        
      map
    })
    

    #Plot Interaction Controls----
    observeEvent(input$click,{ #When clicking on a feature
      clickedbuffs<-nearPoints(data.frame(data$importfile2),input$click,xvar='NADX',yvar='NADY',maxpoints=1,threshold=20,allRows=TRUE)#detect a centerpoint of a gtruth buffer within 20 pixels
      vals$buffs<-rep(TRUE,base::nrow(data$importfile2))
      vals$buffs<-xor(vals$buffs,clickedbuffs$selected_)#then subset the gtruth data based on this selected centroid
      vals$chulls<-data$importfile[featnum(),]#Makes the reactive chulls values under vals= to the subset of the dataframe based on the current iteration being displayed
      vals$truths<-data$importfile2[!vals$buffs,,drop=TRUE]#based on a clicked stem, assign the clicked stem (as a subset from the ground truth) to the reactive value under vals
      })
    
    #Dataframe Controls----
    df<-reactiveValues() #makes a reactive dataframe to store values during matching process
    df$dt<-data.frame(Gtruth_treeID=as.numeric(),
                      Gtruth_stemID=as.numeric(),
                      Gtruth_tag=as.numeric(),
                      Gtruth_StemTag=as.numeric(),
                      Gtruth_dbh_cm=as.numeric(),
                      ULS_ID=as.numeric(),
                      ConvH=as.numeric(),
                      Circle_Pratt=as.numeric(),
                      Circle_LM=as.numeric(),
                      Circle_RANSAC=as.numeric(),
                      Confidence=as.integer(),
                      iteration=as.numeric())
    
    observeEvent(input$iterfeature,{#When the iterfeature button is clicked, and if it is not the first "iteration" (no data present) it will write subsetted reactive data to this row in the dataframe
      try({if(featnum()>1 && !is.na(vals$truths$treeID)){
        row<-data.frame(Gtruth_treeID=vals$truths$treeID,#builds row
                        Gtruth_stemID=vals$truths$stemID,
                        Gtruth_tag=vals$truths$tag,
                        Gtruth_StemTag=vals$truths$StemTag,
                        Gtruth_dbh_cm=vals$truths$dbh_cm,
                        ULS_ID=vals$chulls$StemID,
                        ConvH=vals$chulls$ConvH,
                        Circle_Pratt=vals$chulls$Crcl_Pr,
                        Circle_LM=vals$chulls$Crcl_LM,
                        Circle_RANSAC=vals$chulls$C_RANSA,
                        Confidence=sliderValues()$Value,
                        iteration=featnum()-1)
        df$dt<-rbind(df$dt,row) #binds the row to reactive dataframe
        df$dt<-df$dt[!duplicated(df$dt$iteration),]%>%#removes "duplicate rows"  that are created during this process to ensure only 1 iteration writes 1 row
          filter(!is.na(df$dt))#removes NA rows, if they exist (Primarily used to "realign" the dataframe if the user iterates forwards and backwards multiple times)
        vals$truths<-NULL#resets to NULL when iterating
        }
      },silent=TRUE)#does not return error (when iterating forward and backwards when user does not select buffers and just hits "next")
    })
    output$"MatchedTreesTable"<-renderTable(df$dt) #display the dataframe being built (this will be re-rendered and updated at each iteration)
    
    
    #Misc. tools / features----
    output$info<-renderText({ #output general text for testing purposes primarily
      paste0("Selected Tree ID : x= ",data$importfile2[!vals$buffs,,drop=FALSE]$treeID,
             "\nCurrent Chull Iteration# :",featnum())
    })
    
    output$downloadData<-downloadHandler(
      filename=function(){
        paste("Manually_Matched_Trees_Dataset",".csv",sep="")},
      content=function(file){
        write.csv(df$dt,file,row.names=FALSE)
        })
}


####--------RUN APP-----------------------------------------
shinyApp(ui,server)

exts<-c("C:","Users","SUPER4","Desktop","Apt_APP.pdf")
temp<-paste0(exts,"\\")%>%
  paste(collapse='')
temp<-gsub('\\','/',temp,fixed=TRUE)
temp<-temp[[1]][1:length(temp[[1]])-1]%>%
  paste(collapse="/")
