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
           fluidRow(
             column(12,
                    wellPanel(
                      #Output: Table summarizing the values eneterd--
                      plotOutput("map",click="click"),
                      fluidRow(
                        column(6,actionButton('inv_iterfeature','Previous Convex Hull Feature',position='right'),
                               actionButton('iterfeature','Next Convex Hull Feature',position='left'),offset=6)),
                      fluidRow(
                        column(2,actionButton('nomatch','No Match Present, Skip Feature'),offset=4)),
                      verbatimTextOutput("info")
                      )
                    )
             ),
      #first sidebar panel--
      fluidRow(
        column(8,
               wellPanel(
                 #create widget in first panel--
                 radioButtons("integer","Confidence:",
                              c('1'=1,
                                '2'=2,
                                '3'=3,
                                '4'=4,
                                '5'=5)),
                 tableOutput("MatchedTreesTable"),
                 downloadButton('downloadData','Export Dataframe (.csv)'),
                 useShinyjs()
                 ),offset=2
               )
        )
      )
  )

####--------------SERVER-----------------------------------------------------------------

server<-function(input,output,session){
  data<-reactiveValues(importfile=NULL,importfile2=NULL)
  
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
  
  #Create 'confidence' table to allow user to enter in confidence value----
  sliderValues<-reactive({
    data.frame(Name=c("Confidence"),
               Value=as.character(c(input$integer)),
               stringsAsFactors=FALSE)
  })
  
  #show values in an html table----
  output$values<-renderTable({
    sliderValues() #confidence values
  })
  
  featnum<-reactiveVal(0) #sets up a reactive variable to allow for iteratation through rows in DFs
  observeEvent(input$iterfeature,{
    featnum(featnum()+1)#iterates to next feature
    res<-NULL}) 
  observeEvent(input$inv_iterfeature,{
    if(featnum()-1>0){
      featnum(featnum()-1)#iterates backwards
      if(df$dt$iteration[featnum()]==featnum() && featnum()>1 && nrow(df$dt)>1){
        df$dt<-df$dt[-nrow(df$dt),]
      }
      else({
        df$dt[1,]<-NA
      })
    }
    else(return())
    })
  observeEvent(input$nomatch,featnum(featnum()+1)) #skips feature if clicks no match present
  
  #create reactive groundpoints (to allow for highlighting)
   gtruth<-observeEvent(input$Ground_Truths_Data,
     gtruth$buffs<-rep(TRUE,nrow(data$importfile2))
   )
  
  gtruth<-reactiveValues(buffs=NULL)
   
  #Create a ggplot 'map' that plots the currently selected feature and ground truth w/ centroids
    output$map<-renderPlot({
      if(featnum()==0) return()#will not display a "map" on first iteration where no features exist
      input$map
    
      
      allbuffs<-data$importfile2[gtruth$buffs,,drop=FALSE]#stores within current widget all available gtruth buffers
      selectedbuffs<-data$importfile2[!gtruth$buffs,,drop=FALSE]#stores within current widget all subsetted gtruth buffers
      
      
      map<-ggplot()+
        geom_sf(data=data$importfile[featnum(),],fill='red')+#add layer for convex hulls
        geom_sf(data=allbuffs,fill='blue')+#add layer for gtruth buffers
        geom_sf(data=selectedbuffs,fill='green')+#add layer for highlighted (selected) gtruth buffers
        #geom_point(data=data.frame(gtruthdata),aes(x=gtruthdata$NADY,y=gtruthdata$NADY,colour='pink'))+ #plots centroid (must be clicking within 20pixels of this feature to "match" the gtruth)
        coord_sf(datum=st_crs(data$importfile), #sets projection of the ggplot graph
          xlim=c((st_bbox(data$importfile[featnum(),])$xmin)-2,(st_bbox(data$importfile[featnum(),])$xmax)+2), #edits the extents of the graph
          ylim=c((st_bbox(data$importfile[featnum(),])$ymin)-2,(st_bbox(data$importfile[featnum(),])$ymax)+2))
        
      map
    })

    observeEvent(input$click,{ #When clicking on a feature
      res<-nearPoints(data.frame(data$importfile2),input$click,xvar='NADX',yvar='NADY',maxpoints=1,threshold=20,allRows=TRUE)#detect a centerpoint of a gtruth buffer within 20 pixels
      gtruth$buffs<-rep(TRUE,base::nrow(data$importfile2))
      gtruth$buffs<-xor(gtruth$buffs,res$selected_)#then subset the gtruth data based on this selected centroid
      })
    
    #Build dataframe based off of clicked ground truth stems
    vals<-reactiveValues(
      chulls=NULL,#initializes these to NULL, but will be updated once running
      truths=NULL
    )
    observe({
      vals$chulls<-data$importfile[featnum(),]#Makes the reactive chulls values under vals= to the subset of the dataframe based on the current iteration being displayed
    })
    observeEvent(input$click,{vals$truths<-data$importfile2[!gtruth$buffs,,drop=TRUE]})#based on a clicked stem, assign the clicked stem (as a subset from the ground truth) to the reactive value under vals
    
    df<-reactiveValues() #makes a reactive dataframe to store values during matching process
    df$dt<-data.frame(Gtruth_ID=as.numeric(),
                      Gtruth_SP=as.character(),
                      Gtruth_dbh_cm=as.numeric(),
                      Gtruth_NADX=as.numeric(),
                      GTRUTHNADY=as.numeric(),
                      Chull_diam=as.numeric(),
                      Confidence=as.integer(),
                      iteration=as.numeric())
    
    observeEvent(input$iterfeature,{#When the iterfeature button is clicked, and if it is not the first "iteration" (no data present) it will write subsetted reactive data to this row in the dataframe
      if(featnum()>1 && !is.null(vals$truths)){
        row<-data.frame(Gtruth_ID=vals$truths$treeID,#builds row
                        Gtruth_SP=vals$truths$sp,
                        Gtruth_dbh_cm=vals$truths$dbh_cm,
                        Gtruth_NADX=vals$truths$NADX,
                        Gtruth_NADY=vals$truths$NADY,
                        Chull_diam=vals$chulls$diam,
                        Confidence=sliderValues(),
                        iteration=featnum())
        df$dt<-rbind(df$dt,row) #binds the row to reactive dataframe
        df$dt<-df$dt[!duplicated(df$dt$iteration),]
        #removes "duplicate rows"  that are created during this process to ensure only 1 iteration writes 1 row
      }
    })
    
    output$info<-renderText({ #output general text for testing purposes primarily
      paste0("Selected Tree ID : x=",#gtruthdata[!gtruth$buffs,,drop=FALSE]$treeID, #output window of information
             "\nCurrent Chull Feature # :",featnum(),
             "\nteststuff: ",df$dt$iteration[featnum()-1],
             "\nteststuff2: ")
    })
    
    output$downloadData<-downloadHandler(
      filename=function(){
        paste("Manually_Matched_Trees_Dataset",".csv",sep="")},
      content=function(file){
        write.csv(df$dt,file,row.names=FALSE)
        })
    
    output$"MatchedTreesTable"<-renderTable(df$dt) #display the dataframe being built (this will be re-rendered and updated at each iteration)
}


####--------RUN APP-----------------------------------------
shinyApp(ui,server)

exts<-c("C:","Users","SUPER4","Desktop","Apt_APP.pdf")
temp<-paste0(exts,"\\")%>%
  paste(collapse='')
temp<-gsub('\\','/',temp,fixed=TRUE)
temp<-temp[[1]][1:length(temp[[1]])-1]%>%
  paste(collapse="/")
