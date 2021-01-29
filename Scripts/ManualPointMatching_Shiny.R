#to create a interactive graph which will generate a CSV with the points placed in the graph
library(dplyr)
library(raster)
library(ggplot2)
library(shiny)
library(magrittr)
library(DT)
library(rgdal)
library(rgeos)
library(sf)


#perform initial preprocessing of shapefiles to be loaded into app (hardcoded)
chulldata<-readOGR("C:/Users/SUPER4/Desktop/FilesForColin/Calculated_Chulls.shp") %>%
  spTransform("+proj=utm +zone=17 +ellps=GRS80 +units=m +no_defs")
chulldata$centx<-NA
chulldata$centy<-NA
chulldata$tempid<-1:nrow(chulldata)
for(i in 1:nrow(chulldata)){
  chulldata$centx[[i]]<-as.numeric(data.frame(gCentroid(chulldata))$x)
  chulldata$centy[[i]]<-as.numeric(data.frame(gCentroid(chulldata))$y)
}
chulldata<-st_as_sf(chulldata)
gtruthdata<-readOGR("C:/Users/SUPER4/Desktop/FilesForColin/GroundTruth_DBH_StemBuffers - Copy.shp")%>%
  spTransform("+proj=utm +zone=17 +ellps=GRS80 +units=m +no_defs")
gtruthcents<-gCentroid(gtruthdata,byid=TRUE)
for(i in 1:length(gtruthcents)){
  gtruthcents$xcoord[i]<-data.frame(gtruthcents[i])$x
  gtruthcents$ycoord[i]<-data.frame(gtruthcents[i])$y
}
gtruthdata<-st_as_sf(gtruthdata)
####-------------UI-----------------------------------------------------------------------
ui<-navbarPage(
  #Title for app----
  "Match The Convex Hull With The Respective Ground Truth 'Game'",
  
  #create second tab in UI for data initialization
  tabPanel("Load Data",
           
           #chull UI element
           fileInput('shapefile','Choose LiDAR Convex Hulls',multiple=FALSE,accept="shp"),
           fluidRow(column(4,verbatimTextOutput("chull_fname")))
  ),
  
  #create tab for matching stems game
  tabPanel("Match Stems",
      #first sidebar panel----
      sidebarPanel(
        #create widget in first panel----
        radioButtons("integer","Confidence:",
                     c('1'=1,
                       '2'=2,
                       '3'=3,
                       '4'=4,
                       '5'=5)),
        actionButton('iterfeature','Next Convex Hull Feature'),
        actionButton('inv_iterfeature','Previous Convex Hull Feature'),
        actionButton('nomatch','No Match Present, Skip Feature')
        ),
      #Main panel for displaying outputs----
      mainPanel(
        #Output: Table summarizing the values eneterd----
        plotOutput("map",click="click"),
        verbatimTextOutput("info"),
        tableOutput("MatchedTreesTable"),
        downloadButton('downloadData','Export Dataframe (.csv)')
        )
    )
)

####--------------SERVER-----------------------------------------------------------------

server<-function(input,output,session){
  #Reactive expression to create data frame of allinput values----
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
  observeEvent(input$iterfeature,featnum(featnum()+1)) #iterates to next feature
  observeEvent(input$inv_iterfeature,{
               featnum(featnum()-1)
               df$dt<-df$dt[-nrow(df$dt),]})#iterates backwards
  observeEvent(input$nomatch,featnum(featnum()+1)) #skips feature if clicks no match present
  
  
    #selectedfeat$featid<-selectedfeat$featcent$stemID
  
  #create reactive groundpoints (to allow for highlighting)
  gtruth<-reactiveValues(
    buffs=rep(TRUE,nrow(gtruthdata))
  )
  
  
  #Create a ggplot 'map' that plots the currently selected feature and ground truth w/ centroids
    output$map<-renderPlot({
      if(featnum()==0) return()#will not display a "map" on first iteration where no features exist
      input$map
      
      allbuffs<-gtruthdata[gtruth$buffs,,drop=FALSE]#stores within current widget all available gtruth buffers
      selectedbuffs<-gtruthdata[!gtruth$buffs,,drop=FALSE]#stores within current widget all subsetted gtruth buffers
      
      
      map<-ggplot()+
        geom_sf(data=chulldata[featnum(),],fill='red')+#add layer for convex hulls
        geom_sf(data=allbuffs,fill='blue')+#add layer for gtruth buffers
        geom_sf(data=selectedbuffs,fill='green')+#add layer for highlighted (selected) gtruth buffers
        geom_point(data=data.frame(gtruthcents),aes(x=data.frame(gtruthcents@coords)$x,y=data.frame(gtruthcents@coords)$y,colour='pink'))+ #plots centroid (must be clicking within 20pixels of this feature to "match" the gtruth)
        coord_sf(datum=st_crs(chulldata), #sets projection of the ggplot graph
          xlim=c((st_bbox(chulldata[featnum(),])$xmin)-2,(st_bbox(chulldata[featnum(),])$xmax)+2), #edits the extents of the graph
          ylim=c((st_bbox(chulldata[featnum(),])$ymin)-2,(st_bbox(chulldata[featnum(),])$ymax)+2))
      map
    })

    observeEvent(input$click,{ #When clicking on a feature
      res<-nearPoints(data.frame(gtruthcents),input$click,xvar='x',yvar='y',maxpoints=1,threshold=20,allRows=TRUE)#detect a centerpoint of a gtruth buffer within 20 pixels
      gtruth$buffs<-xor(gtruth$buffs,res$selected_)#then subset the gtruth data based on this selected centroid
    })
    
    
    
    #Build dataframe based off of clicked ground truth stems
    vals<-reactiveValues(
      chulls=NULL,#initializes these to NULL, but will be updated once running
      truths=NULL
    )
    observe({
      vals$chulls<-chulldata[featnum(),]#Makes the reactive chulls values under vals= to the subset of the dataframe based on the current iteration being displayed
    })
    observeEvent(input$click,{vals$truths<-gtruthdata[!gtruth$buffs,,drop=TRUE]})#based on a clicked stem, assign the clicked stem (as a subset from the ground truth) to the reactive value under vals
    
    df<-reactiveValues() #makes a reactive dataframe to store values during matching process
    df$dt<-data.frame(Gtruth_ID=as.numeric(),
                      Gtruth_SP=as.character(),
                      Gtruth_dbh_cm=as.numeric(),
                      Gtruth_NADX=as.numeric(),
                      GTRUTHNADY=as.numeric(),
                      Chull_diam=as.numeric(),
                      Chull_centx=as.numeric(),
                      Chull_centy=as.numeric(),
                      Confidence=as.integer(),
                      iteration=as.numeric())
    
    observeEvent(input$iterfeature,{#When the iterfeature button is clicked, and if it is not the first "iteration" (no data present) it will write subsetted reactive data to this row in the dataframe
      if(featnum()>1){
        row<-data.frame(Gtruth_ID=vals$truths$treeID,#builds row
                        Gtruth_SP=vals$truths$sp,
                        Gtruth_dbh_cm=vals$truths$dbh_cm,
                        Gtruth_NADX=vals$truths$NADX,
                        Gtruth_NADY=vals$truths$NADY,
                        Chull_diam=vals$chulls$diam,
                        Chull_centx=vals$chulls$centx,
                        Chull_centy=vals$chulls$centy,
                        Confidence=sliderValues(),
                        iteration=featnum())
        df$dt<-rbind(df$dt,row) #binds the row to reactive dataframe
        df$dt<-df$dt[!duplicated(df$dt$iteration),]#removes "duplicate rows"  that are created during this process to ensure only 1 iteration writes 1 row
      }
    })
    
    
    
    output$info<-renderText({ #output general text for testing purposes primarily
      paste0("Selected Tree ID : x=",gtruthdata[!gtruth$buffs,,drop=FALSE]$treeID, #output window of information
             "\nCurrent Chull Feature # :",featnum(),
             "\nteststuff: ", vals$chulls$tempid,
             "\nteststuff2: ",vals$truths$treeID)
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
