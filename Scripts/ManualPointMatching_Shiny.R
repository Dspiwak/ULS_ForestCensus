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
           fileInput('chapefile','Choose LiDAR Convex Hulls',multiple=FALSE,accept="shp"),
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
        tableOutput("MatchedTreesTable")
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
    sliderValues() 
  })
  
  featnum<-reactiveVal(0) #sets up a reactive variable to allow for iteratation through rows in DFs
  observeEvent(input$iterfeature,featnum(featnum()+1)) #iterates to next feature
  observeEvent(input$inv_iterfeature,featnum(featnum()-1))
  observeEvent(input$nomatch,featnum(featnum()+1)) #skips feature if clicks no match present
  
  
    #selectedfeat$featid<-selectedfeat$featcent$stemID
  
  #create reactive groundpoints (to allow for highlighting)
  gtruth<-reactiveValues(
    buffs=rep(TRUE,nrow(gtruthdata))
  )
  
  
  #Create a ggplot 'map' that plots the currently selected feature and ground truth w/ centroids
    output$map<-renderPlot({
      if(featnum()==0) return()
      input$map
      
      allbuffs<-gtruthdata[gtruth$buffs,,drop=FALSE]
      selectedbuffs<-gtruthdata[!gtruth$buffs,,drop=FALSE]
      
      
      map<-ggplot()+
        geom_sf(data=chulldata[featnum(),],fill='red')+
        geom_sf(data=allbuffs,fill='blue')+
        geom_sf(data=selectedbuffs,fill='green')+
        geom_point(data=data.frame(gtruthcents),aes(x=data.frame(gtruthcents@coords)$x,y=data.frame(gtruthcents@coords)$y,colour='pink'))+
        coord_sf(datum=st_crs(chulldata))+ #sets projection of the ggplot graph
        xlim((st_bbox(chulldata[featnum(),])$xmin)-2,(st_bbox(chulldata[featnum(),])$xmax)+2)+ #edits the extents of the graph
        ylim((st_bbox(chulldata[featnum(),])$ymin)-2,(st_bbox(chulldata[featnum(),])$ymax)+2)
      map
    })

    observeEvent(input$click,{
      res<-nearPoints(data.frame(gtruthcents),input$click,xvar='x',yvar='y',maxpoints=1,threshold=20,allRows=TRUE)
      gtruth$buffs<-xor(gtruth$buffs,res$selected_)
    })
    
    vals<-reactiveValues(
      chulls=NULL,
      truths=NULL
    )
    observe({
      vals$chulls<-chulldata[featnum(),]$diam
    })
    observeEvent(input$click,{vals$truths<-gtruthdata[!gtruth$buffs,,drop=TRUE]$treeID})
    
    df<-reactiveValues()
    df$dt<-data.frame(ID=as.numeric(),
                      chulldiam=as.numeric(),
                      iteration=as.numeric())
    
    observeEvent(input$iterfeature,{
      if(featnum()>1){
        row<-data.frame(ID=vals$truths,
                        chulldiam=vals$chulls,
                        iteration=featnum())
        df$dt<-rbind(df$dt,row)
        df$dt<-df$dt[!duplicated(df$dt$iteration),]
      }
    })
    
    
    
    output$info<-renderText({
      paste0("Selected Tree ID : x=",gtruthdata[!gtruth$buffs,,drop=FALSE]$treeID, #output window of information
             "\nCurrent Chull Feature # :",featnum(),
             "teststuff: ", vals$chulls,
             "teststuff2: ",vals$truths)
    })
    
    output$"MatchedTreesTable"<-renderTable(df$dt)
    
}
####--------RUN APP-----------------------------------------
shinyApp(ui,server)


temp<-(read.csv("C:/Users/note2/Downloads/Generated (2).csv"))
testemp<-mutate(temp,plotval=1)
spdf<-SpatialPoints(testemp[,1:2],data=as.data.frame(testemp$plotval))
ras<-raster(spdf)
plot(ras)
plot(spdf)