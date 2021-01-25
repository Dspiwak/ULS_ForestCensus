#to create a interactive graph which will generate a CSV with the points placed in the graph
library(dplyr)
library(raster)
library(ggplot2)
library(shiny)
library(magrittr)
library(DT)
library(rgdal)
library(leaflet)
library(rgeos)



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
        actionButton('nomatch','No Match Present, Skip Feature')
        ),
      #Main panel for displaying outputs----
      mainPanel(
        #Output: Table summarizing the values eneterd----
        plotOutput("map",click="click"),
        verbatimTextOutput("info")
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
  
  epsg6346<-leafletCRS(
    crsClass="L.Proj.CRS",
    code="EPSG:6346",
    proj4def="+proj=utm +zone=17 +datum=GRS80 +units=m +no_defs",
    resolutions = 2^(16:7)
  )
  
  #add data to table
  # output$map<-renderLeaflet({
  #   leaflet(options=leafletOptions(nowrap=TRUE, crs=epsg6346))%>%
  #     addPolygons(data=chulldata,fill=TRUE,color='blue')%>%
  #     addPolygons(data=testarea,fill=TRUE,color='red')
  #   })
  
  
  featnum<-reactiveVal(0) #sets up a reactive variable to allow for iteratation through rows in DFs
  observeEvent(input$iterfeature,featnum(featnum()+1)) #iterates to next feature
  observeEvent(input$nomatch,featnum(featnum()+1)) #skips feature if clicks no match present
  
  #Create a ggplot 'map' that plots the currently selected feature and ground truth w/ centroids
    output$map<-renderPlot({
      if(featnum()==0) return()
      #input$map
      input$map
      ggplot()+
        geom_sf(data=chulldata[featnum(),],fill='red')+
        geom_sf(data=gtruthdata,fill='blue')+
        geom_point(data=data.frame(gtruthcents),aes(x=data.frame(gtruthcents@coords)$x,y=data.frame(gtruthcents@coords)$y,colour='pink'))+
        coord_sf(datum=st_crs(chulldata))+ #sets projection of the ggplot graph
        xlim((st_bbox(chulldata[featnum(),])$xmin)-2,(st_bbox(chulldata[featnum(),])$xmax)+2)+ #edits the extents of the graph
        ylim((st_bbox(chulldata[featnum(),])$ymin)-2,(st_bbox(chulldata[featnum(),])$ymax)+2)
      #plot(chulldata)
    })
    
    #ooga<-data.frame(sf::st_coordinates(chulldata))
    output$info<-renderText({
      #paste0(nearPoints(data.frame(chulldata),input$click,'centx','centy'))
      paste0("Nearest Point :",nearPoints(data.frame(gtruthcents),input$click,xvar='x',yvar='y',maxpoints=1,threshold=20),"\nFeature # :",featnum())
    })
}
####--------RUN APP-----------------------------------------
shinyApp(ui,server)


temp<-(read.csv("C:/Users/note2/Downloads/Generated (2).csv"))
testemp<-mutate(temp,plotval=1)
spdf<-SpatialPoints(testemp[,1:2],data=as.data.frame(testemp$plotval))
ras<-raster(spdf)
plot(ras)
plot(spdf)