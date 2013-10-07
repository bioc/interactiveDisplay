#####################################################################################################
###   Main
#####################################################################################################

## declare the display generic
setGeneric("display", function(object, ...) 
  standardGeneric("display")
)

setMethod("display", 
signature(object = "ANY"), 
function(object){
  message("Wrong object")
})

setMethod("display", 
signature(object = "missing"), 
function(object){
  message("Missing object")
})

#####################################################################################################
###   Helper Functions
#####################################################################################################

## helper for JS library tags

.jstags <- function(){
  list(
  HTML("<script type='text/javascript'>"),
  includeHTML(system.file("www", "jquery.min.js", package="interactiveDisplay")),
  HTML("</script>"),
  
  HTML("<script type='text/javascript'>"),
  includeHTML(system.file("www", "d3.v2.js", package="interactiveDisplay")),
  HTML("</script>"),
  
  HTML("<script type='text/javascript'>"),
  includeHTML(system.file("www", "jquery-svgpan.js", package="interactiveDisplay")),
  HTML("</script>"),
  
  HTML("<script type='text/javascript'>"),
  includeHTML(system.file("www", "jscolor.js", package="interactiveDisplay")),
  HTML("</script>"))
}

#.shiny-output-error { visibility: hidden; }
#.shiny-output-error:before { visibility: hidden; }

.csstags <- function(){

  shiny::tags$head(
    shiny::tags$style(type='text/css', "
  
    .span4 {
      width: 370px;
      position: absolute;
      z-index: 50;
    }
  
    .span8 {
      position: absolute;
      left: 400px;
      right: 30px;
      width: auto;
      height: auto;
    }    

    ")
  )
}

## The loading gif/panel
.loading_gif <- function(){
  list(
  conditionalPanel(condition="$('html').hasClass('shiny-busy')", img(src=system.file("www", "ajax-loader.gif", package="interactiveDisplay"))),
  conditionalPanel(condition="!($('html').hasClass('shiny-busy'))", br())
  )
}


## helper for setting up main panel
.GR_GRL_setMainPanel <- function(){
  mainPanel(  
    tabsetPanel(
      tabPanel("Plot", plotOutput("plotname")),
      tabPanel("All Ranges in Object", uiOutput("fulltable")),
      tabPanel("Selected Ranges in Current View", uiOutput("rtable")),
      tabPanel("Deposited Selections", uiOutput("btable"))
    ),
    uiOutput("mcoltabset")
  )
}


#####################################################################################################
###   Additional Functions
#####################################################################################################


ggheat <- function(my_mat,tweak,colorx,colory,hc,hc2,c1,c2,c3,rainbow){
  melted <- melt(my_mat)
  names(melted) <- c("Var1","Var2","value")
  
  melted$Var1 <- factor(melted$Var1, rownames(my_mat)[hc$order])
  melted$Var2 <- factor(melted$Var2, colnames(my_mat)[hc2$order])
  
  if(length(c1)==0){
    c1 <- "F7FF00"
  }
  if(length(c2)==0){
    c2 <- "050505"
  }
  if(length(c3)==0){
    c3 <- "1736FF"
  }
  
  if(rainbow=='default'){
    myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  }
  
  if(rainbow=='tri'){
    myPalette <- colorRampPalette(rev(c(paste("#",toupper(c1),sep=""),paste("#",toupper(c2),sep=""),paste("#",toupper(c3),sep=""))))
  }
    
  gp <- ggplot(melted, aes(x = Var2, y = Var1, fill = value))#, height = 30*dim(my_mat)[1], width =30*dim(my_mat)[2]))
  gp <- gp + geom_tile()
  gp <- gp + coord_fixed()
  gp <- gp + scale_fill_gradientn(colours = myPalette(100))
  gp <- gp + scale_x_discrete(expand = c(0, 0))
  gp <- gp + scale_y_discrete(expand = c(0, 0))
  gp <- gp + coord_equal()
  gp <- gp + theme_bw()
  gp <- gp + theme(axis.ticks = element_blank(), axis.text.x = element_text(size=tweak*150/sqrt(length(my_mat)), angle = -45, hjust = 0, colour=colorx),
                   axis.text.y = element_text(size=tweak*150/sqrt(length(my_mat)), colour=colory))
  gp <- gp + xlab("Samples")
  gp <- gp + ylab("Probes")
  #gp <- gp + ggtitle("Heatmap") + theme(plot.title = element_text(lineheight=.8, face="bold", vjust = 2))
  gp
}

#####################################################################################################

grid2jssvg <- function(gp){

  jscode <- "
              <script type='text/javascript'>
              $(document).ready(function() {
                $('svg').svgPan('viewport');
              });
              </script>
            "
  png(file = "myplot.png", bg = "transparent",height=1000,width=1000)
  print(gp)
  
  mysvg <- gridSVG::grid.export()
  dev.off()
  mysvg2 <- saveXML(mysvg$svg[["g"]])
  mysvg3 <- sub("<g transform=","<g id='viewport' transform=",mysvg2)
  mysvg4 <- sub(">NA<","><",mysvg3)
  htmlxml <- HTML(paste("<svg xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink' version='1.1' width='100%' height='100%'>",jscode,mysvg4,"</svg>",sep=""))
  htmlxml
}

#####################################################################################################

subgr <- function(gr,chr,strand,window1,window2,width1,width2,mcolnames,input){
  temp1 <- gr[seqnames(gr)==as.character(chr)]
  seqlevels(temp1) <- chr
  if(strand=="both"){
    temp2 <- temp1
  }else{
    temp2 <- temp1[strand(temp1)==strand]
  }
  temp3 <- temp2[ranges(temp2)@width > as.numeric(width1)]
  temp4 <- temp3[ranges(temp3)@width < as.numeric(width2)]
  temp5 <- temp4[start(temp4) > as.numeric(window1)]
  temp6 <- temp5[end(temp5) < as.numeric(window2)]  
  for(i in mcolnames){
    temp6 <- temp6[unlist(as.data.frame(temp6)[i]) %in% input[[i]]]
  }
  temp6
}

#####################################################################################################

subgr2 <- function(gr,chr,strand,width,window,mcolnames,input){
  temp1 <- gr[seqnames(gr)==chr]
  seqlevels(temp1) <- chr
  if(strand=="both"){
    temp2 <- temp1
  }else{
    temp2 <- temp1[strand(temp1)==strand]
  }
  temp3 <- temp2[ranges(temp2)@width > width[1]]
  temp4 <- temp3[ranges(temp3)@width < width[2]]
  temp5 <- temp4[start(temp4) > window[1]]
  temp6 <- temp5[end(temp5) < window[2]]
  for(i in mcolnames){
    temp6 <- temp6[unlist(as.data.frame(temp6)[i]) %in% input[[i]]]
  }
  temp6
}


#####################################################################################################

#  Render the UCSC dropdown
  
.choose_gen <- function(){  
  renderUI({
    ucsc_df <- ucscGenomes()
    ucsc_vec <- as.character(ucsc_df$db)
    names(ucsc_vec) <- ucsc_vec
    selectInput("ucscgen","UCSC Genome",ucsc_vec)
  })
}

#####################################################################################################

#myRunApp <- function(...){
#  try(hostname <- suppressWarnings(system2(c("hostname", "-d"), stdout=TRUE, stderr=NULL)), silent=TRUE)
#  if(exists("hostname") && length(hostname) && grepl("ec2\\.internal", hostname)){
#    dots <- list(...)
#    dots[["launch.browser"]] <- FALSE
#    if (is.null(dots[['port']])) dots[['port']] <- 8100L
#    public.dns <- httr::content(GET("http://169.254.169.254/latest/meta-data/public-hostname"))
#    url <- paste0("http://", public.dns, ":", dots[["port"]])
#    cat("If you don't see a new window, ")
#    cat("try disabling your popup blocker and try again.\n")
#    cat("Press ESCAPE in this window when done.\n")
#    browseURL(url)
#    do.call(runApp, dots)
#  }
#  else{
#    runApp(...)
#  }
#}

#####################################################################################################
###   Common UI elements GRanges/GrangesList
#####################################################################################################

#ui = pageWithSidebar(
#  ####if/else?###headerPanel("Genomic Ranges List"),
#  sidebarPanel(
#    ###if/else?###uiOutput("choose_grange"),
#    uiOutput("choose_chrom"),
#    HTML("<hr />"),
#    uiOutput("choose_gen"),
#    uiOutput("gen_text"),
#    HTML("<hr />"),
#    HTML("Ideogram"),
#    checkboxInput("itr_showId", "Show Id", TRUE),
#    checkboxInput("itr_showBandId", "Show BandId", TRUE),
#    HTML("<hr />"),
#    HTML("Gene Region"),
#    checkboxInput("grtr_showId","Show Id", TRUE),
#    HTML("<hr />"),
#    uiOutput("window"),
#    HTML("<hr />"),
#    uiOutput("width"),
#    HTML("<hr />"),
#    selectInput("strand", "Choose a Strand:", choices = c("both","+","-")),
#    HTML("<hr />"),
#    actionButton("savebutton", "Save to Console"),
#    actionButton("bankbutton", "Deposit Ranges in View"),
#    actionButton("clearbutton", "Clear Deposit")
#  ),
#  mainPanel(  
#    tabsetPanel(
#      tabPanel("Plot", plotOutput("plotname")),
#      tabPanel("All Ranges", uiOutput("fulltable")),
#      tabPanel("Selected Ranges", uiOutput("rtable")),
#      tabPanel("Deposited Ranges", uiOutput("btable"))
#    ),
#    uiOutput("mcoltabset")
#  )
#)



#####################################################################################################
###   Common Functions GRanges/GrangesList
#####################################################################################################

##  The subsetted GRanges object converted to a data frame for the perpose of rendering a table in shiny.
#output$rtable <- renderTable({
#  s_object <- s_object()
#  as.data.frame(s_object)
#})

##  Annotation Track     
#atr <- reactive({
#  s_object <- s_object()
#  AnnotationTrack(s_object, chromosome=input$chr, name="Genomic Ranges Annotation Track", fill="black", background.panel = "#f5f5f5", background.title = "#2aa5e5", cex.title = 1.1)
#})

##  Genome Axis Track
#gtr <- reactive({
#  GenomeAxisTrack(chromosome=input$chr,add53 = TRUE, add35 = TRUE, littleTicks = FALSE)
#})

##  Ideogram Track
#itr <- reactive({
#  IdeogramTrack(genome=input$ucscgen, chromosome=input$chr, showId=input$itr_showId, showBandId=input$itr_showBandId)
#})

##  Render the plot range window slider.     
#output$window <- renderUI({
#  max_end <- max_end()
#  min_start <- min_start()
#  sliderInput(inputId = "window",
#              label = "Plot Range Window:",
#              min = min_start, max = max_end, value = c(min_start,max_end), step = 1
#  )
#})

##  Render the width filter slider.
#output$width <- renderUI({
#  max_width <- max_width()
#  min_width <- min_width()
#  sliderInput(inputId = "width",
#              label = "Range Width Window:",
#              min = min_width, max = max_width, value = c(min_width,max_width), step = 1
#  )
#})

##  Render the track plots.
#output$plotname <- renderPlot({
#  itr <- itr()
#  gtr <- gtr()
#  atr <- atr()
#  plotTracks(list(itr, gtr, atr), from=input$window[1], to=input$window[2])
#})

##  Render the UCSC dropdown
#common_choose_gen <- renderUI({
#  ucsc_df <- ucscGenomes()
#  ucsc_vec <- as.character(ucsc_df$db)
#  names(ucsc_vec) <- ucsc_vec
#  selectInput("ucscgen","UCSC Genome",ucsc_vec)
#})

##  Render the text under the UCSC dropdown        
#common_gen_text <- renderText({
#  ucsc_df <- ucscGenomes()
#  ussc_vec <- as.character(ucsc_df$db)
#  i <- which(ussc_vec==input$ucscgen)
#  paste(as.character(ucscGenomes()[i,2:4]),collapse="&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;")
#})
