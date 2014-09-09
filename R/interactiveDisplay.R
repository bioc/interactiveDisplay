################################################################################
###   Main
################################################################################

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

################################################################################
###   Helper Functions
################################################################################

source("http://bioconductor.org/biocLite.R")
.usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])){
    stop(paste("The required package, '",p,"', is missing.  Please install it by
               typing biocLite('",p,"') in the console", sep=""))
  }
  require(p, character.only = TRUE)
}


## helper for JS library tags

.jstags <- function(){
  list(
  tags$script(src="/js-interactiveDisplay/jquery.min.js"),
  tags$script(src="/js-interactiveDisplay/d3.v2.js"),
  tags$script(src="/js-interactiveDisplay/jquery-svgpan.js"),
  tags$script(src="/js-interactiveDisplay/jscolor/jscolor.js"))
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
  conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                   div("Loading...", style = "color:blue")),
  conditionalPanel(condition="!($('html').hasClass('shiny-busy'))", br())
  )
}

#selDataTableOutput <- function (outputId){
#  tagList(singleton(tags$head(tags$link(rel = "stylesheet", 
#    type = "text/css", href = "shared/datatables/css/DT_bootstrap.css"),
#    tags$style(type="text/css", ".rowsSelected td{background-color: rgba(112,164,255,0.2) !important}"),
#    tags$style(type="text/css", ".selectable div table tbody tr{cursor: hand; cursor: pointer;}"),
#    tags$style(type="text/css",".selectable div table tbody tr td{
#      -webkit-touch-callout: none;
#      -webkit-user-select: none;
#      -khtml-user-select: none;
#      -moz-user-select: none;
#      -ms-user-select: none;
#      user-select: none;}"),                          
#    tags$script(src = "shared/datatables/js/jquery.dataTables.min.js"), 
#    tags$script(src = "shared/datatables/js/DT_bootstrap.js"),
#    tags$script(src = "/js-interactiveDisplay/DTbinding.js"))),
#  div(id = outputId, class = "shiny-datatable-output selectable"))
#}

## helper for setting up main panel
.GR_GRL_setMainPanel <- function(sflag){
  mainPanel(
    shiny::tags$head(
      shiny::tags$style(type='text/css', "

    cplot {
      height: 800px;
    }

    svg {
      height: 800px;
    }

    ")
    ),
    .loading_gif(),
    tabsetPanel(
      tabPanel("Static Circle Layout",
               HTML("Use the mouse to drag and pan the circle layout.  Use the 
                     mousewheel to zoom in/out."),
               HTML("<hr />"),
               svgcheckout("cplot",sflag)),
      tabPanel("Interactive Plot", plotOutput("plotname")),
      tabPanel("All Ranges in Object", dataTableOutput("fulltable")),
      tabPanel("Selected Ranges in Current View", dataTableOutput("rtable")),
      tabPanel("Deposited Selections", dataTableOutput("btable"))
    ),
    uiOutput("mcoltabset")
  )
}

################################################################################
###   Additional Functions
################################################################################

ggheat <- function(my_mat,
                   tweak,
                   color_samples,
                   color_probes,
                   hc,
                   hc2,
                   c1,
                   c2,
                   c3,
                   rainbow,
                   flip){
  melted <- melt(my_mat)
  names(melted) <- c("Var1","Var2","value")
  
  melted$Var1 <- factor(melted$Var1, rownames(my_mat)[hc$order])
  melted$Var2 <- factor(melted$Var2, colnames(my_mat)[hc2$order])
  
  if(length(c1)==0){
    c1 <- "EDF8B1"
  }
  if(length(c2)==0){
    c2 <- "7FCDBB"
  }
  if(length(c3)==0){
    c3 <- "2C7FB8"
  }
  
  if(rainbow=='default'){
    myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  }
  
  if(rainbow=='tri'){
    myPalette <- colorRampPalette(rev(c(paste("#",toupper(c1),sep=""),
      paste("#",toupper(c2),sep=""),paste("#",toupper(c3),sep=""))))
  }
    
  gp <- ggplot(melted, aes(x = Var2, y = Var1, fill = value))
  gp <- gp + geom_tile()
  gp <- gp + coord_fixed()
  gp <- gp + scale_fill_gradientn(colours = myPalette(100))
  gp <- gp + scale_x_discrete(expand = c(0, 0))
  gp <- gp + scale_y_discrete(expand = c(0, 0))
  gp <- gp + coord_equal()
  gp <- gp + theme_bw()
  gp <- gp + xlab("Samples")
  gp <- gp + ylab("Probes")
   
  if(flip==TRUE){
    gp <- gp + theme(axis.ticks = element_blank(),
                     axis.text.x = element_text(size = 
                                                tweak*150/sqrt(length(my_mat)),
                                                angle = -45,
                                                hjust = 0,
                                                colour=color_probes),
                     axis.text.y = element_text(size = tweak*150/sqrt(length(my_mat)),
                                                colour = color_samples))
    gp <- gp + coord_flip() 
  }
  else{
    gp <- gp + theme(axis.ticks = element_blank(),
                     axis.text.x = element_text(size = 
                                                tweak*150/sqrt(length(my_mat)),
                                                angle = -45,
                                                hjust = 0,
                                                colour=color_samples),
                     axis.text.y = element_text(size = tweak*150/sqrt(length(my_mat)),
                                                colour = color_probes))
  }
  
  return(gp)
}

################################################################################

grid2jssvg <- function(gp){

  jscode <- "
              <script type='text/javascript'>
              $(document).ready(function() {
                $('svg').svgPan('viewport');
              });
              </script>
            "
  png(filename = "myplot.png", bg = "transparent",height=1000,width=1000)
  print(gp)
  
  mysvg <- gridSVG::grid.export()
  dev.off()
  mysvg2 <- saveXML(mysvg$svg[["g"]])
  mysvg3 <- sub("<g transform=","<g id='viewport' transform=",mysvg2)
  mysvg4 <- sub(">NA<","><",mysvg3)
  htmlxml <- HTML(paste("<svg xmlns='http://www.w3.org/2000/svg' 
  xmlns:xlink='http://www.w3.org/1999/xlink' version='1.1' width='100%' 
  height='100%'>",jscode,mysvg4,"</svg>",sep=""))
  htmlxml
}

#  This pair of functions can be used in cases where it is desirable to
#  give the user a choice between rendering a plot as svg or to use the default
#  Shiny plot function.

svgcheckout <- function(contents,sflag){
  if(sflag==TRUE){
    uiOutput(contents)
  }
  else{
    plotOutput(contents)
  }
}

################################################################################

subgr <- function(gr,chr,strand,window1,window2,width1,width2,mcolnames,input){
  temp1 <- gr[seqnames(gr)==as.character(chr)]
  seqlevels(temp1) <- chr
  if(strand=="both"){
    temp2 <- temp1
  }else{
    temp2 <- temp1[strand(temp1)==strand]
  }
  if(!is.null(width1) && !is.null(width2)){
    temp3 <- temp2[ranges(temp2)@width > as.numeric(width1)]
    temp4 <- temp3[ranges(temp3)@width < as.numeric(width2)]
  }
  else{
    temp4 <- temp2
  }
  temp5 <- temp4[start(temp4) > as.numeric(window1)]
  temp6 <- temp5[end(temp5) < as.numeric(window2)]  
  for(i in mcolnames){
    temp6 <- temp6[unlist(as.data.frame(temp6)[i]) %in% input[[i]]]
  }
  temp6
}

################################################################################

subgr2 <- function(gr,chr,strand,width,window,mcolnames,input){
  temp1 <- gr[seqnames(gr)==chr]
  seqlevels(temp1) <- chr
  if(strand=="both"){
    temp2 <- temp1
  }else{
    temp2 <- temp1[strand(temp1)==strand]
  }
  if(!is.null(width)){
    temp3 <- temp2[ranges(temp2)@width > width[1]]
    temp4 <- temp3[ranges(temp3)@width < width[2]]
  }
  else{
    temp4 <- temp2
  }
  temp5 <- temp4[start(temp4) > window[1]]
  temp6 <- temp5[end(temp5) < window[2]]
  for(i in mcolnames){
    temp6 <- temp6[unlist(as.data.frame(temp6)[i]) %in% input[[i]]]
  }
  temp6
}

################################################################################

#  Render the UCSC dropdown
  
.choose_gen <- function(){  
  renderUI({
    ucsc_df <- ucscGenomes()
    ucsc_vec <- as.character(ucsc_df$db)
    names(ucsc_vec) <- ucsc_vec
    selectInput("ucscgen","UCSC Genome",ucsc_vec)
  })
}

