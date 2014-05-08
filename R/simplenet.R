################################################################################
#  Template
################################################################################

.simplenet <- function(object, ...){
  
  .usePackage('shiny')
  
  app <- list(
    ui =
      bootstrapPage(
        
        sidebarPanel(
          checkboxInput("transpose","Transpose Data"),
          HTML("<hr />"),
          sliderInput(inputId = "charge",
                      label = "Force Layout Charge",
                      min = -2000, max = -10, value = -800, step = 1),
          sliderInput(inputId = "linkDistance",
                      label = "Force Layout Link Distance",
                      min = 10, max = 200, value = 80, step = 1),
          HTML("<hr />"),
          numericInput("con_knum", "Number of Clusters:", 2),
          numericInput("edgenum", "Number of Edges:", 1),
          uiOutput("gen_text"),
          HTML("<hr />"),
          selectInput("hc_method", "Hierarchical Clustering Method",
                      choices = c("ward", "single", "complete", "average", 
                                  "mcquitty", "median", "centroid")),
          selectInput("dist_method", "Distance/Similarity Method",
                      choices = c("euclidean", "maximum", "manhattan", 
                                  "canberra", "minkowski")),
          HTML("<hr />"),
          actionButton("closebutton", "Stop Application")
        ),
        
        mainPanel(
          tags$link(rel="stylesheet", type="text/css",
                    href="/css/interactiveDisplay.css"),
          tags$script(src="http://d3js.org/d3.v2.js"),
          tags$script(src="/js/graph.js"),
          .jstags(),
          .csstags(),
          shiny::tags$head(
            shiny::tags$style(type='text/css', "
              svg {
                height: 95vh;
              }
            ")
          ),
          .loading_gif(),
          uiOutput("svg")
        )
      ),
    
    server = function(input,output) {
      
      tmpdata <- reactive({object})
      
      #  Build the network
      output$net <- reactive({
        data <- data()
        hc <- hc(data)
        data <- data[,hc$order]
        cutoff <- cutoff()
        if(length(data)!=0){
          val <- dm(data())
          if (is.null(val)){
            return(list(names=character(), links=list(source=-1, target=-1)))
          }
          diag(val) <- NA
          val[lower.tri(val)] <- NA
          if(sum(cutoff)==0){
            val[] <- NA
          }else{
            val[val > cutoff] <- NA
          }
          conns <- cbind(source=row(val)[!is.na(val)]-1,
                         target=col(val)[!is.na(val)]-1,
                         weight=val[!is.na(val)])
          if (nrow(conns) == 0){
            conns <- list(source=-1, target=-1, weight=0)
          }
          net <- list()
          net[["names"]] <- colnames(val)
          net[["links"]] <- conns
          net[["groups"]] <- as.numeric(cutree(hc, k=input$con_knum))
          net[["titles"]] <- hc$labels[hc$order]
          net[["colors"]] <- 
            rainbow(input$con_knum,
                    alpha=NULL)[cutree(hc,input$con_knum)[colnames(val)]]
          net[["charge"]] <- input$charge 
          net[["linkDistance"]] <- input$linkDistance
          
          return(net)
        }
        else{
          return()
        }
      })
      
      #  Data for network view, sample or probe
      data <- reactive({
        data <- tmpdata()
        if(length(data)!=0){
          if(input$transpose=="TRUE"){
            data <- data
          }else{
            data <- t(data)
          }
          return(data)
        }
        else{
          return()
        }
      })
      
      #  This determines the distance threshold needed for the desired
      #  number of edges
      cutoff <- reactive({
        data <- data()
        if(length(data)!=0){
          val <- dm(data())
          diag(val) <- NA
          val[lower.tri(val)] <- NA
          cutoff <- sort(val[!is.na(val)],decreasing=FALSE)[input$edgenum]
          if(isNA(cutoff || is.null(cutoff))){
            cutoff <- 0
          }
          return(cutoff)
        }
        else{
          return()
        }
      })

      #  Show Distance Threshold
      output$gen_text <- renderText({
        cutoff <- cutoff()
        if(length(cutoff)!=0){
          paste("Distance Threshold:  ",round(cutoff,4),sep="")
        }
        else{
          return()
        }
      })
      
      #  Clustering
      hc <- function(d){
        hclust(dist(t(d),method = input$dist_method), input$hc_method)
      }
      
      #  Distance Matrix
      dm <- function(d){
        as.matrix(dist(t(d), diag=TRUE, upper=TRUE, method=input$dist_method))
      }
      
      #  The network SVG
      output$svg <- renderUI({
        HTML(paste(
          "<div id=\"net\" class=\"shiny-network-output\"><svg /></div>",
          sep=""))
      })
      
      #  Close Button  
      observe({
        if (input$closebutton == 0)
          return()
        isolate({
          stopApp()
        })
      })
      
    }
  )
  runApp(app, ...)
}

################################################################################

setGeneric("simplenet", function(object, ...)
  standardGeneric("simplenet")
)

setMethod("simplenet", 
          signature(object = c("ANY")),
          function(object, ...){
            .simplenet(object=object, ...)
          })