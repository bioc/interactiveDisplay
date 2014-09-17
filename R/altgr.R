################################################################################
#  altgr
################################################################################

.selDataTableOutput <- function(outputId, ... ){   
  origStyle <- c( 
    '<script src="js-interactiveDisplay/jquery.dataTables.min.js"></script>',
    '<script class="shiny-html-output" 
                src= "/js-interactiveDisplayBase/DTbinding.js"></script>',
    '<link rel = "stylesheet", 
              type = "text/css", 
              href = "shared/datatables/css/DT_bootstrap.css"></link>',
    '<style type="text/css">
              .rowsSelected td{
              background-color: rgba(112,164,255,0.2) 
              !important})  </style>',
    '<style type="text/css"> .selectable div table tbody tr{
              cursor: hand; cursor: pointer;}</style>',
    '<style type="text/css"> .selectable div table tbody tr td{
              -webkit-touch-callout: none;
              -webkit-user-select: none;
              -khtml-user-select: none;
              -moz-user-select: none;
              -ms-user-select: none;
              user-select: none;} </style>',
    '<style type="text/css">
              #myTable tfoot {display:table-header-group;}</style>')
    
    tagList(
      singleton(
        tags$head(HTML(origStyle)
        )
      ),
      div(id = outputId, class = "shiny-datatable-output selectable")
    )
  }

.altgr <- function(object, ...){
  
  .usePackage('shiny')
  .usePackage('ggbio')
  .usePackage('GenomicRanges')
  
  app <- list(
    
    ui = fluidPage(
      .csstags(),
      absolutePanel(
        top = 40, left = 20, width = 240,
        draggable = TRUE,
        style="padding:8px;border-bottom: 1px solid #CCC; background: #F5F5F5;",
        style = "opacity: 0.90",
        h3("Genomic Ranges", align="center"),
        HTML("<hr />"),
        actionButton("btnSend", "Send Rows"),
        em(p("Shift-Click to select multiple rows.")),
        br(),
        tags$button("Select All Rows", class="btn", id="select_all_rows"),
        em(p("Click to select all rows on page")),
        br(),
        tags$button("Deselect All Rows", class="btn", id="deselect_all_rows"),
        em(p("Click to deselect all rows on page")),
        #br(),
        #selectInput("plotchoice", "Plot:",
        #            c("Circle" = "circle",
        #              "Karyogram" = "karyogram")),
        #em(p("Choose plot type")),
        br(),
        uiOutput("choose_meta"),
        em(p("Choose metadata column for coloration"))
      ),
      .loading_gif(),
      plotOutput("plot1", height="800"),
      #dataTableOutput("mytest"),
      .selDataTableOutput(outputId="myTable", ...)
      ),
    
    server = function(input,output) {
            
      obdf <- as.data.frame(object)
      obnames <- names(obdf)
      r <- 1:dim(obdf)[1]
      obdf <- cbind(r,obdf)
      names(obdf) <- c("row",obnames)
      
      #  Metadata based choices
      output$choose_meta <- renderUI({
        mChoices <- names(elementMetadata(object))
        names(mChoices) <- mChoices
        selectInput("meta", label=NULL, mChoices)
      })

      grdf <- reactive({
        dfVec <- input$myTable
        if(length(dfVec)>9 && length(dfVec)!=0){
          df <- as.data.frame(matrix(data=dfVec,
                                     ncol=dim(obdf)[2],
                                     byrow=TRUE))
          names(df) <- c("row",obnames)
          return(df)
        }
      })
      
      mgr <- reactive({
        df <- grdf()
        if(length(df)!=0){       
          ind <- as.numeric(as.character(df[,1]))
          mgr <- object[ind]
          seqlevels(mgr,force=TRUE) <- sort(unique(as.character((df)$seqnames)))
        }
        mgr
      })
            
      output$plot1 <- renderPlot({
        mgr <- mgr()
        if(length(mgr)>1){
          #if(input$plotchoice=="karyogram"){
            p <- eval(parse(text=paste0("
                          autoplot(mgr,
                          layout='karyogram',
                          aes(color = ",input$meta,",
                              fill = ",input$meta,"))
                      ")))
            return(p)
          #}
          #if(input$plotchoice=="circle"){
          #  p <- ggplot()
          #  p <- p + layout_circle(mgr,
          #                         geom = "ideo",
          #                         radius = 30,
          #                         trackWidth = 4,
          #                         aes(fill=seqnames))
          #  p <- p + layout_circle(mgr,
          #                         geom = "scale",
          #                         size = 2, radius = 35,
          #                         trackWidth = 2)
          #  p <- p + layout_circle(mgr, 
          #                         geom = "rect", 
          #                         color = "steelblue4", 
          #                         radius = 23 , 
          #                         trackWidth = 6)
          #  return(plot(p))
          #}
        }
      })
                        
      output$myTable <- renderDataTable({
        obdf
      })
            
      #output$mytest <- renderDataTable({
      #  df <- grdf()
      #  df
      #})
      
      observe({
        if(input$btnSend > 0)
          isolate({
            mgr <- mgr()
            stopApp(returnValue = mgr)
          })
      })
    
    }
  )
  # runApp(app, ...)
  # selectively use the RStudio viewer pane (if available)
  viewer <- getOption("viewer")
  if (!is.null(viewer)){
    runApp(app, launch.browser = rstudio::viewer, ...)
  }else{
    runApp(app, ...)
  }
}

################################################################################

setGeneric("altgr", function(object, ...)
  standardGeneric("altgr")
)

setMethod("altgr", 
          signature(object = c("ANY")),
          function(object, ...){
            .altgr(object=object, ...)
          })