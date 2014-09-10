################################################################################
#  altgr
################################################################################

df2gr <- function(df){

  gr <- GRanges(seqnames = df$seqnames,
                ranges = IRanges(start = df$start, end = df$end),
                strand = df$strand)

  md <- as(df[, setdiff(names(df),c("seqnames", "start", "end", "width", "strand"))],"DataFrame")
  elementMetadata(gr) <- md
  gr
}

.selDataTableOutput <- 
  function(outputId, ... ) 
  {   
    origStyle<- c( 
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
  
  app <- list(
    
    ui = bootstrapPage(
      .csstags(),
      sidebarPanel(
        h3("Genomic Ranges", align="center"),
        HTML("<hr />"),
        actionButton("btnSend", "Send Rows"),
        em(p("Shift-Click to select multiple rows.")),
        br(),
        tags$button("Select All Rows", class="btn", id="select_all_rows"),
        em(p("Click to select all rows on page")),
        br(),
        tags$button("Deselect All Rows", class="btn", id="deselect_all_rows"),
        em(p("Click to deselect all rows on page"))
      ),
      mainPanel(
      .loading_gif(),
      plotOutput("plot1"),
      #dataTableOutput("mytest"),
      .selDataTableOutput(outputId="myTable", ...)
      )
      ),
    
    server = function(input,output) {

      grdf <- reactive({
        dfVec <- input$myTable
        if(length(dfVec)>9 && length(dfVec)!=0){
          df <- as.data.frame(matrix(data=dfVec,
                                     ncol=dim(as.data.frame(object))[2],
                                     byrow=TRUE))
          df <- cbind(df[,1],apply(df[,(2:4)],2,as.numeric),df[,-(1:4)])
          names(df) <- names(as.data.frame(object))   
          return(df)
        }
      })
      
      mygr <- reactive({
        df <- grdf()
        if(length(df)!=0){
          mygr <- df2gr(df)
          seqlevels(mygr,force=TRUE) <- sort(unique(as.character((df)$seqnames)))
          seqlengths(mygr) <- seqlengths(object)[sort(unique(as.character((df)$seqnames)))]
        }
        mygr
      })
      
      output$plot1 <- renderPlot({
        mygr <- mygr()
        if(length(mygr)>1){
          p <- ggplot()
          p <- p + layout_circle(mygr,
                                 geom = "ideo",
                                 radius = 30,
                                 trackWidth = 4,
                                 aes(fill=seqnames))
          p <- p + layout_circle(mygr,
                                 geom = "scale",
                                 size = 2, radius = 35,
                                 trackWidth = 2)
          p <- p + layout_circle(mygr, 
                                 geom = "rect", 
                                 color = "steelblue4", 
                                 radius = 23 , 
                                 trackWidth = 6)
          plot(p)
        }
      })
            
      output$myTable <- renderDataTable({
        as.data.frame(object)
      })
            
      output$mytest <- renderDataTable({
        df <- grdf()
        df
      })
      
      observe({
        if(input$btnSend > 0)
          isolate({
            mygr <- mygr()
            stopApp(returnValue = mygr)
          })
      })
    
    }
  )
  runApp(app, ...)
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