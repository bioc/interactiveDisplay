#####################################################################################################
###   SummarizedExperiment
#####################################################################################################

setMethod("display", 
  signature(object = "SummarizedExperiment"), 
  function(object, ...){
      require(org.Hs.eg.db)
      require(org.Mm.eg.db)
      require(org.Dm.eg.db) ## TODO: We want something more elegant that this
      map <- list(Mouse=org.Mm.eg.db, Human=org.Hs.eg.db, Fly=org.Dm.eg.db)
      app <- list(
      ui = pageWithSidebar(
        
        headerPanel("Summarized Experiment Object"),

        sidebarPanel(
          includeHTML(system.file("www", "graph.js", package="interactiveDisplay")),
          textInput("padjVal","Adjusted p-value cutoff", value=.01), 
          textInput("fcVal","log 2 fold change", value=2),  
          
          ## TODO: organism/annot should be stored/retrieved from object.
          ## Not hard coded as in this example
          
          selectInput("columns", "Annotations", 
            choices = columns(map[["Fly"]]),
            selected = c("SYMBOL", "GENENAME"),
            multiple=TRUE)
        ),
        
        mainPanel(
          tabsetPanel(
            tabPanel("HeatMap", plotOutput("plot")),         
            tabPanel("Genes",tableOutput("view"))
          )
          
        )
      ),
      
      server = function(input, output){
        
        #########################################################
        ## Some helper functions:
        getSec <- reactive({
          ##      load("se.rda")
          ## get rowData
          ## se <- .secretSEVariable
            se <- object
          rd <- rowData(se)
          ## throw out bad data
          se[!is.na(rd$padj),]
        })
        
        getRowDataAtCutoffs <- reactive({
          padjVal <- input$padjVal
          fcVal <- input$fcVal
          sec <- getSec()
          require(GenomicRanges)
          ## get rowData again (clean this time)
          rd <- rowData(sec)
          rd <- rd[rd$padj < padjVal]
          rd <- rd[rd$foldChange >= fcVal]
        })
        
        getIds <- reactive({
          padjVal <- input$padjVal
          fcVal <- input$fcVal
          sec <- getSec()
          getRowDataAtCutoffs()$id
        })
        
        makeColColors <- reactive({
          sec <- getSec()
          f <- colData(sec)[[1]] ## in future, this must be more explicit?
          conds <- levels(f)
          ccs <- rainbow(length(conds))
          map <- cbind(conds,ccs)
          map[match(f, map[,1]),2]
        })
        
        
        ###########################################################
        ## active plots
        
        output$plot <- renderPlot({
          sec <- getSec()
          ids <- getIds()
          
          cc <- makeColColors()
          
          ## heatmap(assays(sec)[[1]][ids,], ColSideColors=cc) #, srt=45)        
          ## ,labRow=character(),
          ## labCol=character())
          
          ## better heatmap?
          require(gplots)
          heatmap.2(assays(sec)[[1]][ids,], ColSideColors=cc, scale="row") 
          
        }, height="auto")
        
        output$view <- renderTable({        
          rd <- getRowDataAtCutoffs()
          ids <- rd$id
          
          cols <- input$columns
          if (!length(cols))
            cols <- head(columns(map[["Fly"]]), 1L)
          
          annot <- select(map[["Fly"]], keys=ids, columns=cols,
                                         keytype="FLYBASE")
          mcs <- mcols(rd)
          mcs <- mcs[, !(colnames(mcs) %in% c("baseMean","foldChange")) ]
          as.data.frame(merge(mcs, annot, by.x = "id", by.y="FLYBASE"))    
        })
      }
    )
    #myRunApp(app, ...)
    runApp(app)
    
  })

