################################################################################
###   GRanges
################################################################################

## helper for setting up sidebar
.setSidebarPanel <- function(){
  sidebarPanel(
    uiOutput("choose_chrom"),
    HTML("<hr />"),
    uiOutput("choose_gen"),
    uiOutput("gen_text"),
    HTML("<hr />"),
    uiOutput("window"),
    HTML("<hr />"),
    uiOutput("width"),
    HTML("<hr />"),
    selectInput("strand", "Choose a Strand:", choices = c("both","+","-")),
    HTML("<hr />"),
    actionButton("bankbutton", "Deposit Ranges in View"),
    HTML("<hr />"),
    actionButton("clearbutton", "Clear Deposit"),
    actionButton("savebutton", "Save to Console")
  )
}

setMethod("display", 
  signature(object = "GRanges"), 
  function(object, ...){
    
    mcolnames <- as.character(names(mcols(object)))
    
    app <- list(
      ui = bootstrapPage(        
        .jstags(),  
        .csstags(),
        h3("Genomic Ranges"),
        .loading_gif(),
        .setSidebarPanel(),
        .GR_GRL_setMainPanel()
      ),
      
      server = function(input, output){
        
        #  This stores parameters for subsetted GRanges per chromosome as a list.
        bank <- list()
        
        t_object <- reactive({
          object[[input$gr]] 
        })
        
        #  This GRanges object is subsetted by user input in the shiny widget.
        s_object <- reactive({
          subgr2(object,input$chr,input$strand,input$width,input$window,
            mcolnames,input)
        })
        
        #  The full submitted GRanges object converted to a data frame for the
        #  purpose of rendering a table in shiny. 
        output$fulltable <- renderTable({
          as.data.frame(object)
        })
        
        #  The subsetted GRanges object converted to a data frame for the
        #  purpose of rendering a table in shiny.
        output$rtable <- renderTable({
          s_object <- s_object()
          as.data.frame(s_object)
        })
        
        #  Annotation Track     
        atr <- reactive({
          s_object <- s_object()
          if(!is.null(input$chr)){
            AnnotationTrack(s_object, chromosome=input$chr,
            name="Genomic Ranges Annotation Track", fill="black",
            background.panel = "#f5f5f5", background.title = "#2aa5e5",
            cex.title = 1.1)
          }      
        })
        
        #  Genome Axis Track
        gtr <- reactive({
          if(!is.null(input$chr)){
            GenomeAxisTrack(chromosome=input$chr,add53 = TRUE, add35 = TRUE, 
              littleTicks = FALSE, showId=TRUE)
          }
        })
        
        #  Ideogram Track
        itr <- reactive({
          if(!is.null(input$chr)){
            IdeogramTrack(genome=input$ucscgen, chromosome=input$chr,
              showId=TRUE, showBandId=TRUE)
          }
        })
        
        #  Render the track plots.
        output$plotname <- renderPlot({
          if(length(input$window)==0){
            return(NULL)
          }
          else{
            itr <- itr()
            gtr <- gtr()
            atr <- atr()
            pt <- plotTracks(list(itr, gtr, atr), from=input$window[1],
              to=input$window[2])
            return(pt)
          }
        })
        
        #  Sets max position for the view window slider for the current
        #  chromosome.
        max_end <- reactive({
          if(length(input$chr)!=1){
            return(NULL)
          }
          else{
            return(max(end(object[seqnames(object)==input$chr])))
          }
        })
        
        #  Sets min position for the view window slider for the current
        #  chromosome.
        min_start <- reactive({
          if(length(input$chr)!=1){
            return(NULL)
          }
          else{
            return(min(start(object[seqnames(object)==input$chr])))
          }
        })
        
        #  Maximum GRange width for filter slider.
        max_width <- reactive({
          if(length(input$chr)!=1){
            return(NULL)
          }
          else{
            return(as.numeric(max(ranges(object[seqnames(
              object)==input$chr])@width)))
          }
        })
        
        #  Minimum GRange width for filter slider.
        min_width <- reactive({
          if(length(input$chr)!=1){
            return(NULL)
          }
          else{
            return(as.numeric(min(ranges(object[seqnames(
              object)==input$chr])@width)))
          }
        })
        
        #  Render the plot range window slider.     
        output$window <- renderUI({
          if(is.null(max_end())){
            return()
          }
          else{
            max_end <- max_end()
            min_start <- min_start()
            sliderInput(inputId = "window",
                        label = "Plot Window:",
                        min = min_start, max = max_end,
                        value = c(min_start,max_end), step = 1
            )
          }
        })
        
        #  Render the width filter slider.
        output$width <- renderUI({
          max_width <- max_width()
          min_width <- min_width()
          
          if(length(max_width)==0 || length(min_width)==0){
            return(NULL)
          }
          
          else{
            sliderInput(inputId = "width",
                        label = "Range Length Filter:",
                        min = min_width, max = max_width,
                        value = c(min_width,max_width), step = 1
            )
          }
        })

        #  Render the choose chromosome dropdown.
        cl <- levels(seqnames(object))
        output$choose_chrom <- renderUI({
          chromChoices <- cl
          names(chromChoices) <- cl
          selectInput("chr", "Chromosome", chromChoices)
        })
        
        #  Render the UCSC dropdown
        output$choose_gen <- .choose_gen()

        #  Render the text under the UCSC dropdown        
        output$gen_text <- renderText({
          ucsc_df <- ucscGenomes()
          ucsc_vec <- as.character(ucsc_df$db)
          i <- which(ucsc_vec==input$ucscgen)
          ucsc_text <- paste(as.character(ucscGenomes()[i,2:4]),
            collapse="&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;")
          if(nchar(ucsc_text)==0){
            return()
          }
          ucsc_text
        })
        
        #  Dynamically build tabs for checkbox groups for metadata subsetting
        args <- reactive({
          args <- list()
          for(i in mcolnames){
            tx <- as.data.frame(object)
            tx <- sort(unique(tx[,i]))
            args[[i]] <- tabPanel(i, checkboxGroupInput(i, 
              paste("Select ", i,sep=""),tx,selected=tx))
          }
          args
        })
        output$mcoltabset <- renderUI({
          args <- args()
          do.call(tabsetPanel, args)
        })

        #  Save Button  
        observe({
          if (input$savebutton == 0)
            return()
          isolate({
            temp <- list()
            for(i in names(bank)){
              p <- unlist(bank[[i]])
              sgr <- subgr(object,i,p[1],p[2],p[3],p[4],p[5],mcolnames,input)
              temp[[i]] <- sgr
            }
            subgrl <- GRangesList(temp)
            subgr <- do.call(c, unname(as.list(subgrl)))
            stopApp(returnValue=subgr)
          })
        })
        
        #  Deposit Button
        observe({
          if (input$bankbutton == 0)
            return()
          isolate({
            bank[[input$chr]] <<- c(input$strand,input$window[1],
              input$window[2],input$width[1],input$width[2])
            output$btable <- renderTable({
              df <- t(as.data.frame(bank))
              colnames(df) <- c("Strand","Min Position","Max Position",
                "Min Width","Max Width")
              df
            })
          })
        })
        
        #  Clear Button 
        observe({
          if (input$clearbutton == 0)
            return()
          isolate({
            bank <<- list()
            output$btable <- renderTable({
              return()
            })
          })
        })
        
      }
    )
    #myRunApp(app, ...)
    runApp(app)
  })

