################################################################################
###   GRangesList
################################################################################

require(GenomicRanges)
require(Gviz)
require(rtracklayer)

## helper for setting up sidebar
.GL_setSidebarPanel <- function(){
  sidebarPanel(
    uiOutput("choose_grange"),
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
  signature(object = "GRangesList"), 
  function(object,  sflag = TRUE, ...){
        
    app <- list(
      ui = bootstrapPage(        
        .jstags(),  
        .csstags(),
        h3("Genomic Ranges List"),
        .loading_gif(),
        .GL_setSidebarPanel(),
        .GR_GRL_setMainPanel(sflag)
      ),
      
      server = function(input, output, session){
        
        # This stores parameters for subsetted GRanges per chromosome as a list.
        bank <- list()
        
        t_object <- reactive({
          object[[input$gr]] 
        })
        
        mcolnames <- reactive({
          t_object <- t_object()
          as.character(names(mcols(t_object)))
        })
        
        #  This GRanges object is subsetted by user input in the shiny widget.
        s_object <- reactive({
          t_object <- t_object()
          subgr2(t_object,input$chr,input$strand,input$width,input$window,
            mcolnames(),input)
        })
        
        #  The full submitted GRanges object converted to a data frame for the 
        #  perpose of rendering a table in shiny. 
        output$fulltable <- renderDataTable({
          t_object <- t_object()
          as.data.frame(t_object)
        })
        
        #  The subsetted GRanges object converted to a data frame for the 
        #  perpose of rendering a table in shiny.
        output$rtable <- renderDataTable({
          s_object <- s_object()
          as.data.frame(s_object)
        })
        
        #  Annotation Track     
        atr <- reactive({
          s_object <- s_object()
          if(!is.null(input$chr)){
            AnnotationTrack(s_object, chromosome=input$chr,
            name="Genomic Ranges Annotation Track",
            fill="black", background.panel = "#f5f5f5",
            background.title = "#2aa5e5", cex.title = 1.1)
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
        
        
        # The circle plot
        
        cplot <- reactive({
          
          p <- ggplot() + layout_circle(object,
                                        geom = "ideo",
                                        fill = "gray70",
                                        radius = 30,
                                        trackWidth = 4)
          p <- p + layout_circle(object,
                                 geom = "scale",
                                 size = 2, radius = 35,
                                 trackWidth = 2)
          #p <- p + layout_circle(object,
          #                       geom = "text", 
          #                       aes(label = seqnames), 
          #                       vjust = 0, 
          #                       radius = 38, 
          #                       trackWidth = 7)
          p <- p + layout_circle(object, 
                                 geom = "rect", 
                                 color = "steelblue", 
                                 radius = 23 , 
                                 trackWidth = 6)
          
          gr <- object
          gr <- sample(gr,input$linkn)
          values(gr)$to.gr <- sample(gr,input$linkn)
          
          
          values(gr)$rearrangements <- ifelse(as.character(seqnames(gr)) == as.character(seqnames((values(gr)$to.gr))), "intrachromosomal", "interchromosomal")
          
          
          p <- p + layout_circle(gr, 
                                 geom = "link", 
                                 linked.to = "to.gr", 
                                 aes(color = rearrangements), 
                                 radius = 22, 
                                 trackWidth = 1)
          p
        })
        
        if(sflag==TRUE){
          output$cplot <- renderUI({
            #progress <- Progress$new(session, min=1, max=10)
            #on.exit(progress$close())
            
            #progress$set(message = 'Calculation in progress',
            #             detail = 'This may take a while...')
            
            #for (i in 1:10) {
            #  progress$set(value = i)
            #  Sys.sleep(0.1)
            #}
            grid2jssvg(cplot())
          })
        }
        else{
          output$cplot <- renderPlot({
            #progress <- Progress$new(session, min=1, max=10)
            #on.exit(progress$close())
            
            #progress$set(message = 'Calculation in progress',
            #             detail = 'This may take a while...')
            
            #for (i in 1:10) {
            #  progress$set(value = i)
            #  Sys.sleep(0.1)
            #}
            cplot()
          })        
        }
        
        
        
        #  Sets max position for the view window slider for the current
        #  chromosome.
        max_end <- reactive({
          t_object <- t_object()
          max(end(t_object[seqnames(t_object)==input$chr]))
        })
        
        #  Sets min position for the view window slider for the current
        #  chromosome.
        min_start <- reactive({
          t_object <- t_object()
          min(start(t_object[seqnames(t_object)==input$chr]))
        })
        
        #  Maximum GRange width for filter slider.
        max_width <- reactive({
          t_object <- t_object()
          as.numeric(max(ranges(t_object[seqnames(t_object)==input$chr])@width))
        })
        
        #  Minimum GRange width for filter slider.
        min_width <- reactive({
          t_object <- t_object()
          as.numeric(min(ranges(t_object[seqnames(t_object)==input$chr])@width))
        })
        
        #  Render the plot range window slider.     
        output$window <- renderUI({
          max_end <- max_end()
          min_start <- min_start()
          sliderInput(inputId = "window",
                      label = "Plot Window:",
                      min = min_start,
                      max = max_end,
                      value = c(min_start,max_end),
                      step = 1)
        })
        
        #  Render the width filter slider.
        output$width <- renderUI({
          max_width <- max_width()
          min_width <- min_width()
          sliderInput(inputId = "width",
                      label = "Range Length Filter:",
                      min = min_width,
                      max = max_width,
                      value = c(min_width,max_width),
                      step = 1)
        })
        
        gl <- names(object)
        output$choose_grange <- renderUI({
          grangeChoices <- gl
          names(grangeChoices) <- gl
          selectInput("gr", "GRange", grangeChoices)
        })
        
        #  Render the choose chromosome dropdown.
        #cl <- levels(seqnames(object()))
        output$choose_chrom <- renderUI({
          t_object <- t_object()
          cl <- unique(as.character(seqnames(t_object)))
          chromChoices <- cl
          names(chromChoices) <- cl
          selectInput("chr", "Chromosome", chromChoices)
        })
        

        
        #  Render the text under the UCSC dropdown        
        output$choose_gen <- .choose_gen()
        
        #  Dynamically build tabs for checkbox groups for metadata subsetting
        args <- reactive({
          args <- list()
          mcolnames <- mcolnames()
          for(i in mcolnames){
            t_object <- t_object()
            tx <- as.data.frame(t_object)
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
            tempy <- list()
            for(i in names(bank)){
              temp <- list()
              tbank <- bank[[i]]
              for(j in names(tbank)){
                p <- unlist(tbank[[j]])
                sgr <- subgr(object[[i]],p[2],p[3],p[4],p[5],p[6],p[7],
                  mcolnames(),input)
                temp[[j]] <- sgr
              }     
              subgrly <- GRangesList(temp)
              tempy[[i]] <- do.call(c, unname(as.list(subgrly)))
            }
            subgrl <- GRangesList(tempy)
            stopApp(returnValue=subgrl)
          })
        })
        
        #  Deposit Button
        observe({
          if (input$bankbutton == 0)
            return()
          isolate({
            bank[[input$gr]][[input$chr]] <<- c(input$gr,input$chr,input$strand,
              input$window[1],input$window[2],input$width[1],input$width[2])
            output$btable <- renderDataTable({
              df <- ldply(lapply(bank, ldply))[,-1]
              colnames(df) <- c("GRange","Chromosome","Strand","Min Position",
                "Max Position","Min Width","Max Width")
              df <- as.data.frame(df)
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
            output$btable <- renderDataTable({
              return()
            })
          })
        })
        
      }
    )
    #myRunApp(app, ...)
    runApp(app, ...)
  })
