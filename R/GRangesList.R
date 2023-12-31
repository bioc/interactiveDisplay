################################################################################
###   GRangesList
################################################################################

## helper for setting up sidebar
.GL_setSidebarPanel <- function(){
  sidebarPanel(
    h3("Genomic Ranges List", align="center"),
    HTML("<hr />"),
    uiOutput("choose_grange"),
    uiOutput("choose_chrom"),
    HTML("<hr />"),
    uiOutput("choose_gen"),
    uiOutput("gen_text"),
    HTML("<hr />"),
    checkboxInput("suppress","Suppress Ideogram"),
    HTML("<hr />"),
    #dummy slider until shiny bug gets fixed
    conditionalPanel
    (
      condition = '0==1',
      sliderInput("dummyslider", "", min=0, max=1, value=0)
    ),
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
  function(object,  sflag = FALSE, ...){
    
    .usePackage('GenomicRanges')
    .usePackage('Gviz')
    .usePackage('rtracklayer')
    .usePackage('ggbio')
        
    app <- list(
      ui = bootstrapPage(        
        .jstags(),  
        .csstags(),
        shiny::tags$head(
          shiny::tags$style(type='text/css', "

        cplot {
          height: 800px;
        }

        ")),
        .GL_setSidebarPanel(),
        .GR_GRL_setMainPanel(sflag)
      ),
      
      server = function(input, output, session){
        
        # This stores parameters for subsetted GRanges per chromosome as a list.
        bank <- list()
        
        # Allow for more types of GRanges data
        options(ucscChromosomeNames=FALSE)
        
        t_object <- reactive({
          if(!is.null(input$gr)){
            return(object[[input$gr]])
          }
        })
        
        mcolnames <- reactive({
          t_object <- t_object()
          if(!is.null(t_object)){
            mcolnames <- as.character(names(mcols(t_object)))
            #quick fix for S4 metadata
            temp_new <- c()
            for(i in mcolnames){
              if(dim(as.data.frame(mcols(t_object)[i]))[2] == 1){
                temp_new <- c(temp_new,i)
              }
            }
            return(temp_new)
          }
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
          if(!is.null(input$chr) && input$suppress == FALSE){
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
            if(input$suppress == FALSE){
              pt <- plotTracks(list(itr, gtr, atr), from=input$window[1],
                to=input$window[2])
            }
            else{
              pt <- plotTracks(list(gtr, atr), from=input$window[1],
                               to=input$window[2])
            }
            return(pt)
          }
        })
        outputOptions(output, "plotname", suspendWhenHidden = FALSE)
        
        # The circle plot 
        cplot <- reactive({
          if(length(input$gr)==0 || length(t_object)==0){
            return(NULL)
          }
          else{
            t_object <- t_object()
            p <- ggplot() + layout_circle(t_object,
                                          geom = "ideo",
                                          radius = 30,
                                          trackWidth = 4,
                                          aes(fill=seqnames))
            p <- p + layout_circle(t_object,
                                   geom = "scale",
                                   size = 2, radius = 35,
                                   trackWidth = 2)
                        
            p <- p + layout_circle(t_object, 
                                   geom = "rect", 
                                   color = "steelblue4", 
                                   radius = 23 , 
                                   trackWidth = 6)
            
            return(p)
          }
        })

        #Dev option for suppressing svg output
        if(sflag==TRUE){
          output$cplot <- renderUI({
            if(length(cplot())==0){
              return(NULL)
            }
            else{
              cplot <- cplot()
              g <- grid2jssvg(cplot)
              return(g)
            }
          })
        }
        else{
          output$cplot <- renderPlot({
            if(length(cplot())==0){
              return(NULL)
            }
            else{
              cplot <- cplot()
              return(plot(cplot))
            }
          })
        }
        
        #  Sets max position for the view window slider for the current
        #  chromosome.
        max_end <- reactive({
          t_object <- t_object()
          
          if(length(input$chr)!=1){
            return(NULL)
          }
          else{
            return(suppressWarnings(
              max(end(t_object[seqnames(t_object)==input$chr]))))
          }
        })
        
        #  Sets min position for the view window slider for the current
        #  chromosome.
        min_start <- reactive({
          t_object <- t_object()
          if(length(input$chr)!=1){
            return(NULL)
          }
          else{
            return(suppressWarnings(
              min(start(t_object[seqnames(t_object)==input$chr]))))
          }
        })
        
        #  Maximum GRange width for filter slider.
        max_width <- reactive({
          t_object <- t_object()
          if(length(input$chr)!=1){
            return(NULL)
          }
          else{
            return(as.numeric(suppressWarnings(max(ranges(t_object[seqnames(
              t_object)==input$chr])@width))))
          }
        })
        
        #  Minimum GRange width for filter slider.
        min_width <- reactive({
          t_object <- t_object()
          if(length(input$chr)!=1){
            return(NULL)
          }
          else{
            return(as.numeric(suppressWarnings(min(ranges(t_object[seqnames(
              t_object)==input$chr])@width))))
          }
        })
        
        #  Render the plot range window slider.     
        output$window <- renderUI({
          max_end <- max_end()
          min_start <- min_start()
          if(is.null(max_end()) || is.infinite(max_end())){
            return()
          }
          else{
            return(sliderInput(inputId = "window",
                        label = "Plot Window:",
                        min = min_start,
                        max = max_end,
                        value = c(min_start,max_end),
                        step = 1))
          }
        })
        
        #  Render the width filter slider.
        output$width <- renderUI({
          max_width <- max_width()
          min_width <- min_width()
          if(is.null(max_width()) || is.infinite(max_width())){
            return()
          }
          else{
          return(sliderInput(inputId = "width",
                      label = "Range Length Filter:",
                      min = min_width,
                      max = max_width,
                      value = c(min_width,max_width),
                      step = 1))
          }
        })
        
        gl <- names(object)
        output$choose_grange <- renderUI({
          grangeChoices <- gl
          names(grangeChoices) <- gl
          selectInput("gr", "GRange", grangeChoices)
        })
        
        #  Render the choose chromosome dropdown.
        output$choose_chrom <- renderUI({
          t_object <- t_object()
          if(!is.null(t_object)){
            cl <- unique(as.character(seqnames(t_object)))
            chromChoices <- cl
            names(chromChoices) <- cl
            return(selectInput("chr", "Chromosome", chromChoices))
          }
        })
        
        #  Render the text under the UCSC dropdown        
        output$choose_gen <- .choose_gen(...)
        
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
            if(length(bank)!=0){
              templ <- list()
              for(i in names(bank)){
                tbank <- bank[[i]]
                temp <- list()
                for(j in names(tbank)){
                  p <- unlist(tbank[[j]])
                  sgr <- subgr(object[[i]],p[2],p[3],p[4],p[5],p[6],p[7],
                    mcolnames(),input)
                  temp[[j]] <- as.data.frame(sgr)
                }
                mergegr <- do.call(rbind,temp)
                mc <- apply(mergegr[names(mergegr)[-(1:5)]],2,list)
                gr <- GRanges(seqnames=mergegr$seqnames,
                              IRanges(mergegr$start,mergegr$end),
                              mergegr$strand,
                              mc)
                seqlengths(gr) <- seqlengths(object[[i]])[levels(seqnames(gr))]
                templ[[i]] <- gr
              }
              grl <- GRangesList(templ)
              grl <- lapply(grl,unname)
              stopApp(returnValue=grl)
            }
            else{
              stopApp()
            }
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
    interactiveDisplayBase::.runApp(app, ...)
  })
