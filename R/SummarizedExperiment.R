################################################################################
###   SummarizedExperiment
################################################################################

setMethod("display", 
          signature(object = "SummarizedExperiment"), 
          function(object, ...){

            app <- list(
              ui =
                bootstrapPage(
                .jstags(),
                
                sidebarPanel(
                  h3("Summarized Experiment", align = "center"),
                  HTML("<hr />"),
                  uiOutput("choose_chrom"),
                  HTML("<hr />"),
                  #dummy slider until shiny bug gets fixed
                  conditionalPanel
                  (
                    condition = '0==1',
                    sliderInput("dummyslider", "", min=0, max=1, value=0)
                  ),
                  uiOutput("binsui"),
                  HTML("<hr />"),
                  uiOutput("window"),
                  HTML("<hr />"),
                  HTML("Use the mouse to drag and pan the plot. Use the 
                        mousewheel to zoom in/out."),
                  HTML("<hr />"),
                  actionButton("closebutton", "Stop Application")
                ),
                
                mainPanel(
                  .csstags(), 
                  shiny::tags$head(
                    shiny::tags$style(type='text/css', "

                    heat {
                      height: 800px;
                    }
                
                    net {
                      height: 800px;
                    }
                
                    svg {
                      height: 800px;
                    }
                
                    ")),
                  .loading_gif(),
                  uiOutput("heat")
                )
              ),
              
              server = function(input, output){
                
                #  Sets max position for the view window slider for the current
                #  chromosome.
                max_end <- reactive({
                  if(length(input$chr)==1){
                    return(max(end(object[seqnames(object)==input$chr])))
                  }
                  else{
                    return(NULL)
                  }
                })
                
                #  Sets min position for the view window slider for the current 
                #  chromosome.
                min_start <- reactive({
                  if(length(input$chr)==1){
                    min(start(object[seqnames(object)==input$chr]))
                  }
                  else{
                    return(NULL)
                  }
                })
                
                #  Render the plot range window slider.     
                output$window <- renderUI({
                  max_end <- max_end()
                  min_start <- min_start()
                  if(is.numeric(max_end) && is.numeric(min_start)){
                    sliderInput(inputId = "window",
                                label = "Chromosome Range:",
                                min = min_start,
                                max = max_end,
                                value = c(min_start,max_end),
                                step = 1)
                  }
                  else{
                  return(NULL)                  
                  }
                })
                
                #  Render the bins slider.     
                output$binsui <- renderUI({
                  sliderInput(inputId = "bins",
                              label = "Number of Bins",
                              min = 10, max = 100, value = 30, step = 1)
                })
                        
                smaller <- reactive({
                  if(length(input$chr)!=0 &&
                     length(input$bins)!=0 &&
                     length(input$window)!=0){
                    
                    si <- which(as.character(seqnames(rowData(
                      object)))==input$chr)
                    subr <- rowData(object)[si]
                    
                    subr <- subr[start(subr) > input$window[1]]
                    subr <- subr[end(subr) < input$window[2]]
                    
                    orn <- subr$id[order(start(subr))]
                    rfh <- assays(object)[[1]][orn,]
                    ng <- dim(rfh)[1]
                    gs <- split(1:ng,round(as.numeric(
                          cut(1:ng,as.numeric(input$bins)))))
                    
                    smaller <- c()
                    for(i in 1:length(gs)){
                      new <- apply(rfh[gs[[i]],],2,mean)
                      smaller <- rbind(smaller,new)
                    }
                    rownames(smaller) <- 1:length(gs)
                    
                    return(smaller)
                  }          
                })
                
                #  Render the choose chromosome dropdown.
                cl <- as.character(unique(seqnames(rowData(object))))
                output$choose_chrom <- renderUI({
                  chromChoices <- cl
                  names(chromChoices) <- cl
                  selectInput("chr", "Chromosome", chromChoices)
                })
                              
                #  Close Button  
                observe({
                  if (input$closebutton == 0)
                    return()
                  isolate({
                    stopApp()
                  })
                })
                  
                #  Heat Plot
                output$heat <- renderUI({
                  smaller <- smaller()

                  if(is.null(smaller)){
                    return(NULL)
                  }
                  else{  
                    smaller <- smaller()
                        
                    melted <- melt(smaller)
                    names(melted) <- c("Var1","Var2","value")
                    
                    myPalette <- colorRampPalette(rev(brewer.pal(11, 
                    "Spectral")))
                    
                    gp <- ggplot(melted, aes(x = Var1, y = Var2, fill = value))
                    gp <- gp + geom_tile()
                    gp <- gp + coord_fixed()
                    gp <- gp + scale_fill_gradientn(colours = myPalette(100))
                    gp <- gp + scale_x_discrete(expand = c(0, 0))
                    gp <- gp + scale_y_discrete(expand = c(0, 0))
                    gp <- gp + coord_equal()
                    gp <- gp + theme_bw()
                    gp <- gp + xlab(paste(
                      "Region:  ",
                      input$window[1],
                      "  -  ",
                      input$window[2],
                      sep=""))
                    gp <- gp + ylab("Samples")
                    gp <- gp + 
                          ggtitle("Binned Mean Counts by Position") + 
                          theme(plot.title = element_text(lineheight=.8,
                          face="bold", vjust = 2))
                    gp <- gp + theme(axis.text.x = element_blank(),
                                     axis.ticks = element_blank(),
                                     axis.title.y = element_text(vjust = -1),
                                     panel.background = element_blank(),
                                     panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(),
                                     panel.border = element_blank())

                    svgjs <- grid2jssvg(gp)
                    return(svgjs)
                  }
                })
                 
              }
            )
            runApp(app, ...)
          })
