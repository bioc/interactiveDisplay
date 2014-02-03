
## The loading gif/panel
.loading_gif <- function(){
  list(
    conditionalPanel(condition="$('html').hasClass('shiny-busy')", img(src = 
      "http://upload.wikimedia.org/wikipedia/commons/d/de/Ajax-loader.gif")),
    conditionalPanel(condition="!($('html').hasClass('shiny-busy'))", br())
  )
}

setGeneric("bicgo2", function(object, ...)
  standardGeneric("bicgo2")
)

setMethod("bicgo2", 
          signature(object = c("ANY")),
          function(object, ...){
            
            .usePackage('shiny')
            .usePackage('GOstats')
            .usePackage('biclust')
            .usePackage('hgu95av2.db')
            .usePackage('GO.db')
            .usePackage('gplots')
            .usePackage('mixOmics')
            
            app <- list(
              ui =
                bootstrapPage(
                  
                  sidebarPanel(
                    HTML("Hit \"View/Update GO Summary\" after making selections."),
                    HTML("<hr />"),
                    HTML("Probes/Samples in overview with a red sidebar are selected for the GO summary"),
                    HTML("<hr />"),
                    uiOutput("clusterncol"),
                    uiOutput("clusternrow"),
                    uiOutput("selclustcol"),
                    uiOutput("selclustrow"),
                    HTML("<hr />"),
                    uiOutput("cutoff"),
                    HTML("<hr />"),
                    numericInput(inputId="setsize",
                                 label= "Min probe pop for GO term",
                                 min=1,max=1000,value=10,step=1),
                    HTML("<hr />"),
                    #submitButton("Update View"),
                    actionButton("gobutton", "View/Update GO Summary"),
                    HTML("<hr />"),
                    actionButton("stopbutton", "Stop App"),
                    HTML("<hr />"),
                    HTML("List of Selected Probes/Samples"),
                    tableOutput("probeT"),
                    tableOutput("sampleT")
                    
                  ),
                  mainPanel(
                    .loading_gif(),
                    HTML("Overview"),
                    plotOutput("heat"),
                    HTML("Selected Probes/Samples"),
                    plotOutput("zoomed"),
                    HTML("GO summary of Selected Probes"),
                    dataTableOutput("GOtable")
                  )
                ),
              
              server = function(input,output) {
                
                
                ################################################################################
                # Some nonreactive processing to do just once.
                
                pkgName <- paste(annotation(object),".db",sep="")
                pkg <- get(pkgName)
                
                res <- suppressWarnings(
                  AnnotationDbi::select(pkg, keys(hgu95av2.db),
                                        c("ENTREZID","GENENAME","GO"), "PROBEID"))
                resa <- cbind(res$PROBEID,res$GO)
                resb <- resa[!duplicated(resa),]
                resc <- resb[!is.na(resb[,2]),]
                
                
                map <- (suppressWarnings(
                  AnnotationDbi::select(GO.db, resc[,2], "TERM")))
                map <- cbind(resc[,1],map)
                names(map) <- c("PROBEID","GOID","TERM")
                
                ################################################################################
                
                # Subset expression data
                ex <- reactive({
                  if(length(input$cutoff)!=1){
                    return(NULL)
                  }
                  else{
                    ex <- exprs(object)
                    tmpdata <- as.matrix(as.data.frame(ex))
                    p <- rev(order(apply(tmpdata,1,var)))
                    p <- p[seq_len(input$cutoff)]
                    ex <- ex[p,]     
                    return(ex)
                  }
                })
         
                output$clusterncol <- renderUI({
                  maxsamples <- dim(exprs(object))[1]
                  numericInput(inputId="clusterncol",
                               label= "Number of sample clusters",
                               min=1,max=maxsamples,value=2,step=1)
                })
            
                output$clusternrow <- renderUI({
                  maxprobes <- dim(exprs(object))[1]
                  numericInput(inputId="clusternrow",
                               label= "Number of probe clusters",
                               min=1,max=maxprobes,value=2,step=1)
                })
                
                colcut <- reactive({
                  ex <- ex()
                  if(length(ex)==0){
                    return(NULL)
                  }
                  else{
                    return(cutree(hclust(dist(t(ex))),input$clusterncol))
                  }
                })
                
                rowcut <- reactive({
                  ex <- ex()
                  if(length(ex)==0){
                    return(NULL)
                  }
                  else{
                    return(cutree(hclust(dist(ex)),input$clusternrow))
                  }
                })
            
                output$selclustcol <- renderUI({
                  tx <- seq(input$clusterncol)
                  if(length(tx) == 0){
                    tx <- 1:2
                  } 
                  checkboxGroupInput("selclustcol",
                                     "Select Sample Clusters",
                                     tx,
                                     selected=tx)
                })
                output$selclustrow <- renderUI({
                  tx <- seq(input$clusternrow)
                  if(length(tx) == 0){
                    tx <- 1:2
                  } 
                  checkboxGroupInput("selclustrow",
                                     "Select Probe Clusters",
                                     tx,
                                     selected=tx)
                })
                
                hprobes <- reactive({
                  ex <- ex()
                  if(length(ex)==0){
                    return(NULL)
                  }
                  else{
                    rowcut <- rowcut()
                    return((rownames(ex))[(rowcut %in% as.numeric(input$selclustrow))])
                  }
                })
                
                hsamples <- reactive({
                  ex <- ex()
                  if(length(ex)==0){
                    return(NULL)
                  }
                  else{
                    colcut <- colcut()
                    return((colnames(ex))[(colcut %in% as.numeric(input$selclustcol))])
                  }
                })
                
                # Current highlighted probes
                output$probeT <- renderTable({
                  hprobes <- hprobes()
                  if(length(hprobes)==0){
                    return(NULL)
                  }
                  else{
                    results <- as.data.frame(hprobes)
                    names(results) <- c("Probes")
                    return(results)
                  }
                })
                
                # Current highlighted samples
                output$sampleT <- renderTable({
                  hsamples <- hsamples()
                  if(length(hsamples)==0){
                    return(NULL)
                  }
                  else{
                    results <- as.data.frame(hsamples)
                    names(results) <- c("Samples")
                    return(results)
                  }
                })
                
                #  GO Button
                output$GOtable <- renderDataTable({
                  if (input$gobutton == 0)
                    return()
                  isolate({
                  # GO data
                    hprobes <- hprobes()
                    if(length(hprobes)==0){
                      return(NULL)
                    }
                    else{
                      sets <- Filter(function(x) length(x) >= input$setsize,
                                     split(map$PROBEID, map$TERM))
                      universe <- unlist(sets, use.names=FALSE)
                      siggenes <- hprobes
                      sigsets <- Map(function(x, y) x[x %in% y], sets,
                                     MoreArgs=list(y=siggenes))
                      result <- as.data.frame(hyperg(sets, sigsets, universe))
                      result <- result[rev(order(as.numeric(result[,6]))),]
                      result <- cbind(rownames(result),result)
                      result <- as.data.frame(result)
                      return(result)
                    }
                  })
                })
                
                #  Simple heatplot with highlighting
                output$heat <- renderPlot({
                  hprobes <- hprobes()
                  hsamples <- hsamples()
                  ex <- ex()
                  hprobes <- hprobes()
                  if(length(hprobes)==0 || length(hprobes)==0 || length(hprobes)==0){
                    return(NULL)
                  }
                  else{
                    color.map <- function(x,i){
                      if(x %in% i){
                        return("#FF0000")
                      }
                      else{
                        return("#0000FF")
                      }
                    }
                    colcolors <- unlist(lapply(colnames(ex), color.map, hsamples))
                    rowcolors <- unlist(lapply(rownames(ex), color.map, hprobes))
                    return(cim(ex,trace="none",ColSideColors=colcolors,RowSideColors=rowcolors))
                  }
                })
                
                output$zoomed <- renderPlot({ 
                  ex <- ex()
                  if(length(ex)==0){
                    return(NULL)
                  }
                  else{
                    colcut <- colcut()
                    rowcut <- rowcut()
                    cim(ex[(rowcut %in% as.numeric(input$selclustrow)),(colcut %in% as.numeric(input$selclustcol))],trace="none")
                  }
                })
                                
                #  Subset probes by variance of expression
                output$cutoff <- renderUI({
                  maxprobes <- dim(exprs(object))[1]
                  numericInput(inputId="cutoff",
                               label= "Number of top probes by variance",
                               min=2,max=maxprobes,value=500,step=1)
                })
                
                #  Stop Button  
                observe({
                  if (input$stopbutton == 0)
                    return()
                  isolate({
                    #results <- results()
                    stopApp()
                  })
                })
                
              }
            )
            runApp(app, ...)
          })
