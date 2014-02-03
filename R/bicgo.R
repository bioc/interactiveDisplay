
## The loading gif/panel
.loading_gif <- function(){
  list(
    conditionalPanel(condition="$('html').hasClass('shiny-busy')",
      img(src = 
      "http://upload.wikimedia.org/wikipedia/commons/d/de/Ajax-loader.gif")),
    conditionalPanel(condition="!($('html').hasClass('shiny-busy'))", br())
  )
}

setGeneric("bicgo", function(object, ...)
  standardGeneric("bicgo")
)

setMethod("bicgo", 
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
                    selectInput("bicmethod", "Choose Bicluster Method:",
                                       c("BCPlaid")),
                                         #"BCCC",
                                         #"BCXmotifs",
                                         #"BCSpectral",
                                         #"BCBimax",
                                         #"BCQuest")),
                    HTML("<hr />"),
                    uiOutput("clustern"),
                    uiOutput("cutoff"),
                    HTML("<hr />"),
                    numericInput(inputId="setsize",
                                label= "Min probe pop for GO term",
                                min=1,max=1000,value=10,step=1),
                    HTML("<hr />"),
                    actionButton("savebutton", "Save to Console"),
                    HTML("<hr />"),
                    tableOutput("probeT"),
                    tableOutput("sampleT")
                    
                  ),
                  mainPanel(
                    .loading_gif(),
                    plotOutput("heat"),
                    plotOutput("zoomed"),
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
                  ex <- exprs(object)
                  tmpdata <- as.matrix(as.data.frame(ex))
                  p <- rev(order(apply(tmpdata,1,var)))
                  s <- rev(order(apply(tmpdata,2,var)))
                  p <- p[seq_len(input$cutoff)]
                  s <- s[seq_len(26)]
                  ex <- ex[p,s]     
                  ex <- ex[hclust(dist(ex))$order,hclust(dist(t(ex)))$order]
                  return(ex)
                })
                
                # Bicluster expression data
                bic <- reactive({
                  ex <- ex()
                  bn <- 0
                  while(bn == 0){
                  bic <- biclust(as.matrix(ex),
                          method=input$bicmethod,
                          back.fit = 2,
                          shuffle = 3,
                          fit.model = ~m + a + b,
                          iter.startup = 5,
                          iter.layer = 30,
                          verbose = TRUE)
                  bn <- bic@Number
                  }
                  bic
                })
                
                # Halfway to results object
                out <- reactive({
                  ex <- ex()
                  bic <- bic()
                  out <- list()
                  for(i in seq(bic@Number)){
                    out[[i]] <- list(rownames(ex)[bic@RowxNumber[,i]],
                                     colnames(ex)[bic@NumberxCol[i,]])
                  }
                  out
                })
                
                # Render bicluster dropdown
                output$clustern <- renderUI({
                  if(length(bic())!=1){
                    return(NULL)
                  }
                  else{
                    bic <- bic()
                    return(selectInput("clusterN", "Choose Cluster:",
                              seq(bic@Number)))
                  }
                })
                
                # Results to be displayed and saved back to console
                results <- reactive({
                  bic <- bic()
                  out <- out()
                  for (i in seq(bic@Number)){
                    sets <- Filter(function(x) length(x) >= input$setsize,
                                   split(map$PROBEID, map$TERM))
                    universe <- unlist(sets, use.names=FALSE)
                    siggenes <- out[[i]][[1]]
                    sigsets <- Map(function(x, y) x[x %in% y], sets,
                                   MoreArgs=list(y=siggenes))
                    result <- as.data.frame(hyperg(sets, sigsets, universe))
                    result <- result[rev(order(as.numeric(result[,6]))),]
                    result <- cbind(rownames(result),result)
                    out[[i]][[3]] <- result
                  }
                  return(out)
                })
                
                # Current highlighted probes
                output$probeT <- renderTable({
                  if(length(input$clusterN)!=1){
                    return(NULL)
                  }
                  else{
                    results <- as.data.frame(results()[[as.numeric(as.character(
                      input$clusterN))]][[1]])
                    names(results) <- c("Probes")
                    results
                    return(results)
                  }
                })
                
                # Current highlighted samples
                output$sampleT <- renderTable({
                  if(length(input$clusterN)!=1){
                    return(NULL)
                  }
                  else{
                    results <- as.data.frame(results()[[as.numeric(as.character(
                      input$clusterN))]][[2]])
                    names(results) <- c("Samples")
                    return(results)
                  }
                })
                
                # GO data
                output$GOtable <- renderDataTable({
                  if(length(input$clusterN)!=1){
                    return(NULL)
                  }
                  else{
                    results <- as.data.frame(results()[[as.numeric(as.character(
                      input$clusterN))]][[3]])
                    return(results)
                  }
                })
                
                #  Simple heatplot with highlighting
                output$heat <- renderPlot({
                  if(length(input$clusterN)!=1){
                    return(NULL)
                  }
                  else{  
                    ex <- ex()
                    bic <- bic()
                    color.map <- function(x,i){
                      if(x %in% results()[[as.numeric(as.character(
                        input$clusterN))]][[i]]){
                        return("#FF0000")
                      }
                      else{
                        return("#0000FF")
                      }
                    }
                    colcolors <- unlist(lapply(colnames(ex), color.map,2))
                    rowcolors <- unlist(lapply(rownames(ex), color.map,1))
                    cim(ex,trace="none",ColSideColors=colcolors,RowSideColors=rowcolors)
                    #image(ex,ColSideColors=colcolors,RowSideColors=rowcolors)
                    #heatmap.2(ex,trace="none",ColSideColors=colcolors,RowSideColors=rowcolors)
                    #heatmapBC(ex,bic,trace="none",ColSideColors=colcolors,RowSideColors=rowcolors)
                    #heatmap(ex,ColSideColors=colcolors,RowSideColors=rowcolors)
                    #drawHeatmap2(ex,bicResult=bic,number=NA,plotAll=FALSE)#,
                                 #ColSideColors=colcolors,
                                 #RowSideColors=rowcolors)
                    
                  }
                })

                output$zoomed <- renderPlot({
                  if(length(input$clusterN)!=1){
                    return(NULL)
                  }
                  else{  
                    ex <- ex()
                    p <- results()[[as.numeric(as.character(input$clusterN))]][[1]]
                    s <- results()[[as.numeric(as.character(input$clusterN))]][[2]]
                    cim(ex[p,s],trace="none")
                  }
                })
                
                maxprobes <- dim(exprs(object))[1]
                
                #  Subset probes by variance of expression
                output$cutoff <- renderUI({
                  numericInput(inputId="cutoff",
                              label= "Number of top probes by variance",
                              min=1,max=maxprobes,value=100,step=1)
                })
                
                #  Save Button  
                observe({
                  if (input$savebutton == 0)
                    return()
                  isolate({
                    results <- results()
                    stopApp(returnValue=results)
                  })
                })
                
              }
)
runApp(app, ...)
})
