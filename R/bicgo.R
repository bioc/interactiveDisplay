
## The loading gif/panel
#.loading_gif <- function(){
#  list(
#    conditionalPanel(condition="$('html').hasClass('shiny-busy')", img(src = 
#      "http://upload.wikimedia.org/wikipedia/commons/d/de/Ajax-loader.gif")),
#    conditionalPanel(condition="!($('html').hasClass('shiny-busy'))", br())
#  )
#}

.bicgo <- function(object, ...){
            
  .usePackage('shiny')
  .usePackage('GOstats')
  #.usePackage('biclust')
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
          actionButton("nextbutton", "Next Pass"),
          actionButton("resetbutton", "Reset"),
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
          HTML("Each GO Summary submission will be bundled and returned to the console"),
          HTML("<hr />"),
          actionButton("stopbutton", "Stop/Save")      
        ),
        mainPanel(
          
          tags$head(
            tags$style(type='text/css', 
                       "#geneT { background-color: #DCE8BE; } 
                        #GOtable { background-color: #CABEE8 }")),                    
          .loading_gif(),
          
          tabsetPanel(
            tabPanel("Plots",
                     HTML("Overview"),
                     plotOutput("heat"),
                     HTML("Selected Probes/Samples"),
                     plotOutput("zoomed")
            ),
            tabPanel("Selected Probes/Genes",
                     dataTableOutput("geneT")
            ),
            tabPanel("Selected Samples",
                     tableOutput("sampleT")
            )
          ),
          tabsetPanel(
            #HTML("<hr />"),
            #HTML("GO Summary of Selected Probes"),
            tabPanel("GO Summary of Selected Probes", dataTableOutput("GOtable"))
          )
        )
      ),
    
    server = function(input,output) {
      
      
      ################################################################################
      # Some nonreactive processing to do just once.
      
      #pkgName <- paste(annotation(object),".db",sep="")
      pkgName <- "hgu95av2.db"
      pkg <- get(pkgName)
      
      res <- suppressWarnings(
        select(pkg, keys(hgu95av2.db),
                              c("ENTREZID","GENENAME","GO"), "PROBEID"))
      resa <- cbind(res$PROBEID,res$GO)
      resb <- resa[!duplicated(resa),]
      resc <- resb[!is.na(resb[,2]),]
      
      
      map <- (suppressWarnings(
        select(GO.db, resc[,2], "TERM")))
      map <- cbind(resc[,1],map)
      names(map) <- c("PROBEID","GOID","TERM")
      
      ################################################################################
      record <- list()
      lex <- c()
      r <- 0
      ################################################################################             
  
      observe({
        if(input$resetbutton == 0){
          return()
        }
        else{
          r <<- 0
        }
      })
      
      observe({
        if(input$nextbutton == 0){
          return()
        }
        else{
          r <<- 1
        }
      })
      
      # Subset expression data
      tex <- reactive({
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
      
      ex <- reactive({
        input$resetbutton
        input$nextbutton
        if (r==0){
          lex <<- tex()
          return(tex())
        }
        else{
          isolate({
            if(length(ex)==0){
              return(NULL)
            }
            else{
              if(r==0 || length(input$clusterncol) == 0 || length(input$clusternrow) == 0){
                lex <<- tex()
                return(tex())
              }
              else{
                colcut <- cutree(hclust(dist(t(lex))),input$clusterncol)
                rowcut <- cutree(hclust(dist(lex)),input$clusternrow)
                
                temp <- lex[(rowcut %in% as.numeric(input$selclustrow)),(colcut %in% as.numeric(input$selclustcol))]
                if(sum(dim(lex))>3){
                  lex <<- temp
                }
                return(lex)
              }
            }
          })
        }
      })
  
  
      output$clusterncol <- renderUI({
        ex <- ex()
        maxsamples <- dim(ex)[2]
        if(length(maxsamples)!=1){
          return(NULL)
        }
        else{
        numericInput(inputId="clusterncol",
                     label= "Number of sample clusters",
                     min=1,max=maxsamples,value=2,step=1)
        }
      })
  
      output$clusternrow <- renderUI({
        ex <- ex()
        maxprobes <- dim(ex)[1]
        if(length(maxprobes)!=1){
          return(NULL)
        }
        else{
        numericInput(inputId="clusternrow",
                     label= "Number of probe clusters",
                     min=1,max=maxprobes,value=2,step=1)
        }
      })
      
      colcut <- reactive({
        ex <- ex()
        if(length(ex)==0 || length(input$clusterncol) == 0){
          return(NULL)
        }
        else{
          return(cutree(hclust(dist(t(ex))),input$clusterncol))
        }
      })
      
      rowcut <- reactive({
        ex <- ex()
        if(length(ex)==0 || length(input$clusternrow) == 0){
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
      
      hgenes <- reactive({
        ex <- ex()
        if(length(ex)==0){
          return(NULL)
        }
        else{
          hprobes <- hprobes()
          return(as.data.frame(select(hgu95av2.db, hprobes, c("ENTREZID","GENENAME"), "PROBEID")))
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
      
      # Current highlighted genes
      output$geneT <- renderDataTable({
        hgenes <- hgenes()
        if(length(hprobes)==0){
          return(NULL)
        }
        else{
          results <- as.data.frame(hgenes)
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
            
            hgenes <- hgenes()
            hsamples <- hsamples()
            
            
            record[[input$gobutton]] <<- list(hprobes,
                                     hsamples,
                                     hgenes,
                                     result)
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
        if(length(hprobes)==0 || length(input$clusterncol) == 0 || length(input$clusternrow) == 0){
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
        if(length(ex)==0 || length(input$clusterncol) == 0 || length(input$clusternrow) == 0){
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
        if(length(maxprobes)!=1){
          return(NULL)
        }
        else{
        numericInput(inputId="cutoff",
                     label= "Number of top probes by variance",
                     min=2,max=maxprobes,value=500,step=1)
        }
      })
      
      #  Stop Button  
      observe({
        if (input$stopbutton == 0)
          return()
        isolate({
          stopApp(returnValue=record)
        })
      })
      
    }
  )
  runApp(app, ...)
}

setGeneric("bicgo", function(object, ...)
  standardGeneric("bicgo")
)

setMethod("bicgo", 
          signature(object = c("ANY")),
          function(object, ...){
            .bicgo(object=object, ...)
          })
