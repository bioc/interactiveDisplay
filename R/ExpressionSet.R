#####################################################################################################
###   ExpressionSet
#####################################################################################################

heatcolor1 <- function(inputId1) {
  tagList(
    shiny::tags$input(id = inputId1, class = "color", value = "F7FF00", onchange = "window.Shiny.onInputChange('color1', this.color.toString())")
  )
}

heatcolor2 <- function(inputId2) {
  tagList(
    shiny::tags$input(id = inputId2, class = "color", value = "050505", onchange = "window.Shiny.onInputChange('color2', this.color.toString())")
  )
}

heatcolor3 <- function(inputId3) {
  tagList(
    shiny::tags$input(id = inputId3, class = "color", value = "1736FF", onchange = "window.Shiny.onInputChange('color3', this.color.toString())")
  )
}


.ES_setSidebarPanel <- function(){
  sidebarPanel(
    tableOutput("expinfo"),
    HTML("<hr />"),
    selectInput("either", "Network/Dendrogram View:  Sample or Probe:", choices = c("probe","sample")),
    HTML("<hr />"),
    uiOutput("choose_probe"),
    HTML("<hr />"),
    selectInput("order", "Show Top or Bottom Ranked,", choices = c("top","bottom")),
    selectInput("measure", "Based On:", choices = c("variance","average")),
    uiOutput("pmeancutoff"),
    uiOutput("smeancutoff"),
    sliderInput(inputId = "tweak",
                label = "Tweak Axis Label Font Size",
                min = .1, max = 2, value = 1, step = .1),
    HTML("<hr />"),
    #sliderInput(inputId = "con_knum",
    #            label = "Number of Clusters",
    #            min = 1, max = 100, value = 4, step = 1),
    uiOutput("edge"),
    uiOutput("gen_text"),
    HTML("<hr />"),
    numericInput("con_knum", "Number of Clusters:", 5),
    selectInput("hc_method", "Hierarchical Clustering Method", choices = c("ward", "single", "complete", "average", "mcquitty", "median", "centroid")),
    selectInput("dist_method", "Distance/Similarity Method", choices = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")),
    HTML("<hr />"),
    radioButtons('rainbow', 'Color Scale',
                 c('Rainbow'='default',
                   'Three Color'='tri'),
                 'Rainbow'),
    heatcolor1('exampleTextarea1'),
    heatcolor2('exampleTextarea2'),
    heatcolor3('exampleTextarea3'),
    HTML("<hr />"),
    actionButton("closebutton", "Stop Application")
  )
}

.ES_setMainPanel <- function(){
  mainPanel(
    includeHTML(system.file("www", "graph.js", package="interactiveDisplay")),
    .jstags(),
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

    ")
    ),
    tabsetPanel(
      tabPanel("Heat Plot",uiOutput("heat")),
      tabPanel("Network View",uiOutput("svg")),
      tabPanel("Dendrogram",plotOutput("dendro"))
    ),
    tabsetPanel(
      tabPanel("GO Summary",
        tableOutput("gotest1"),
        tableOutput("gotest2")
      )
    )
  )
}

setMethod("display",  
  signature(object = c("ExpressionSet")), 
  function(object, ...){
        
    app <- list(
      ui =
        bootstrapPage(
          #.jstags(),
          h3("Expression Set"),
          .loading_gif(),
          .ES_setSidebarPanel(),
          .ES_setMainPanel()
        ),
      
      server = function(input, output){
        
        #  Subset the submitted data
        tmpdata <- reactive({
          tmpdata <- as.matrix(as.data.frame(exprs(object)))
          if(is.null(input$pmeancutoff)){
            return()
          }
          else{
            
            if(input$measure=="average"){
              p <- order(apply(tmpdata,1,mean))
              s <- order(apply(tmpdata,2,mean))
            }
            if(input$measure=="variance"){
              p <- order(apply(tmpdata,1,var))
              s <- order(apply(tmpdata,2,var))
            }
            
            if(input$order=="top"){
              p <- rev(p)
              s <- rev(s)
            }
            if(input$order=="bottom"){

            }
            
            p <- p[1:input$pmeancutoff]
            s <- s[1:input$smeancutoff]
            
            #p <- rev(order(apply(tmpdata,1,mean)))[1:input$pmeancutoff]
            #s <- rev(order(apply(tmpdata,2,mean)))[1:input$smeancutoff]
            tmpdata <- tmpdata[p,s]
            return(tmpdata)
          }
        })
        
        #  Table of experiment info if present
        output$expinfo <- renderTable({
          expinfo <- as.data.frame(expinfo(experimentData(object)))
          names(expinfo) <- c("Experiment Info")
          expinfo
        })
        
        #  Probe Selection dropdown for GO summary
        output$choose_probe <- renderUI({
          prChoices <- sort(rownames(tmpdata()))
          names(prChoices) <- prChoices
          selectInput("probe", "Chosen Probe for GO Summary", prChoices)
        })
     
        #  GO summary
        output$gotest1 <- renderTable({
          if(length(input$probe)!=1){
            return()
          }
          else{
            pkgName <- paste(annotation(object),".db",sep="")
            try(require(pkgName,character.only=TRUE),silent=TRUE)
            if(exists(pkgName)==FALSE){
              return(as.data.frame("No annotation package available"))
            }
            else{
              pkg <- get(pkgName)
              if(class(pkg)=="ChipDb"){
                res <- select(pkg, input$probe, c("ENTREZID","GENENAME"), "PROBEID")
                return(res)
              }
              else{
                return("Object does not have a ChipDb annotation")
              }
            }
          }
        })
        
        #  GO summary
        output$gotest2 <- renderTable({
          if(length(input$probe)!=1){
            return()
          }
          else{
            pkgName <- paste(annotation(object),".db",sep="")
            try(require(pkgName,character.only=TRUE),silent=TRUE)
            if(exists(pkgName)==FALSE){
              return(as.data.frame("No annotation package available"))
            }
            else{
              pkg <- get(pkgName)
              if(class(pkg)=="ChipDb"){
                res <- select(pkg, input$probe, c("ENTREZID","GENENAME","GO"), "PROBEID")
                res2 <- head(select(GO.db, res$GO, "TERM", "GOID"))
                return(res2)
              }
              else{
                return(as.data.frame("Object does not have a ChipDb annotation"))
              }
            }
          }
        })
        
        #  Data for network view, sample or probe
        data <- reactive({
          data <- tmpdata()
          if(input$either=="sample"){
            data <- data
          }else{
            data <- t(data)
          }
          data
        })
        
        #  Subset probes by average expression
        output$pmeancutoff <- renderUI({
          textInput("pmeancutoff", "Number of Probes to Display", 20)   
          #textInput("pmeancutoff", "Top Probes by Average Expression", dim(exprs(object))[1])
        })
        
        #  Subset samples by average expression
        output$smeancutoff <- renderUI({   
          textInput("smeancutoff", "Number of Samples to Display", dim(exprs(object))[2])
        })
        
        #  This determines the distance threshold needed for the desired number of edges
        cutoff <- reactive({
          data <- data()
          val <- cor(data())
          diag(val) <- 0
          val[lower.tri(val)] <- 0
          cutoff <- sort(val[val>0],decreasing=TRUE)[input$edgenum]
          cutoff
        })
        
        #  Build the network
        output$net <- reactive({
          data <- data()
          hc <- hc()
          val <- cor(data())
          if (is.null(val)){
            return(list(names=character(), links=list(source=-1, target=-1)))
          }
          diag(val) <- 0
          val[lower.tri(val)] <- 0
          cutoff <- sort(val[val>0],decreasing=TRUE)[input$edgenum]
          if(sum(cutoff)==0){
            val[] <- 0
          }else{
            val[val < cutoff] <- 0
          }
          conns <- cbind(source=row(val)[val>0]-1, target=col(val)[val>0]-1, weight=val[val>0])
          if (nrow(conns) == 0){
            conns <- list(source=-1, target=-1, weight=0)
          }
          net <- list()
          net[["names"]] <- colnames(val)
          net[["links"]] <- conns
          #net[["groups"]] <- as.numeric(cutree(hclust(as.dist(1-val)) , k=input$con_knum))
          #net[["groups"]] <- as.numeric(cutree(hclust(as.dist(1-val)) , k=1))
          net[["groups"]] <- as.numeric(cutree(hc, k=input$con_knum))
          net[["titles"]] <- colnames(val)
          if(input$either=="sample"){
            net[["colors"]] <- colorx()
          }
          if(input$either=="probe"){
            net[["colors"]] <- colory()
          }
          
          net
        })
        
        #  Edge number slider with animation and threshold shown underneith.
        output$edge <- renderUI({
          data <- data()
          sliderInput(inputId = "edgenum",
                      label = "Edge Number:",
                      min = 0, max = (((length(data[1,]))^2)/2 - length(data[1,])/2), value = 1, step = 1,
                      animate=animationOptions(interval=2000, loop=FALSE))
        })
        output$gen_text <- renderText({
          cutoff <- cutoff()
          paste("Distance Threshold:  ",round(cutoff,4),sep="")
        })
        
        #  The network SVG
        output$svg <- renderUI({
          HTML(paste("<div id=\"net\" class=\"shiny-network-output\"><svg /></div>", sep=""))
        })
        
        #  The heatmap SVG
        output$heat <- renderUI({
          my_mat <- tmpdata()
          
          if(is.null(my_mat)){
            return()
          }
          else{
            gp <- ggheat(my_mat,input$tweak,colorx(),colory(),hc(),hc2(),input$color1,input$color2,input$color3,input$rainbow)
            svgjs <- grid2jssvg(gp)
            return(svgjs)
          }
          
        })
        
        hc <- reactive({
          hclust(dist(t(data()),method = input$dist_method), input$hc_method)
        })
        
        hc2 <- reactive({
          hclust(dist(data(),method = input$dist_method), input$hc_method)
        })
        
        output$dendro <- renderPlot({
          nc <- input$con_knum
          hc <- hc() 
          ## a smallish simple dendrogram
          dhc <- as.dendrogram(hc)
          cut <- cutree(hc,nc)[hc$labels[hc$order]]
          
          ## toy example to set colored leaf labels :
          local({
            colLab <<- function(n) {
              if(is.leaf(n)) {
                a <- attributes(n)
                i <<- i+1
                attr(n, "nodePar") <- c(a$nodePar, list(lab.col = mycols[cut[i]]))
                attr(n, "edgePar") <- c(a$nodePar, list(col = mycols[cut[i]]))
              }
              n
            }
            mycols <- grDevices::rainbow(nc)
            i <- 0
          })
          dL <- dendrapply(dhc, colLab)
          plot(dL) ## --> colored labels!
          #par(op)
        })
        
        colorx <- reactive({
          hc <- hc2()
          nc <- input$con_knum
          rainbow(nc, alpha=NULL)[cutree(hc,nc)[hc$labels[hc$order]]]
        })
        
        colory <- reactive({
          hc <- hc()
          nc <- input$con_knum
          rainbow(nc, alpha=NULL)[cutree(hc,nc)[hc$labels[hc$order]]]
        })
        
        colorxa <- reactive({
          hc <- hc2()
          nc <- input$con_knum
          rainbow(nc,alpha=NULL)[cutree(hc,nc)]
        })
        
        colorya <- reactive({
          hc <- hc()
          nc <- input$con_knum
          rainbow(nc,alpha=NULL)[cutree(hc,nc)]
        })
        
        output$color1_text <- renderText({
          input$color1
        })
        
        #  Close Button  
        observe({
          if (input$closebutton == 0)
            return()
          isolate({
            stopApp()
          })
        })
        
      })

    #myRunApp(app, ...)
    runApp(app)
  })
  


