################################################################################
###   ExpressionSet
################################################################################

heatcolor1 <- function(inputId1) {
  tagList(
    shiny::tags$input(id = inputId1, class = "color", value = "EDF8B1",
    onchange = "window.Shiny.onInputChange('color1', this.color.toString())")
  )
}

heatcolor2 <- function(inputId2) {
  tagList(
    shiny::tags$input(id = inputId2, class = "color", value = "7FCDBB",
    onchange = "window.Shiny.onInputChange('color2', this.color.toString())")
  )
}

heatcolor3 <- function(inputId3) {
  tagList(
    shiny::tags$input(id = inputId3, class = "color", value = "2C7FB8",
    onchange = "window.Shiny.onInputChange('color3', this.color.toString())")
  )
}


.ES_setSidebarPanel <- function(){
  sidebarPanel(
    tableOutput("expinfo"),
    HTML("<hr />"),
    selectInput("either", "Network/Dendrogram View:  Sample or Probe:",
      choices = c("probe","sample")),
    HTML("<hr />"),
    uiOutput("choose_probe"),
    HTML("<hr />"),
    selectInput("order", "Show Top or Bottom Ranked,",
      choices = c("top","bottom")),
    selectInput("measure", "Based On:", choices = c("variance","average")),
    uiOutput("pmeancutoff"),
    uiOutput("smeancutoff"),
    HTML("<hr />"),
    sliderInput(inputId = "tweak",
                label = "Tweak Axis Label Font Size",
                min = -1, max = 1, value = 0, step = .1),
    HTML("<hr />"),
    #sliderInput(inputId = "con_knum",
    #            label = "Number of Clusters",
    #            min = 1, max = 100, value = 4, step = 1),
    uiOutput("edge"),
    uiOutput("gen_text"),
    HTML("<hr />"),
    numericInput("con_knum", "Number of Clusters:", 5),
    selectInput("hc_method", "Hierarchical Clustering Method",
      choices = c("ward", "single", "complete", "average", 
                  "mcquitty", "median", "centroid")),
    selectInput("dist_method", "Distance/Similarity Method",
    choices = c("euclidean", "maximum", "manhattan", 
                "canberra", "minkowski")),
    HTML("<hr />"),
    radioButtons('rainbow', 'Heatmap Color Scale',
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
      tabPanel("Heat Plot",
               HTML("Use the mouse to drag and pan the heatmap.  Use the 
                     mousewheel to zoom in/out."),
               HTML("<hr />"),
               uiOutput("heat")),
      tabPanel("Network View",uiOutput("svg")),
      tabPanel("Dendrogram",plotOutput("dendro"))
    ),
    tabsetPanel(
      tabPanel("GO Summary",
        HTML("Use the sidebar drop-down or simply click on probe nodes in the 
              force layout plot in the Network View tab to view a gene ontology 
              summary if available."),
        HTML("<hr />"),
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
            
            p <- p[seq_len(input$pmeancutoff)]
            s <- s[seq_len(input$smeancutoff)]
            
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
          tmpdata <- tmpdata()
          if(length(tmpdata)!=0){
            prChoices <- sort(rownames(tmpdata))
            names(prChoices) <- prChoices
            selectInput("probe", "Chosen Probe for GO Summary", prChoices)
          }
          else{
            return(NULL)
          }
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
              if(is(pkg,"ChipDb")){
                res <- select(pkg, input$probe,
                  c("ENTREZID","GENENAME"),
                  "PROBEID")
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
                res <- select(pkg, input$probe,
                  c("ENTREZID","GENENAME","GO"), "PROBEID")
                res2 <- head(select(GO.db, res$GO, "TERM", "GOID"))
                return(res2)
              }
              else{
                return(as.data.frame(
                  "Object does not have a ChipDb annotation"))
              }
            }
          }
        })
        
        #  Data for network view, sample or probe
        data <- reactive({
          data <- tmpdata()
          if(length(data)!=0){
            if(input$either=="sample"){
              data <- data
            }else{
              data <- t(data)
            }
            return(data)
          }
          else{
            return()
          }
        })
        
        #  Subset probes by average expression
        output$pmeancutoff <- renderUI({
          textInput("pmeancutoff", "Number of Probes to Display", 20)   
          #textInput("pmeancutoff", "Top Probes by Average Expression",
          #dim(exprs(object))[1])
        })
        
        #  Subset samples by average expression
        output$smeancutoff <- renderUI({   
          textInput("smeancutoff", "Number of Samples to Display",
            dim(exprs(object))[2])
        })
        
        #  This determines the distance threshold needed for the desired
        #  number of edges
        cutoff <- reactive({
          data <- data()
          if(length(data)!=0){
            #val <- cor(data())
            val <- dm(data())
            diag(val) <- NA
            val[lower.tri(val)] <- NA
            cutoff <- sort(val[!is.na(val)],decreasing=FALSE)[input$edgenum]
            return(cutoff)
          }
          else{
            return()
          }
        })
        
        #  This determines the distance threshold needed for the desired 
        #  number of edges
        cutoff_max <- reactive({
          data <- data()
          if(length(data)!=0){
            #val <- cor(data())
            val <- dm(data())
            diag(val) <- NA
            val[lower.tri(val)] <- NA
            cutoff <- 1
            edgenum <- 1
            while(!is.na(cutoff)){
              cutoff <- sort(val[!is.na(val)],decreasing=FALSE)[edgenum]
              edgenum <- edgenum + 1
            }
            return(edgenum - 2)
          }
          else{
            return()
          }
        })
        
        #  Build the network
        output$net <- reactive({
          data <- data()
          hc <- hc(data)
          data <- data[,hc$order]
          cutoff <- cutoff()
          if(length(data)!=0){
            #val <- cor(data())
            val <- dm(data())
            if (is.null(val)){
              return(list(names=character(), links=list(source=-1, target=-1)))
            }
            diag(val) <- NA
            val[lower.tri(val)] <- NA
            #cutoff <- sort(val[!is.na(val)],decreasing=TRUE)[input$edgenum]
            if(sum(cutoff)==0){
              val[] <- NA
            }else{
              val[val > cutoff] <- NA
            }
            conns <- cbind(source=row(val)[!is.na(val)]-1,
                           target=col(val)[!is.na(val)]-1,
                           weight=val[!is.na(val)])
            if (nrow(conns) == 0){
              conns <- list(source=-1, target=-1, weight=0)
            }
            net <- list()
            net[["names"]] <- colnames(val)#hc$labels[hc$order]
            net[["links"]] <- conns
            #net[["groups"]] <- as.numeric(cutree(hclust(as.dist(1-val)) ,
            #  k=input$con_knum))
            #net[["groups"]] <- as.numeric(cutree(hclust(as.dist(1-val)) , k=1))
            net[["groups"]] <- as.numeric(cutree(hc, k=input$con_knum))
            net[["titles"]] <- hc$labels[hc$order]
            net[["colors"]] <- 
              rainbow(input$con_knum,
                      alpha=NULL)[cutree(hc,input$con_knum)[colnames(val)]]
            #if(input$either=="sample"){
            #  net[["colors"]] <- colorx()
            #}
            #if(input$either=="probe"){
            #  net[["colors"]] <- colory()
            #}
            
            return(net)
          }
          else{
            return()
          }
        })
        
        #  Edge number slider with animation and threshold shown underneith.
        output$edge <- renderUI({
          data <- data()
          cutoff_max <- cutoff_max()
          if(length(data)!=0 && length(cutoff_max)!=0 ){
            sliderInput(inputId = "edgenum",
                        label = "Edge Number:",
                        #min = 0, 
                        #max = (((length(data[1,]))^2)/2 - length(data[1,])/2),
                        #value = 1, step = 1,
                        min = 0, max = cutoff_max, value = 1, step = 1,
                        animate=animationOptions(interval=2000, loop=FALSE))
          }
          else{
            return(NULL)
          }
        })
        output$gen_text <- renderText({
          cutoff <- cutoff()
          if(length(cutoff)!=0){
            paste("Distance Threshold:  ",round(cutoff,4),sep="")
          }
          else{
            return()
          }
        })
        
        #  The network SVG
        output$svg <- renderUI({
          HTML(paste(
            "<div id=\"net\" class=\"shiny-network-output\"><svg /></div>",
            sep=""))
        })
        
        #  The heatmap SVG
        output$heat <- renderUI({
          my_mat <- tmpdata()
          
          if(is.null(my_mat)){
            return()
          }
          else{
            gp <- ggheat(my_mat,exp(input$tweak),color_samples(),color_probes(),
                         hc(t(tmpdata())),hc(tmpdata()),
                         input$color1,input$color2,input$color3,input$rainbow)
            svgjs <- grid2jssvg(gp)
            return(svgjs)
          }
        })
        
        #Clustering
        hc <- function(d){
          hclust(dist(t(d),method = input$dist_method), input$hc_method)
        }
        
        #Distance Matrix
        
        dm <- function(d){
          as.matrix(dist(t(d), diag=TRUE, upper=TRUE, method=input$dist_method))
        }
                       
        #Color Dendrogram
        output$dendro <- renderPlot({
          nc <- input$con_knum
          hc <- hc(data())
          
          dhc <- as.dendrogram(hc)
          cut <- cutree(hc,nc)[hc$labels[hc$order]]
          
          colLab <- local({
              mycols <- grDevices::rainbow(nc)
              i <- 0
              function(n) {
                  if(is.leaf(n)) {
                      a <- attributes(n)
                      i <<- i+1
                      attr(n, "nodePar") <- c(a$nodePar,
                                              list(lab.col = mycols[cut[i]]))
                      attr(n, "edgePar") <- c(a$nodePar,
                                              list(col = mycols[cut[i]]))
                  }
                  n
              }
          })
          dL <- dendrapply(dhc, colLab)
          plot(dL)
        })
        
        #More cluster group coloring for x/y axis and network nodes
        color_samples <- reactive({
          d <- tmpdata()
          hc <- hc(d)
          nc <- input$con_knum
          rainbow(nc, alpha=NULL)[cutree(hc,nc)[hc$labels[hc$order]]]
        })
        
        color_probes <- reactive({
          d <- tmpdata()
          hc <- hc(t(d))
          nc <- input$con_knum
          rainbow(nc, alpha=NULL)[cutree(hc,nc)[hc$labels[hc$order]]]
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
  


