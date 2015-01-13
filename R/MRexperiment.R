################################################################################
###   metagenomeSeq
################################################################################

.MGS_setSidebarPanel <- function(){
  sidebarPanel(
    h3("MRexperiment", align = "center"),
    HTML("<hr />"),
    checkboxInput("norm","Normalize data (CSS-Normalization)",value=TRUE),
    # checkboxInput("usePD","Include")
    numericInput("pd","Phenotype column:",1,min=1),
    selectInput("plot","Options for:",choices=c("Feature abundance","Heatmap","PCA/MDS","Diversity","OTU Richness"),
      selected="Feature abundance"),
    HTML("<hr />"),
    conditionalPanel(condition='input.plot=="Feature abundance"',
      h4("Feature abundance options", align = "center"),
      numericInput("featureIndex","Display feature (index):",1,min=1)
    ),
    conditionalPanel(condition='input.plot=="Heatmap"',
      h4("Heatmap options", align = "center"),
      numericInput("noFeatures","Number of features:",10,min=1,max =200),
      radioButtons("heatMethod","Choose features by:",c("Median Absolute Deviation"="mad","Variability"="sd"))
    ),
    conditionalPanel(condition='input.plot=="PCA/MDS"',
      h4("PCA/MDS options", align = "center"),
      radioButtons("pcaOrMds","PCA or MDS:",c("PCA"="TRUE","MDS"="FALSE"),selected="TRUE"),
      radioButtons("useDist","Make use of count distances:",c("False"="FALSE","True"="TRUE"),selected="FALSE"),
      conditionalPanel(condition = "input.useDist == 'TRUE'",
        selectInput("distance", "Distance:", 
          choices=c("euclidean","manhattan","canberra","bray",
            "kulczynski","jaccard","gower","altGower","morisita",
            "horn","raup","binomial","chao","cao"),selected="euclidean")
      ),
      numericInput('dimensionx', 'X-axis dimension:', 1,
         min = 1, max = 4),
      numericInput('dimensiony', 'Y-axis dimension:', 2,
         min = 1, max = 4)
    ),
    conditionalPanel(condition='input.plot=="Diversity"',
      h4("Diversity options", align = "center"),
      selectInput("diversity","Diversity index:",choices=c("shannon", "simpson", "invsimpson"))
    )
    # HTML("<hr />"),
    # actionButton("closebutton", "Stop/Save")
  )
}


.MGS_setMainPanel <- function(){
  mainPanel(
    tags$link(rel="stylesheet", type="text/css",
        href="/css-interactiveDisplay/interactiveDisplay.css"),
    tags$script(src="http://d3js.org/d3.v2.js"),
    .jstags(),
    .csstags(),
    shiny::tags$head(
      #  Insert what css you need in here, I just left in an example.
      shiny::tags$style(type='text/css', "
  
      heat {
        height: 800px;
      }
  
      net {
        height: 800px;
      }
  
      svg {
        height: 100vh;
      }
  
      ")
    ),
    .loading_gif(),
    tabsetPanel(
      tabPanel("Feature abundance",
               HTML("Manhattan plots of bacterial abundances where each bar represents the abundance in a sample."),
               HTML("<hr />"),
               plotOutput("plot2")),
      tabPanel("Heatmap",
               HTML("Heatmap of the abundances sorted by either MAD or SD"),
               HTML("<hr />"),
               plotOutput("plot1")),
      tabPanel("PCA/MDS",
               HTML("PCA or MDS plots (methods for dimensionality-reduction). Projecting samples onto two dimensions and coloring the samples by phenotype."),
               HTML("<hr />"),
               plotOutput("plot3")),
      tabPanel("Diversity",
               HTML("Shannon diversity index boxplots of the various stratified samples."),
               HTML("<hr />"),
               plotOutput("plot4")),
      tabPanel("OTU Richness",
               HTML("."),
               HTML("<hr />"),
               plotOutput("plot5"))
    ),
    tabsetPanel(
      tabPanel("Phenotype data",
               HTML("Phenotype information"),
               HTML("<hr />"),
               dataTableOutput("table2")
      ),
      tabPanel("Feature data",
               HTML("Feature information"),
               HTML("<hr />"),
               dataTableOutput("table1")
      )
    )
  )
}

setMethod("display",  
  signature(object = c("MRexperiment")), 
  function(object, ...){
    
    .usePackage('metagenomeSeq')
    .usePackage('vegan')
    .usePackage('gplots')

    app <- list(
      ui =
        bootstrapPage(
          .MGS_setSidebarPanel(),
          .MGS_setMainPanel()
        ),
      
      server = function(input, output){

        data = reactive({
                MRcounts(object,norm=input$norm,log=input$norm)
          })
        rawdata = reactive({
                MRcounts(object)
          })
        pd   = reactive({
                pd = pData(object)
                if(ncol(pd)==0){
                  return(as.matrix(rep(1,nrow(pd))))
                } else {
                  return(pd)
                }
          })

        output$plot1 <- renderPlot({
          
          #  Get data ready for ggplot
          mat = data()
          otusToKeep = which(rowSums(mat) > 0)
          otuStats = apply(mat[otusToKeep, ], 1, input$heatMethod)
          otuIndices = otusToKeep[order(otuStats, decreasing = TRUE)[1:input$noFeatures]]
          my_mat <- mat[otuIndices,]
          # melted <- melt(my_mat)
          # names(melted) <- c("Var1","Var2","value")
          
          # melted$Var1 <- factor(melted$Var1, rownames(my_mat))
          # melted$Var2 <- factor(melted$Var2, colnames(my_mat))
          trials = pd()[,input$pd]
          heatmapColColors=brewer.pal(12,"Set3")[as.integer(factor(trials))];
          heatmapCols = colorRampPalette(brewer.pal(9, "RdBu"))(50)
          gplots::heatmap.2(my_mat,trace="none",cexRow=.8,cexCol=.8,col = heatmapCols,ColSideColors = heatmapColColors)
          legend("left",fill=unique(heatmapColColors),legend=unique(trials))
        },height=850)
                
        output$plot2 <- renderPlot({
          mat = data()
          pdata  = factor(pd()[,input$pd])
          ylabel = ifelse(input$norm,yes=expression(bold("Log"[2]*" Abundance")),no="No. raw reads")
          plotFeature(mat,otuIndex = input$featureIndex,ylab=ylabel,main=rownames(mat)[input$featureIndex],
          classIndex = list(All_samples=1:ncol(mat)),col=pdata,font.lab=2,font.axis=1,sort=FALSE,xaxt="n")
          legend("topleft",legend=unique(pdata),fill=unique(pdata),box.col="NA")
        })
        output$plot3 <- renderPlot({
          useDist = input$useDist
          pdata  = factor(pd()[,input$pd])
          plotOrd(data(),n=min(100,nrow(data())),pch=21,bg=pdata,usePCA=input$pcaOrMds,
            comp=c(input$dimensionx,input$dimensiony),
            useDist=useDist,distfun=vegan::vegdist,dist.method=input$distance)
          legend("topleft",levels(pdata),fill=factor(levels(pdata)),box.col="NA")
        })
        output$plot4 <- renderPlot({
          pdata  = factor(pd()[,input$pd])
          mat = t(rawdata())
          H = vegan::diversity(mat,index=input$diversity)
          boxplot(H~pdata,ylab=paste(input$diversity,"diversity index"))
        })
        output$plot5 <- renderPlot({
          totalCounts = libSize(object)
          numFeatures=colSums(data()!=0)
          pdata = factor(pd()[,input$pd])
          plot(totalCounts, numFeatures,pch=21,xlab="Depth of coverage",
            ylab= "Number of detected features",bg=pdata)
          legend("topleft",levels(pdata),fill=factor(levels(pdata)),box.col="NA")
        })

        #  Example Table
        output$table1 <- renderDataTable({
          fd = fData(object)
          Index = 1:nrow(fd)
          Rownames = rownames(fd)
          fd = cbind(Index,Rownames,fd)
          return(fd)
        },options = list(iDisplayLength = 10))
                          
        #  Example Table
        output$table2 <- renderDataTable({
          return(pData(object))
        },options = list(iDisplayLength = 10))

        #  Close/Save Button, depending on your needs.
        # observe({
        #   if (input$closebutton == 0)
        #     return()
        #   isolate({
        #     stopApp(returnValue="Something you want to return to the R console")
        #   })
        # })     
      })
    runApp(app, ...)
  })
