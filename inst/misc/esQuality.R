#####################################################################################################
###   SummarizedExperiment 2
#####################################################################################################

setGeneric("esQuality", function(object, ...) {
  standardGeneric("esQuality")
})

setMethod("esQuality", 
          signature(object = "ExpressionSet"), 
          function(object, ...){
            
            require(altcdfenvs)
            require("hgu95av2cdf")
            require("hgu95av2.db")
            
            
            app <- list(
              ui =
                bootstrapPage(
                  
                  tagList(

                  ),
                  
                  h3("Expression Set Quality Check"),
                  
                  sidebarPanel(
                    uiOutput("choose_sample"),
                    HTML("<hr />"),
                    uiOutput("cutoff")
                  ),
                  
                  mainPanel(
                    
                    shiny::tags$head(
                      #shiny::tags$style(type='text/css', ".span8 { width: 800px; height: 800px }"),
                      #shiny::tags$style(type='text/css', ".levels { width: 800px; height: 800px }")
                    ),
                    
                    plotOutput("levels")
                  )
                ),
              
              server = function(input, output){
                
                placeholder <- wrapCdfEnvAffy(hgu95av2cdf, 640, 640, "hgu95av2")
                probe_list <- rownames(exprs(object))
                
                
                #  Sets max position for the silder
                max_end <- reactive({
                  mat <- mat()
                  round(max(mat))
                })
                
                #  Render the plot range window slider.     
                output$cutoff <- renderUI({
                  max_end <- max_end()
                  sliderInput(inputId = "cutoff",
                              label = "Cutoff:",
                              min = 0, max = max_end, value = 0, step = 1
                  )
                })
                

                #  Render the choose sample dropdown.
                sl <- colnames(exprs(object))
                output$choose_sample <- renderUI({
                  sampleChoices <- sl
                  names(sampleChoices) <- sl
                  selectInput("sample", "Sample", sampleChoices)
                })
                
                
                
                mat <- reactive({
                  
                  mat <- matrix(0,nrow=640,ncol=640)
                  
                  for(i in probe_list){
                    current_probe <- i
                    probe_i <- indexProbes.CdfEnvAffy(placeholder, "pm", i)[[1]]
                    xy <- index2xy (placeholder,probe_i)
                    
                    for(j in 1:length(xy[,1])){
                      mat[as.numeric(xy[j,][1]),as.numeric(xy[j,][2])] <- exprs(object)[current_probe,input$sample]
                    }
                  }
                  print(sum(mat > 0))
                  mat
                  
                })
                
                
                output$levels <- renderPlot({
                  
                  mat <- mat()
                  plot(which(mat>input$cutoff, arr.ind = T),asp=1)
                  
                })
                
              }
            )
            #runApp(app)
            myRunApp(app, ...)
          })
