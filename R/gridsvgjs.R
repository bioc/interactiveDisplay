setGeneric("gridsvgjs", function(object)
  standardGeneric("gridsvgjs")
)

setMethod("gridsvgjs", 
          signature(object = c("ANY")),
          function(object){
            
            app <- list(
              ui =
                bootstrapPage(
                  
                  .jstags(),
                  
                  mainPanel(
                    
                    shiny::tags$head(
                      shiny::tags$style(type='text/css', 
                        ".span8 { width: 100%; align: right; }")
                    ),
                                        
                    uiOutput("svgplot")
                  )
                ),
              
              server = function(input,output) {
                
                
                
                output$svgplot <- renderUI({    
                  
                  
                  jscode <- "
                            <script type='text/javascript'>
                            $(document).ready(function() {
                              $('svg').svgPan('viewport');
                            });
                            </script>
                            "
                  
                  print(object)
                  mysvg <- grid.export()
                  mysvg2 <- saveXML(mysvg$svg[["g"]])
                  mysvg3 <- sub("<g transform=",
                                "<g id='viewport' transform=",
                                mysvg2)                 
                  mysvg4 <- sub(">NA<","><",mysvg3)
                  htmlxml <- HTML(paste("<svg xmlns='http://www.w3.org/2000/svg'
                               xmlns:xlink='http://www.w3.org/1999/xlink' 
                               version='1.1' width='100%' height='100%'>",
                               jscode,mysvg4,"</svg>",sep=""))
                  htmlxml
                  
                })
              }
            )
            runApp(app)
          })
