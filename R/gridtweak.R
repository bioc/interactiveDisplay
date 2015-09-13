
setGeneric("gridtweak", function(...)
  standardGeneric("gridtweak")
)

setMethod("gridtweak", 
          signature(),
          function(...){   
            app <- list(
              
              ui =
                sidebarLayout(
                  sidebarPanel(
                    shiny::tags$textarea(id="inbox", rows=4, cols=80, "gp <- ggplot() +
geom_point(data=iris,aes(x=Sepal.Length,y=Petal.Length))
print(gp)"),
                    actionButton("gobutton", "Refresh Plot"),
                    shiny::tags$style(type='text/css', "#inbox { width: 1000px; }")
                  ),
                  mainPanel( 
                    shiny::tags$head(
                      shiny::tags$style(type='text/css',
                                        ".span8 { width: 100%; align: right; }")
                    ),
                    .jstags(),
                    uiOutput("svgplot")
                  )
                ),
              
              server = function(input,output) {
                output$svgplot <- renderUI({
                  if (input$gobutton == 0)
                    return()
                  isolate({
                    jscode <- "
                    <script type='text/javascript'>
                    $(document).ready(function() {
                    $('svg').svgPan('viewport');
                    });
                    </script>
                    "
                    eval(parse(text=input$inbox))
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
                })
              }
                  )
            interactiveDisplayBase::.runApp(app, ...)
          })
