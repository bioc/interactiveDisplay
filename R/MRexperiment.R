################################################################################
###   metagenomeSeq
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

.MGS_setSidebarPanel <- function(){
  sidebarPanel(
    h3("MRexperiment", align = "center"),
    HTML("<hr />"),
    checkboxInput("noheat","Suppress Heatmap"),
    checkboxInput("flip","Transpose Heatmap"),
    HTML("<hr />"),
    radioButtons('rainbow', 'Heatmap Color Scale',
                 c('Rainbow'='default',
                   'Three Color'='tri'),
                 'Rainbow'),
    heatcolor1('exampleTextarea1'),
    heatcolor2('exampleTextarea2'),
    heatcolor3('exampleTextarea3'),
    HTML("<hr />"),
    actionButton("closebutton", "Stop/Save")
  )
}

.MGS_setMainPanel <- function(){
  mainPanel(
    tags$link(rel="stylesheet", type="text/css",
        href="/css/interactiveDisplay.css"),
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
      tabPanel("Plot 1",
               HTML("On-screen Documentation"),
               HTML("<hr />"),
               plotOutput("plot1")),
      tabPanel("Plot 2",
               HTML("On-screen Documentation"),
               HTML("<hr />"),
               plotOutput("plot2")),
      tabPanel("Plot 3",
               HTML("On-screen Documentation"),
               HTML("<hr />"),
               plotOutput("plot3"))
    ),
    tabsetPanel(
      tabPanel("Table 1",
               HTML("On-screen Documentation"),
               HTML("<hr />"),
               dataTableOutput("table1")
      ),
      tabPanel("Table 2",
               HTML("On-screen Documentation"),
               HTML("<hr />"),
               dataTableOutput("table2")
      )
    )
  )
}

setMethod("display",  
  signature(object = c("MRexperiment")), 
  function(object, ...){
    
    .usePackage('metagenomeSeq')
       
    app <- list(
      ui =
        bootstrapPage(
          .MGS_setSidebarPanel(),
          .MGS_setMainPanel()
        ),
      
      server = function(input, output){
        
        tmpdata <- reactive({
          #  Generate some fake example data.
          m <- matrix(rnorm(2600),nrow=26)
          rownames(m) <- replicate(26,paste0(sample(letters,1),sample(1000,1)))
          colnames(m) <- replicate(100,paste0(sample(LETTERS,1),sample(100000,1)))
          return(m)
        })
                
        output$plot1 <- renderPlot({
          
          #  Get data ready for ggplot
          my_mat <- tmpdata()
          melted <- melt(my_mat)
          names(melted) <- c("Var1","Var2","value")
          
          melted$Var1 <- factor(melted$Var1, rownames(my_mat))
          melted$Var2 <- factor(melted$Var2, colnames(my_mat))
          
          #  For that suppression option.
          if(is.null(my_mat) || input$noheat == TRUE){
            return()
          }
          #  This looks weird but it's a little trick to keep the color picker
          #  from refreshing/loading unnecessarily...
          else{    
            if(input$rainbow == 'default'){
              isolate({
                c1 <- input$color1
                c2 <- input$color2
                c3 <- input$color3
              })
            }
            else{
              c1 <- input$color1
              c2 <- input$color2
              c3 <- input$color3
            }
          }
          
          #  Default colorblind friendly colors.
          if(length(c1)==0){
            c1 <- "EDF8B1"
          }
          if(length(c2)==0){
            c2 <- "7FCDBB"
          }
          if(length(c3)==0){
            c3 <- "2C7FB8"
          }
          
          #  Full "rainbow" spectrum
          if(input$rainbow=='default'){
            myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
          }
          
          #  Chosen tricolor palette
          if(input$rainbow=='tri'){
            myPalette <- colorRampPalette(rev(c(paste("#",toupper(c1),sep=""),
                                                paste("#",toupper(c2),sep=""),
                                                paste("#",toupper(c3),sep=""))))
          }
          
          #  A lot of ggplot formatting...
          #  There are many excellent heatplot libraries, but to have more
          #  control I decided to built it from the ground up.
          gp <- ggplot(melted, aes(x = Var2, y = Var1, fill = value))
          gp <- gp + geom_tile()
          gp <- gp + coord_fixed()
          gp <- gp + scale_fill_gradientn(colours = myPalette(100))
          gp <- gp + scale_x_discrete(expand = c(0, 0))
          gp <- gp + scale_y_discrete(expand = c(0, 0))
          gp <- gp + coord_equal()
          gp <- gp + theme_bw()
          gp <- gp + xlab("X label")
          gp <- gp + ylab("Y label")
          gp <- gp + theme(axis.ticks = element_blank(),
                           axis.text.x = element_text(angle = -45,
                                                      hjust = 0))
          if(input$flip==TRUE){
            gp <- gp + coord_flip() 
          }
          
          return(print(gp))
        })
                
        output$plot2 <- renderPlot({
          #Empty!
        })
        
        output$plot3 <- renderPlot({
          #Empty!
        })
        
        #  Example Table
        output$table1 <- renderDataTable({
          return(mtcars)
        })
        
        #  Example Table
        output$table2 <- renderDataTable({
          return(iris)
        })
                             
        #  Close/Save Button, depending on your needs.
        observe({
          if (input$closebutton == 0)
            return()
          isolate({
            stopApp(returnValue="Something you want to return to the R console")
          })
        })     
      })
    runApp(app, ...)
  })
