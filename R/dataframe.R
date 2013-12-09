#needHack <- TRUE # change to false if problem referenced here is solved:
# https://groups.google.com/forum/#!topic/shiny-discuss/Lhdi7A_csR4

# Note that if needHack is true, you'll have to (I think) restart the shiny app
# in order for changes to the js file to take effect. 

#source("global.R")


.selDataTableOutput <- function (outputId) 
{
    ## Temp. fragment of HTML
    FRAG = tagList(
      HTML("<script type='text/javascript'>"),
      includeHTML(system.file("www", "DTbinding.js",
                              package="interactiveDisplay")),
      HTML("</script>")
      )
    
  tagList(singleton(tags$head(tags$link(rel = "stylesheet", 
    type = "text/css", href = "shared/datatables/css/DT_bootstrap.css"),
    tags$style(type="text/css", ".rowsSelected td{
               background-color: rgba(112,164,255,0.2) !important}"),
    tags$style(type="text/css", ".selectable div table tbody tr{
               cursor: hand; cursor: pointer;}"),
    tags$style(type="text/css",".selectable div table tbody tr td{
               -webkit-touch-callout: none;
               -webkit-user-select: none;
               -khtml-user-select: none;
               -moz-user-select: none;
               -ms-user-select: none;
               user-select: none;}"),                          
    tags$script(src = "shared/datatables/js/jquery.dataTables.min.js"), 
    tags$script(src = "shared/datatables/js/DT_bootstrap.js"),
    FRAG)),
    div(id = outputId, class = "shiny-datatable-output selectable"))
}


.dataFrame <- function(df, ...){
colNames <- colnames(df)
app <- list(ui=
pageWithSidebar(
  
  headerPanel("Data Tables binding"),
  sidebarPanel(
    p("Shift-Click to select multiple rows."),
    actionButton("btnSend", "Send Rows"),
    tags$button("Select All Rows", class="btn", id="select_all_rows"),
    tags$button("Deselect All Rows", class="btn", id="deselect_all_rows")
  ),
  
  mainPanel(
    .selDataTableOutput("myTable")
  )
),

server=
function(input, output) {
  
  output$myTable <- renderDataTable({df}, options = list(bSortClasses = TRUE))
  
  # observe({
  #   print("Click Event")
  #   print(input$myTable)
  #   })

 observe({
    if (input$btnSend > 0)
      isolate({
          #print(input$myTable)
          dfVec <- input$myTable
          df <- as.data.frame(matrix(data=dfVec, ncol=dim(df)[2], byrow=TRUE))
          names(df) <- colNames
          stopApp(return = df)
      })
    })
})


## here is how we run it  (wrap in method for data.frame when done)
runApp(app, ...)
}

setMethod("display",
          signature(object = c("data.frame")),
          function(object, ...){.dataFrame(df=object, ...)})


#################################################
## testing:
## library(interactiveDisplay); df <- mtcars;
## foo <- interactiveDisplay:::.dataFrame(df)
## foo <- display(df)

## TODO: add support for trapping last usage (for cases where user
## accidently calls it without assignment like this : display(df)
