.selDataTableOutput <- 
    function(outputId) 
{
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
    tags$style(type="text/css", '#myTable tfoot {display:table-header-group;}'), 
    tags$script(src = "shared/datatables/js/jquery.dataTables.min.js"), 
    tags$script(src = "shared/datatables/js/DT_bootstrap.js"),
    tags$script(src = "js/DTbinding.js"))),
    div(id = outputId, class = "shiny-datatable-output selectable"))
}


.dataFrame <- 
    function(df, ...)
{
    colNames <- colnames(df)
    colwidths <- c("10%", "20%", "30%", "2%", "10%","30%")
    aoColumnDefs <- list(NULL)
    for(i in 1:ncol(df)){
        column <- list(sWidth=colwidths[i],  aTargets=list(i-1))
        aoColumnDefs[[i]] <- column
    }
    app <- list(ui=pageWithSidebar(
        headerPanel("Data Tables binding"),
        sidebarPanel(
            tags$head(
                tags$style(type="text/css", "select { max-width: 200px; }"),
                tags$style(type="text/css", "textarea { max-width: 185px; }"),
                tags$style(type="text/css", ".jslider { max-width: 200px; }"),
                tags$style(type='text/css', ".well { max-width: 310px; }"),
                tags$style(type='text/css', ".span4 { max-width: 310px; }")
            ), 
            helpText("Aim to display summary statistic here! "),
            actionButton("btnSend", "Send Rows"),
            p("Ctrl-Click to select multiple rows."),
            br(),
            tags$button("Select All Rows", class="btn", id="select_all_rows"),
            br(),br(),
            tags$button("Deselect All Rows", class="btn", id="deselect_all_rows")
        ),
        mainPanel(
            .selDataTableOutput("myTable")
        )
    ),
    server=
        function(input, output) 
    {  
        output$myTable <- renderDataTable({df}, 
	    options = list(
                "sDom" = '<"top"i>rt<"top"f>lt<"clear">',
                bSortClasses = TRUE,
                bRetrieve= TRUE,
                #for fixing width of columns  
                bAutoWidth =  FALSE,
                aoColumnDefs=aoColumnDefs,
                aoColumns=NULL,                               
                #for pagination
                aLengthMenu = c(1000, 5000, "All"),
                iDisplayLength = 1000,
                #for searching   
                oSearch = list(
                    sSearch= "",
                    bSmart = TRUE,
                    bRegex = FALSE,
                    bCaseInsensitive=TRUE
                )  
            )
        )  
        observe({
            if(input$btnSend > 0)
                isolate({
                   #print(input$myTable)
                   dfVec <- input$myTable
                   df <- as.data.frame(matrix(data=dfVec, ncol=dim(df)[2],
                       byrow=TRUE))
                   names(df) <- colNames
                   stopApp(returnValue = df)
                })
            })
    })
## then run it...
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

