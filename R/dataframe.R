.selDataTableOutput <- 
     function(outputId, ...,style="" ) 
{   
     origStyle<- c( 
         '<script src="shared/datatables/js/jquery.dataTables.min.js"></script>',
         '<script class="shiny-html-output" 
                  src= "/js/DTbinding.js"></script>',
         '<link rel = "stylesheet", 
                type = "text/css", 
                href = "shared/datatables/css/DT_bootstrap.css"></link>',
         '<style type="text/css">
                .rowsSelected td{
                background-color: rgba(112,164,255,0.2) 
                !important})  </style>',
         '<style type="text/css"> .selectable div table tbody tr{
                cursor: hand; cursor: pointer;}</style>',
         '<style type="text/css"> .selectable div table tbody tr td{
                -webkit-touch-callout: none;
                -webkit-user-select: none;
                -khtml-user-select: none;
                -moz-user-select: none;
                -ms-user-select: none;
                user-select: none;} </style>',
         '<style type="text/css">
                #myTable tfoot {display:table-header-group;}</style>')     
                  
    if(missing(style)){
        #mtcars
        tagList(
            singleton(
                tags$head(HTML(origStyle)
                )
            ),
            div(id = outputId, class = "shiny-datatable-output selectable")
        )        
    }else
    {
        #AnnotationHub
        tagList(
            singleton(
                tags$head(
                    HTML(origStyle),
                    HTML(style)
                )
            ),
            div(id = outputId, class = "shiny-datatable-output selectable")
        )        
    }
}

.dataFrame <- 
function(df, ..., summaryMessage = "", serverOptions = list(bSortClasses=TRUE,bRetrieve=TRUE),style="")
{  
    colNames <- colnames(df)
    app <- list(ui=pageWithSidebar(
        headerPanel("Data Tables binding"),
        sidebarPanel(
            tags$head(
                tags$style(type='text/css', ".span4 { max-width: 330px; }")
            ), 
            conditionalPanel(
                condition= "output.summary",
                strong(uiOutput('summary'))
            ),
            br(),
            actionButton("btnSend", "Send Rows"),
            em(p("Shift-Click to select multiple rows.")),
            br(),
            tags$button("Select All Rows", class="btn", id="select_all_rows"),
            em(p("Click to select all rows on page")),
            br(),
            tags$button("Deselect All Rows", class="btn", id="deselect_all_rows"),
            em(p("Click to deselect all rows on page"))
        ),
        mainPanel(
                if(length(style)!=1)
                    .selDataTableOutput(outputId="myTable", style=style)
                else{
                    .selDataTableOutput(outputId="myTable",...)}
        )
    ), server=function(input, output) {  
        output$myTable <- 
            renderDataTable({df}, 
                options = serverOptions                   
            )
        if (length(summaryMessage)!=1){
        output$summary <- renderUI({
            HTML(paste0(
                '<span class="shiny-html-output" >',summaryMessage[1],'</span> ',
                '<br>',
                '<span class="shiny-html-output" >',summaryMessage[2],'</span> ',
                '<br>',
                '<span class="shiny-html-output" >',summaryMessage[3],'</span> ',
                '<br>',
                '<span class="shiny-html-output" >',summaryMessage[4],'</span> ' ,
                '<br>',
                '<span class="shiny-html-output" >',summaryMessage[5],'</span> ' ,
                '<br>',
                '<span class="shiny-html-output" >',summaryMessage[6],'</span> ' ,
                '<br>'
                ))    
            })
        }
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
    runApp(app, ...)
}



setMethod("display", signature(object = c("data.frame")),
    function(object, ...)
{
    .dataFrame(df=object, ...)
})


#################################################
## testing:
## library(interactiveDisplay); df <- mtcars;
## foo <- interactiveDisplay:::.dataFrame(df)
## foo <- display(df)

## TODO: add support for trapping last usage (for cases where user
## accidently calls it without assignment like this : display(df)

