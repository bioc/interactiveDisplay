#####################################################################################################
###   Annotation Hub
#####################################################################################################

setMethod("display", 
          signature(object = "AnnotationHub"), 
          function(object, ...){
            
            app <- list(
              ui = 
                #includeHTML("inst/www/js/myselect2.js"),
                bootstrapPage(
                #includeHTML("inst/www/js/myselect2.js"),
                h3("AnnotationHub"),
                sidebarPanel(
                  #includeHTML(system.file("www/js", "tools.js", package="interactiveDisplay")),
                  includeHTML("inst/www/js/tools.js"),
                  uiOutput("tool"),
                  tableOutput("info"),
                  tableOutput("info2"),
                  #tableOutput("info3"),
                  #tableOutput("info4"),
                  actionButton("savebutton", "Stop App")
                ),
                mainPanel(
                
                )
              ),
              
              server = function(input, output){
                
                #  Render the UCSC dropdown
                output$tool <- renderUI({
                  ucsc_vec <- names(object)
                  names(ucsc_vec) <- ucsc_vec
                  selectInput("tool", label = "Tool:",choices = ucsc_vec, selected = NULL, multiple = FALSE)
                })
                
                #  Table of experiment info if present
                output$info <- renderTable({
                  info <- t(as.data.frame(list(
                    #"filters" = filters(object),
                    "length" = length(object),
                    "hubUrl" = hubUrl(object),
                    "hubCache" = hubCache(object),
                    "hubResource" = hubResource(object),
                    #"possibleDates" = possibleDates(object),
                    "snapshotDate" = snapshotDate(object),
                    "snapshotUrl" = snapshotUrl(object),
                    #"snapshotUrls" = snapshotUrls(object),
                    "snapshotVersion" = snapshotVersion(object)
                    #"cols" = cols(object),
                    #"keytypes" = keytypes(object)
                    )))
                  colnames(info) <- c("Info")
                  info
                })
                
                #  Table of experiment info if present
                output$info2 <- renderTable({
                  info <- as.data.frame(
                    keytypes(object)
                  )
                  info
                })
                
                ##  Table of experiment info if present
                #output$info3 <- renderTable({
                #    info <- as.data.frame(
                #    head(snapshotUrls(object))
                #  )
                #  #colnames(info) <- c("Info")
                #  info
                #})
                
                ##  Table of experiment info if present
                #output$info4 <- renderTable({
                #  info <- as.data.frame(
                #    possibleDates(object)
                #  )
                #  #colnames(info) <- c("Info")
                #  info
                #})
                
                #  Save Button  
                observe({
                  if (input$savebutton == 0)
                    return()
                  isolate({
                    stopApp(returnValue=NULL)
                  })
                })
                
                
              }  
            )
            myRunApp(app, ...)
          })
