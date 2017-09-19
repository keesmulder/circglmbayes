#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

shinyServer(function(input, output) {

  ### Argument names:
  ArgNames <- reactive({
    Names <- names(formals(input$readFunction)[-1])
    Names <- Names[Names!="..."]
    return(Names)
  })

  # Argument selector:
  output$ArgSelect <- renderUI({
    if (length(ArgNames())==0) return(NULL)

    selectInput("arg","Argument:",ArgNames())
  })

  ## Arg text field:
  output$ArgText <- renderUI({
    fun__arg <- paste0(input$readFunction,"__",input$arg)

    if (is.null(input$arg)) return(NULL)

    Defaults <- formals(input$readFunction)

    if (is.null(input[[fun__arg]]))
    {
      textInput(fun__arg, label = "Enter value:", value = deparse(Defaults[[input$arg]]))
    } else {
      textInput(fun__arg, label = "Enter value:", value = input[[fun__arg]])
    }
  })


  ### Data import:
  Dataset <- reactive({
    if (is.null(input$file)) {
      # User has not uploaded a file yet
      return(data.frame())
    }

    args <- grep(paste0("^",input$readFunction,"__"), names(input), value = TRUE)

    argList <- list()
    for (i in seq_along(args))
    {
      argList[[i]] <- eval(parse(text=input[[args[i]]]))
    }
    names(argList) <- gsub(paste0("^",input$readFunction,"__"),"",args)

    argList <- argList[names(argList) %in% ArgNames()]

    Dataset <- as.data.frame(do.call(input$readFunction,c(list(input$file$datapath),argList)))
    return(Dataset)
  })

  # Select variables:
  output$varselect <- renderUI({

    if (identical(Dataset(), '') || identical(Dataset(),data.frame())) return(NULL)

    # Variable selection:
    selectInput("vars", "Variables to use:",
                names(Dataset()), names(Dataset()), multiple =TRUE)
  })

  # Show table:
  output$table <- renderTable({

    if (is.null(input$vars) || length(input$vars)==0) return(NULL)

    return(Dataset()[,input$vars,drop=FALSE])
  })


  ### Download dump:

  output$downloadDump <- downloadHandler(
    filename = "Rdata.R",
    content = function(con) {

      assign(input$name, Dataset()[,input$vars,drop=FALSE])

      dump(input$name, con)
    }
  )

  ### Download save:

  output$downloadSave <- downloadHandler(
    filename = "Rdata.RData",
    content = function(con) {

      assign(input$name, Dataset()[,input$vars,drop=FALSE])

      save(list=input$name, file=con)
    }
  )

})
