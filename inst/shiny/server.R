#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(foreign)
library(CircGLMBayes)
library(shinyjs)
library(dplyr)

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

    # Check if we should load the example
    if (input$datasource == 'example') {
      return(essbhv[1:100, ])

    # Check if user has not uploaded a file yet
    } else if (is.null(input$file)) {
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


  # Select Outcome:
  output$outcomeselect <- renderUI({

    if (identical(Dataset(), '') || identical(Dataset(),data.frame())) return(NULL)

    # Variable selection:
    selectInput("outcome", "Select the outcome:",
                names(Dataset()), names(Dataset()), multiple = FALSE)
  })

  # Select predictors:
  output$predictorselect <- renderUI({

    if (identical(Dataset(), '') || identical(Dataset(),data.frame())) return(NULL)

    # Predictors can be any variable that is not the outcome.
    predOpts <- names(Dataset())[names(Dataset()) != input$outcome]

    selectInput("predictors", "Select the predictors:",
                choices = predOpts, multiple = TRUE)
  })


  # Show full data:
  output$showdata <- renderTable({

    if (identical(Dataset(), '') || identical(Dataset(),data.frame())) return(NULL)

    return(Dataset())
  })







  ########### ANALYSIS ####################
  # Model run:
  getModel <- reactive({

    # If the run button is pressed, invalidate the model and run this again.
    input$run

    # However, isolate so that the model is not re-run when anything else changes.
    isolate({

      dat <- Dataset()

      if (length(input$predictors) == 0) return("No predictors were selected.")

      mod <- circGLM(th = dat[, input$outcome], X = dat[, input$predictors, drop = FALSE])
    })

    return(mod)
  })



  output$coeftable <- renderTable(coef(getModel()))
  output$bftables  <- renderTable(BF.circGLM(getModel()))


  # R Verbatim text output.
  output$basetextprint <- renderPrint(print(getModel(),                       digits = input$digits))
  output$mcmctextprint <- renderPrint(print(mcmc_summary.circGLM(getModel()), digits = input$digits))
  output$bftextprint   <- renderPrint(print(BF.circGLM(getModel()),           digits = input$digits))
  output$alltextprint  <- renderPrint(print(getModel(), type = 'all',         digits = input$digits))





  ### Download dump:
#
#   output$downloadDump <- downloadHandler(
#     filename = "Rdata.R",
#     content = function(con) {
#
#       assign(input$name, Dataset()[,input$vars,drop=FALSE])
#
#       dump(input$name, con)
#     }
#   )
#
#   ### Download save:
#
#   output$downloadSave <- downloadHandler(
#     filename = "Rdata.RData",
#     content = function(con) {
#
#       assign(input$name, Dataset()[,input$vars,drop=FALSE])
#
#       save(list=input$name, file=con)
#     }
  # )

})
