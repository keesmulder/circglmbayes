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

shinyServer(function(input, output, session) {




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
                names(Dataset()), names(Dataset())[1], multiple = FALSE)
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
  output$showdata <- renderDataTable({

    if (identical(Dataset(), '') || identical(Dataset(),data.frame())) return(NULL)

    return(Dataset())
  })



  output$outputmenu <- renderMenu({
    if (rvs$modelWasRun) {

      sidebarMenu(
        id = "outputsidebar",
        h3(" Output"),
        menuItem("Results",     tabName = "results"),
        menuItem("Plots",       tabName = "plotout"),
        menuItem("Text Output", tabName = "rtext"),

        h3(" Output options"),
        numericInput("digits", "Digits", 2, 0, 8)
      )
    } else {
      sidebarMenu("")
    }
  })

  # Set up a reactive value which check if the model was run.
  rvs <- reactiveValues(modelWasRun = FALSE)
  observeEvent(input$run, {

    # Save the fact that the model was run so we can open the output menu.
    rvs$modelWasRun <- TRUE

    updateTabItems(session, inputId = "outputsidebar", selected = "results")

    # Run the model to make sure that the output is generated.
    getModel()
  })

  ########### ANALYSIS ####################
  # Model run:
  getModel <- reactive({

    if (!rvs$modelWasRun) return("Model was not yet run.")

    # If the run button is pressed, invalidate the model and run this again.
    input$run

    # However, isolate so that the model is not re-run when anything else changes.
    isolate({

      thisFormula <- as.formula(paste(input$outcome, "~", paste(c(1, input$predictors), collapse = " + ")))

      print(thisFormula)

      mod <- circGLM(formula = thisFormula, data = Dataset())
      # mod <- circGLM(th = dat[, input$outcome, drop = FALSE], X = dat[, input$predictors, drop = FALSE])
    })

    return(mod)
  })

  output$textoverview <- renderText(paste("MCMC run for", getModel()$TotalIts, "iterations, of which", getModel()$SavedIts, "were used."))

  output$coeftable <- renderTable(coef(getModel()), rownames = TRUE, digits = reactive(input$digits))
  output$bftables  <- renderTable(BF.circGLM(getModel()), rownames = TRUE, digits = reactive(input$digits))
  output$ICtable  <- renderTable(IC_compare.circGLM(getModel()), rownames = TRUE, digits = reactive(input$digits))



  # PLOT OUTPUTS
  # output$traceplot      <- renderPlot(plot(getModel(), type = 'trace'))
  output$tracestackplot <- renderPlot(plot(getModel(), type = 'tracestack'))
  output$predictplot    <- renderPlot(plot(getModel(), type = 'predict'))
  output$meancompplot   <- renderPlot(plot(getModel(), type = 'meancompare'))
  output$meanboxplot    <- renderPlot(plot(getModel(), type = 'meanboxplot'))

  # R Verbatim text output.
  output$basetextprint <- renderPrint(print(getModel(),                       digits = input$digits))
  output$mcmctextprint <- renderPrint(print(mcmc_summary.circGLM(getModel()), digits = input$digits))
  output$bftextprint   <- renderPrint(print(BF.circGLM(getModel()),           digits = input$digits))
  output$alltextprint  <- renderPrint(print(getModel(), type = 'all',         digits = input$digits))

})
