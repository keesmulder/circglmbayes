#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

shinyUI(pageWithSidebar(

  # Header:
  headerPanel("Circular GLM"),

  # Input in sidepanel:
  sidebarPanel(

    tabsetPanel(

      tabPanel(
        "Load data",

        tags$style(type='text/css', ".well { max-width: 20em; }"),
        # Tags:
        tags$head(
          tags$style(type="text/css", "select[multiple] { width: 100%; height:10em}"),
          tags$style(type="text/css", "select { width: 100%}"),
          tags$style(type="text/css", "input { width: 19em; max-width:100%}")
        ), br(),


        radioButtons("datasource", "Choose data source",
                     choices = c("Use the example" = "example", "Load your own data" = "user")),
        # checkboxInput("useexample", "Use example dataset instead.", value = FALSE),


        conditionalPanel(
          "input.datasource == 'user'",
          # Select filetype:
          selectInput("readFunction", "Function to read data:", c(
            # Base R:
            "read.table",
            "read.csv",
            "read.csv2",
            "read.delim",
            "read.delim2",

            # foreign functions:
            "read.spss",
            "read.arff",
            "read.dbf",
            "read.dta",
            "read.epiiinfo",
            "read.mtp",
            "read.octave",
            "read.ssd",
            "read.systat",
            "read.xport",

            # Advanced functions:
            "scan",
            "readLines"
          ), selected = "read.csv"),

          # Upload data:
          fileInput("file", "Upload data-file:"),

          # Select whether to show all the advanced options.
          checkboxInput("advOptions", "Show advanced options", value = FALSE),

          conditionalPanel(
            "input.advOptions",

            # Argument selecter:
            htmlOutput("ArgSelect"),

            # Argument field:
            htmlOutput("ArgText"),

            # Variable selection:
            htmlOutput("varselect")
            #
            #       br(),
            #
            #       textInput("name","Dataset name:","Data")
            #
            #       downloadLink('downloadDump', 'Download source'),
            #       downloadLink('downloadSave', 'Download binary')
          )
        ), br()


      ),

      tabPanel(
        "Analysis",

        br(),

        # Outcome selection:
        htmlOutput("outcomeselect"),

        # Predictor selection:
        htmlOutput("predictorselect"),

        actionButton("run", "Run analysis")


      )
    )
  ),


  # Main:
  mainPanel(

    tabsetPanel(
      tabPanel("Full data",
        tableOutput("showdata")
      )
      ,
      tabPanel("Main results",
               verbatimTextOutput("textout")
      )
    )
  )
))
