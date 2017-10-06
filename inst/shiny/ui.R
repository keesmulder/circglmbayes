
library(shiny)
library(shinydashboard)




shinyUI(dashboardPage(

  # Header:
  dashboardHeader(title = "Circular GLM"),

  # Input in sidepanel:
  dashboardSidebar(

    sidebarMenu(

      h3(" Input"),
      menuItem("Load data", tabName = "loaddata"),
      menuItem("Analysis options", tabName = "analysisopts"),
      menuItem("View dataset",     tabName = "fulldata"),
      actionButton("run", "Run analysis")

    ),
    sidebarMenuOutput("outputmenu")
  ),


  # Main:
  dashboardBody(

    tabItems(
      tabItem(
        tabName = "loaddata",

              fluidRow(
                column(width = 6,
              h3("Data selection"), br(),

              tags$style(type='text/css', ".well { max-width: 20em; }"),
              # Tags:
              tags$head(
                tags$style(type="text/css", "select[multiple] { width: 100%; height:10em}"),
                tags$style(type="text/css", "select { width: 100%}"),
                tags$style(type="text/css", "input { width: 19em; max-width:100%}")
              ),


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

                )
                       )
              ),
                column(width = 6,

                # Outcome selection:
                htmlOutput("outcomeselect"),

                # Predictor selection:
                htmlOutput("predictorselect")
              )
                       )

      ),



      tabItem(
        tabName = "analysisopts"
      ),


      tabItem(tabName = "fulldata",
              br(), br(),
              dataTableOutput("showdata")
      ),

      tabItem(tabName = "results",
              # h3("Overview"),
              textOutput("textoverview"),
              h3("Coefficients"),
              tableOutput("coeftable"),
              h3("Model fit (Information Criteria)"),
              tableOutput("ICtable")
              # h3("Hypothesis tests"),
              # tableOutput("bftables")

      ),
      tabItem(tabName = "plotout",
              # box(
              #   title = "MCMC Chains",
              #   plotOutput("traceplot")
              # ),
              box(
                title = "Prediction plot",
                plotOutput("predictplot")
              ),
              box(
                title = "MCMC Results",
                plotOutput("tracestackplot")
              ),
              box(
                title = "Meancompare",
                plotOutput("meancompplot")
              ),
              box(
                title = "Mean Boxplot",
                plotOutput("meanboxplot")
              )
      ),
      tabItem(tabName = "rtext",
              fluidRow(
                column(
                  width = 6,
                  box(
                    width = NULL,
                    title = "Overview",
                    verbatimTextOutput("basetextprint")
                  ),
                  box(
                    width = NULL,
                    title = "MCMC Summary",
                    verbatimTextOutput("mcmctextprint")
                  ),
                  box(
                    width = NULL,
                    title = "Bayes Factors and posterior model probabilities",
                    verbatimTextOutput("bftextprint")
                  )
                ),
                column(
                  width = 6,
                  box(
                    width = NULL,
                    title = "Full results object",
                    verbatimTextOutput("alltextprint")
                  )
                )
              )
      )
    )
  )
))
