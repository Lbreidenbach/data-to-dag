library(shiny)
library(DiagrammeR)
library(HydeNet)
library(rjags)
library(MatchIt)
library(ggplot2)
library(plyr)
library(EmpiricalBrownsMethod)
library(dplyr)
library(gtable)
library(grid)
library(GGally)
library(ggridges)
library("reshape2")
library(gridExtra)
source("plot_dag_ui.R")

node_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(3, 
           selectInput(ns("node_dist"), 
                       label = paste0("Choose distribution of ", id), 
                       choices = c("binary", "continuous"))
    ),
    column(9,
           uiOutput(ns("ui_placeholder"))
    )
  )
} 

node_server <- function(input, output, session) {
  return_value <- reactive({input$inner_element})
  ns <- session$ns
  output$ui_placeholder <- renderUI({
    type <- req(input$node_dist)
    if(type == "binary") {
      numericInput(ns("inner_element"), "Prevalence (decimal from 0-1):", 0)
    } else if (type == "continuous") {
      numericInput(ns("inner_element"), "Mean:", 0)
    }
  })
  
  ## if we later want to do some more sophisticated logic
  ## we can add reactives to this list
  list(return_value = return_value) 
}

ui <- fluidPage(
  titlePanel(h1("Bias Simulator")),
  
  sidebarLayout(position = "right",
    sidebarPanel("sidebar panel"),
    mainPanel("main panel",
              # Copy the line below to make a text input box
              textInput("dag_text", label = h3("Text input"), value = "~a|b + c|b"),
              p("Nodes in the DAG"),
              selectInput("inSelect", "Select input",
                          c("Item A")),
              p("details of node"),
              div(id="placeholder"),
              verbatimTextOutput("out"),
              hr(),
              grVizOutput("value"))
  
  )
)

# Define server logic ----
server <- function(input, output, session) {
  handler = reactiveVal(list())
  
  output$value = renderGrViz({ dag_ui(paste(input$dag_text)) })
  observe({
    x <- get_nodes(input$dag_text)
    
    # Can use character(0) to remove all choices
    if (is.null(x))
      x <- character(0)
    
    # Can also set the label and select items
    updateSelectInput(session, "inSelect",
                      label = paste("Select node to define"),
                      choices = x,
                      selected = tail(x, 1)
    )
  })
  observeEvent(input$inSelect, {
    new_id <- paste("node", input$inSelect, sep = "_")
    insertUI(
      selector = "#placeholder",
      where = "beforeBegin",
      ui = node_ui(new_id)
    )
    handler_list <- isolate(handler())
    new_handler <- callModule(node_server, new_id)
    handler_list <- c(handler_list, new_handler)
    names(handler_list)[length(handler_list)] <- new_id
    handler(handler_list)
  })
  output$out <- renderPrint({
    lapply(handler(), function(handle) {
      handle()
    })
  })  
}

# Run the app ----
shinyApp(ui = ui, server = server)

