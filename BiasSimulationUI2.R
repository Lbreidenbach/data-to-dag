library(shiny)
library(DiagrammeR)
library(HydeNet)
library(rjags)
library(MatchIt)
library(ggplot2)
library(plyr)
#library(EmpiricalBrownsMethod)
library(dplyr)
library(gtable)
library(grid)
library(GGally)
library(ggridges)
library("reshape2")
library(gridExtra)
source("plot_dag_ui.R")
library(purrr)

textInput3<-function (inputId, label, value = "",...) {
  div(style="display:inline-block",
      tags$label(label, `for` = inputId),
      tags$input(id = inputId, type = "text", value = value,...))
}

node_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(5,
           selectInput(ns("node_dist"),
                       label = paste0(id,"'s distribution"),
                       choices = c("binary", "continuous"))
    ),
    column(7,
           uiOutput(ns("ui_placeholder"))

    )
  )
}

node_server <- function(input, output, session) {
  return_value <- reactive({list(distribution = input$node_dist,c(prevalence = input$prev,
                              mean = input$mean,
                              std_dev = input$stdev
                              ))})
  ns <- session$ns
  node_name = sub("-", "", session$ns(""))
  output$ui_placeholder <- renderUI({
    type <- req(input$node_dist)
    if(type == "binary") {
      numericInput(ns("prev"), paste0(node_name,"'s prevalence (between 0-1):"), 0.5)

    } else if (type == "continuous") {
      tagList(numericInput(ns("mean"), paste0("Mean of ", node_name,":"),0),
              numericInput(ns("stdev"), paste0("Std dev of ", node_name,":"), 1),
              )


    }
  })

  ## if we later want to do some more sophisticated logic
  ## we can add reactives to this list
  # list(return_value = reactive({input$node_dist}), try = reactive({input$inner_element}))
  return_value
}

ui <- fluidPage(
  titlePanel(h1("Bias Simulator", align ="center")),
  fluidRow(column(12, mainPanel(p("format"),
                                textInput("dag_text", label = h3("Input DAG formula"),
                                          value = "~exposure|confounder + outcome|confounder*exposure + collider|exposure*outcome",
                                          width = "100%"),
                                verbatimTextOutput("bn_info"),
                                width = 12)
                  )
           ),
  sidebarLayout(position = "right",

                fluidRow(column(6, sidebarPanel(h4("---Set Distributions---", align = "center"),
                                                uiOutput("node_box"),
                                                div(id="placeholder"),
                                                h4("---Set Beta Values---", align = "center"),
                                                uiOutput("beta_box"),
                                                h4("---Set Analysis---", align = "center"),
                                                uiOutput("exp_out"),
                                                uiOutput("conf_sel"),
                                                width = 12
                )),
                ),
                column(6, align = "center", mainPanel(h4("Dag plot"),
                                    grVizOutput("value"),


                                    width = 12
                )),



  ),
  fluidRow(column(12, mainPanel(verbatimTextOutput("out"))))

)

# Define server logic ----
server <- function(input, output, session) {

  handler = reactiveVal(list())

  output$value = renderGrViz({ dag_ui(paste(input$dag_text)) })


  observeEvent(input$dag_text, {



    #beta interface
    beta_id = get_arrows(input$dag_text)
    output$beta_box = renderUI(
      map(beta_id, ~textInput3(.x, paste0( .x, ":"), value = 0, class = "input-small"))

    )

    #node distribution interface
    x <- get_nodes(input$dag_text)
    ui_num = length(x)
    new_id = unlist(lapply(c(1:ui_num), function(y) paste(x[y])))

    output$node_box = renderUI({
      map(new_id, ~node_ui(.x))
    })

    handler_list <- isolate(handler())

    new_handler = lapply(1:ui_num, function(y)
      callModule(node_server, new_id[y])
    )
    handler_list <-unlist(new_handler)
    names(handler_list) <- new_id
    handler(handler_list)


    #maybe set exposure first and render different ui for outcome so exp and out cant be set to same node
    output$exp_out = renderUI(


      fluidRow(column(6,
                      selectInput("exp", "Exposure", choices = x)
                      ),
               column(6,
                      selectInput("out", "Outcome", choices = x)
                      )
               )

    )

    output$conf_sel = renderUI(

      fluidRow(column(6,
                      checkboxGroupInput("conf", "Confounders", x[!x==input$exp & !x==input$out])
      ),
      column(6,
             selectInput("sel", "Selection Node", choices = handler_df(handler)[!handler_df(handler)==input$exp & !handler_df(handler)==input$out])
      )
      )

    )

    output$bn_info = renderPrint({
      # test_df = data.frame(unlist(lapply(handler(), function(handle) {
      #   handle()[["distribution"]]
      # })), holder = 1)
      # rownames(test_df[test_df[,1]=="binary",])

      unlist(lapply(x, function(y) {paste0("dag = setNode(dag, ", y, ", nodeType = ",get_dist(y, handler),")")}))

      # as.character(lapply(handler(), function(handle) {
      #   handle()[["distribution"]]})[["exposure"]])

      # test_dag = function(ate, exp_p, out_p){
      #   oop = out_p
      #   oppy = exp_p
      #   dag = HydeNetwork(~age
      #                     + sex
      #                     + depression| sex * age
      #                     + diabetes|sex *  age
      #                     + in_biovu| sex * age * depression * diabetes)
      #
      #   plot(dag)
      #
      #   dag = setNode(dag, age, nodeType = "dnorm", mu =1, tau = 1 )
      #   dag = setNode(dag, sex, nodeType = "dbern", prob = 0.505)
      #   dag = setNode(dag, depression, nodeType = "dbern",
      #                 prob = paste0("ilogit(", 0.1, "* age +",
      #                               0.9, "* sex +",
      #                               set_p(0.16, 0.1 * 1
      #                                     + 0.505 * 0.9)
      #                               ,")")
      #   )
      #   dag = setNode(dag, diabetes, nodeType = "dbern",
      #                 prob = paste0("ilogit(", 0.05, "* age + ",
      #                               -0.094, "* sex +",
      #                               ate, "* depression +",
      #                               set_p(0.074, 0.05 * 1
      #                                     + 0.505 * -0.094
      #                                     + ate * 0.16)
      #                               ,")")
      #   )
      #   dag = setNode(dag, in_biovu, nodeType = "dbern",
      #                 prob = paste0("ilogit(", 0.1, "* age + ",
      #                               0.25, "* sex +",
      #                               0.23, "* depression +",
      #                               0.92, "* diabetes +",
      #                               set_p(0.1, 0.1 * 1
      #                                     + 0.25 * 0.505
      #                                     + 0.92 * 0.074
      #                                     + 0.23 * 0.16)
      #                               ,")")
      #   )
      # }


    })


  })


}

# Run the app ----
shinyApp(ui = ui, server = server)
