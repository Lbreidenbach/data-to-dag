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

set_p = function(p,model){
  p2 = log(p/(1-p))
  b0 = p2-model
  return(b0)
}

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
  return_value <- reactive({list(distribution = input$node_dist, prevalence = input$prev,
                              mean = input$mean,
                              std_dev = input$stdev
                              )})
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
      map(beta_id, ~textInput3(.x, paste0( .x), value = 0))

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

      ###CURRENT BUG, TEST DF FAILS WHEN ALL NODES CONTINUOUS


      test_df = data.frame(unlist(lapply(handler(), function(handle) {
        handle()[["distribution"]]
      })), holder = 1)

      test_df$prevalence = rep(0.5, length(x))
      test_df$mean = rep(0, length(x))
      test_df$std_dev = rep(1, length(x))

      rownames(test_df) = x

      # test_df[,3] = data.frame(unlist(lapply(handler(), function(handle) {
      #   handle()[["prevalence"]]
      # })))

      if("binary" %in% test_df[[1]]){
        prev_vals = test_df[,3] = unlist(lapply(handler(), function(handle) {
          handle()[["prevalence"]]
        }))
        test_df[names(prev_vals), 3] = prev_vals

      }

      if("continuous" %in% test_df[[1]]){
        mean_vals = unlist(lapply(handler(), function(handle) {
          handle()[["mean"]]
        }))
        test_df[names(mean_vals), 4] = mean_vals

        std_vals = unlist(lapply(handler(), function(handle) {
          handle()[["std_dev"]]
        }))
        test_df[names(std_vals), 5] = std_vals
      }



      colnames(test_df) = c("distribution", "holder", "prevalence", "mean", "std_dev")
      beta_df = data.frame(parent = unlist(lapply(beta_id, function(x){strsplit(x, " -> ")[[1]][1]})),
                           child = unlist(lapply(beta_id, function(x){strsplit(x, " -> ")[[1]][2]})),
                           beta_vals = unlist(lapply(beta_id, function(x) input[[paste(x)]]))
                           )



      bi_nodes_df = test_df[test_df$distribution=="binary",]
      cont_nodes_df = test_df[test_df$distribution=="continuous",]



      parent_list = lapply(c(1:length(rownames(test_df))), function(y) beta_df[beta_df[2]==rownames(test_df)[[y]], c(1,3)]) #child == y from the apply functions
      names(parent_list) = rownames(test_df)


      get_type =function(i,y){
        if(length(parent_list[[i]][[1]])==0){
          return(0)
        }
        if(test_df[parent_list[[i]][y,1],1]=="binary"){
          col_num = 3
        }else if(test_df[parent_list[[i]][y,1],1]=="continuous"){
          col_num = 4
        }
        as.numeric(test_df[parent_list[[i]][y,1],col_num])

      }
      get_lgm = function(i){
        if(length(parent_list[[i]][[1]]) == 0){
          return("")
        }
        out_lgm = paste(unlist(lapply(c(1:length(parent_list[[i]][[1]])), function(y){
          paste0(
            parent_list[[i]][y,2],
            '," * ',
            parent_list[[i]][y,1],
            ' + ",')
        })), collapse = " " )

        return(out_lgm)
      }

      get_linker = function(i){
        if(length(parent_list[[i]][[1]])==0){
          return(set_p(as.numeric(test_df[rownames(test_df)[[i]],3]), 0))
        }
        set_p(as.numeric(test_df[rownames(test_df)[[i]],3]),
              sum(unlist(lapply(c(1:length(parent_list[[i]][[1]])), function(y){
                as.numeric(parent_list[[i]][y,2]) *
                  get_type(i,y)
              }))))
      }

      if(nrow(bi_nodes_df)>0){
        bi_code = unlist(lapply(c(1:length(rownames(bi_nodes_df))), function(i){paste0('\ndag = setNode(dag, ',
                                                                                       rownames(bi_nodes_df)[[i]],
                                                                                       ', nodeType = "dbern", prob = paste0("ilogit(",',
                                                                                       paste(get_lgm(which(rownames(test_df) == rownames(bi_nodes_df)[[i]])), get_linker(i), collapse = " "),',")"))')}))
      }else{
        bi_code = NULL
      }

      if(nrow(cont_nodes_df)>0){
        cont_code = unlist(lapply(c(1:length(rownames(cont_nodes_df))), function(i)
        {paste0('\ndag = setNode(dag, ',
                rownames(cont_nodes_df)[[i]], ', nodeType = "dnorm", mu = paste0(',
                get_lgm(which(rownames(test_df) == rownames(cont_nodes_df)[[i]])),
                cont_nodes_df[i,4],
                "), tau = ",
                cont_nodes_df[i,5],
                ")")
        }))
      }else{
        cont_code = NULL
      }





      cat(c(bi_code, cont_code))







      ###Bugged lines



      # lapply(c(1:length(test_df[[1]])), function(i) paste(get_lgm(i), get_linker(i), collapse = " "))


      # lapply(c(1:length(parent_list[[4]][[1]])), function(y){
      #   as.numeric(parent_list[[4]][y,2]) *
      #     get_type(y,4)
      # })




      # paste0("illogit(",
      #        parent_list[[1]][1,2],
      #        " * ",
      #        parent_list[[1]][1,1],
      #        " + ",
      # set_p(as.numeric(test_df[names(parent_list[1]),3]),
      #       as.numeric(parent_list[[1]][1,2]) *
      #         get_type(1)), ")")





      #c(bi_code, cont_code)


      #strsplit(beta_id, " -> ")[[1]][1]


      # test_df[5] = data.frame(unlist(lapply(handler(), function(handle) {
      #   handle()[["sumstat"]][[3]]
      # })), holder = 1)

      # test_df$mean = data.frame(unlist(lapply(handler(), function(handle) {
      #   handle()[["sumstat"]][2]
      # })), holder = 1)


      #rownames(test_df[test_df[,1]=="binary",])

      #lapply(x, function(y) {paste0("dag = setNode(dag, ", y, ", nodeType = ",get_dist(y, handler),")")})

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
