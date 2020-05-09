library(shiny)
library(DT)
library(shinythemes)


options(shiny.maxRequestSize=70*1024^2)

# Define UI for application that draws a histogram
shinyUI(
  navbarPage("GCNA Hub",theme=shinytheme("yeti"),
             
    ###gene query panel###         
    tabPanel(
      "Gene Query",
      sidebarPanel(
        fileInput("file","File Input (.csv only):",
                  multiple = FALSE,
                  accept = c("text/csv",".csv")),
          
        radioButtons("disp","Display",
                       choices = c(Head="Head",All="All"),
                       selected = "Head"),
          
        hr(),
          
        textInput("Target_Name",label = "Module Query","Input Gene to Query"),
        actionButton("run","Query", class = "btn-primary"),
          
        hr(),
          
        downloadButton("download","Download Result",class = "btn-primary")
        ),
      
    mainPanel(
      
      tabsetPanel(
        
        tabPanel(
          "Data Preview",
          h4("RawData Table"),
          #tableOutput("contents")
          DT::dataTableOutput("contents")
        ),
        
        tabPanel(
          "Query Result",
          h4("Gene in the Same Module"),
          #tableOutput("result"))
          DT::dataTableOutput("result")
        ),
        tabPanel(
          "Module Tree Map",
          h4("Tree Map Based on Module"),
          plotOutput("Treemap")
        )
      )
    )),
    ###WGCNA panel###
    tabPanel(
      "WGCNA",
      sidebarPanel(
        fileInput("rdata","File Input (.Rdata only):",
                  multiple = FALSE,
                  accept = c(".Rdata")),
        actionButton("run1","Run", class = "btn-primary"),
        
        hr(),
        
        sliderInput("slider", "R Square Selection:", 0, 1, 10),
        
        hr(),
        
        textInput("Target_Name",label = "Sotfthreshold Selection","Input softthreshold"),
        actionButton("run2","Run WGCNA", class = "btn-primary")),
      
      mainPanel(
        
        tabsetPanel(
          
          tabPanel(
            "Outlier detection",
            h4("Sample clustering to detect outliers"),
            plotOutput("outlier")
          ),
          
          tabPanel(
            "Softthreshold Selection",
            h4("Pick SoftThreshold"),
            plotOutput("sft")
          ),
          
          tabPanel(
            "Module identification",
            h4("Module indentification"),
            plotOutput("mod")
          ),
          
          tabPanel(
            "Module Correlation",
            h4("Eigengene adjacency heatmap"),
            plotOutput("heatmap")
          ),
          
          tabPanel(
            "Gene Frequency in Modules",
            h4("Gene Frequency in Modules"),
            plotOutput("genefre")
          )
        )
      )
    )
  )
)
