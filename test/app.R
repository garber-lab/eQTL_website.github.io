#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
load(file="./temp_test.RData")

# Define UI for application
ui <- fluidPage(
  
  # Application title
  titlePanel("Plotting on-demand"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    sidebarPanel(
      # Text input
      textInput("input_gene", "Enter gene name:", ""),
      textInput("input_snp", "Enter SNP ID:",""),
      
      # Plot button
      actionButton("plot_button", "Plot")
    ),
    
    # Output: Plot
    mainPanel(
      plotOutput("text_plot")
    )
  )
)

# Define server logic
server <- function(input, output) {
  #print("loading environment")
  #dataEnv <- new.env()
  #load(file="temp_test.RData", envir=dataEnv)
  #print("loaded environment")
  
  # Render plot when plot button is clicked
  observeEvent(input$plot_button, {
    print("plot button clicked")
    gene <- input$input_gene
    snp <- input$input_snp
    output$text_plot <- renderPlot({
      print("Calling make_reQTL_plot_CPM_3cts")
      plot_gene_expression(table, gene_of_interest = gene_of_interest)
    })
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
