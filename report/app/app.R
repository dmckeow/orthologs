library(shiny)
library(ggtree)
library(ggplot2)
library(dplyr)

print(getwd())
options(shiny.trace = TRUE)
options(shiny.fullstacktrace = TRUE)

# Load data from other scripts
load("../data/og_data.RData")

SUPERGROUP_COLS <- c(
  Choanoflagellata = "#A6CEE3",
  Cnidaria = "#1F78B4",
  Ctenophora = "#B2DF8A",
  Deuterostomia = "#33A02C",
  Filasterea = "#FB9A99",
  Placozoa = "#E31A1C",
  Porifera = "#FF7F00",
  Protostomia = "#CAB2D6",
  Teretosporea = "#6A3D9A"
)

OG_CALLER_COLS <- c(
  orthofinder = "#D95F02",
  broccoli = "#1B9E77"
)

# Define UI
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      #tree_plot img {
        width: 100%;
        height: auto;
      }
    "))
  ),
  titlePanel("Species tree"),
  
  fluidRow(
    column(2,  # Reduce the width of the sidebar
      selectInput("heatmap_df", "Select Heatmap Data Frame:",
                  choices = c("Gene Class" = "hm_total_genes_in_id_by_geneclass", 
                              "Orthogroup" = "hm_total_genes_in_id_by_og")),
      
      uiOutput("column_select")
    ),
    column(12,  # Increase the width of the main panel
      plotOutput("tree_plot", height = "800px")  # Increase the height of the plot
    )
  )
)

# Define Server Logic
server <- function(input, output, session) {
  
  # Reactive expression to get the selected heatmap data
  selected_heatmap_data <- reactive({
    req(input$heatmap_df)  # Ensure the heatmap data frame is selected
    if(input$heatmap_df == "hm_total_genes_in_id_by_geneclass") {
      return(tree_meta$hm_total_genes_in_id_by_geneclass)
    } else {
      return(tree_meta$hm_total_genes_in_id_by_og)
    }
  })
  
  # Render the selectInput for heatmap columns dynamically based on the selected data frame
  output$column_select <- renderUI({
    data <- selected_heatmap_data()  # Get the current heatmap data
    selectInput("heatmap_columns", "Select Columns for Heatmap:",
                choices = colnames(data), 
                selected = colnames(data)[1], 
                multiple = TRUE,
                width = "100%")
  })
  
  # Reactive expression for filtered heatmap data
  filtered_heatmap_data <- reactive({
    req(input$heatmap_columns)
    data <- selected_heatmap_data()
    data[, input$heatmap_columns, drop = FALSE]
  })
  
  # Reactive plot rendering
  output$tree_plot <- renderPlot({
    req(filtered_heatmap_data())
    
    # Plot the tree with ggplot
    sp_tree <- ggtree(tree, aes(color = supergroup),
                      linewidth = 1.0, layout = "rectangular") %<+% tree_meta$tips +
      geom_tiplab() +
      scale_color_manual(values = SUPERGROUP_COLS)
    
    # Add heatmap with filtered data
    sp_tree <- gheatmap(sp_tree,
                        filtered_heatmap_data(),
                        offset = 8,
                        width = 2,
                        colnames_angle = 90,
                        colnames_offset_y = -1,
                        hjust = 1,
                        low = "white",
                        high = "blue",
                        color = "white")
    
    sp_tree + ggtree::vexpand(.1, -1) +
      coord_cartesian(clip = "off")
  })
}


# Run the application 
shinyApp(ui = ui, server = server)
