# HackBio 2026
# This is a Shiny web application. 
#  Run locally with: shiny::runApp()

library(shiny)
library(dplyr)
#------------------------- Boilerplate functions ------------------------------

compute_gene_stats <- function( expr_mat, meta_df, target_ct){
  #identify cells inside and outside the selected cell type
  in_cells <- meta_df$cell_id[meta_df$cell_type == target_ct]
  out_cells <- meta_df$cell_id[meta_df$cell_type != target_ct] #subset expression matrix
  xin <- expr_mat[in_cells,, drop = FALSE]
  xout <- expr_mat[out_cells,, drop = FALSE] # Detection rates
  det_in <- colMeans(xin>0)
  det_out <- colMeans(xout>0) # Mean expression
  mean_in <- colMeans(xin)
  mean_out <- colMeans(xout) # Specificity difference
  diff <- mean_in - mean_out
  
  return(data.frame(
    gene = colnames(expr_mat),
    det_in = det_in,
    det_out = det_out,
    mean_in = mean_in,
    mean_out = mean_out,
    diff = diff,
    stringsAsFactors = FALSE)
  )
}

 ##select marker gene
pick_marker_gene <- function(gene_stats_df){
  gene_stats_df <- gene_stats_df%>%
    arrange(desc(diff), desc(det_in))
  gene_stats_df$gene[1]
}

##scale marker expression 0 - 100
scale_0_100 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if(rng[1] == rng[2]) return(rep(50, length(x)))
  (x - rng[1]) / (rng[2] - rng[1]) * 100
}

## -------------------------- Data-set ---------------------------------------

expr_matrix <- read.csv("https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/sc_synthetic/expression_matrix.csv", header = T,
                        row.names = 1)
cell_metadata <- read.csv("https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/sc_synthetic/cell_metadata.csv", header = T)
umap_coordinates <- read.csv("https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/sc_synthetic/umap_coordinates.csv", header = T)

 # Make sure IDs match and are aligned
common_ids <- Reduce(
  intersect, 
  list(rownames(expr_matrix), cell_metadata$cell_id, umap_coordinates$cell_id)
)

 #reorder metadata and UMAP to match expr-matrx rows
cell_metadata <- cell_metadata[match(common_ids, cell_metadata$cell_id), ]
umap_coordinates <- umap_coordinates[match(common_ids, umap_coordinates$cell_id),]

#Merged dataset of Metadata and UMAP by cell id

plot_df <- merge(cell_metadata,
                 umap_coordinates,
                 by = "cell_id")

# ----------------------- UI ---------------------------------------

ui <- fluidPage(
  titlePanel("Simulated Cell Type Viewer"),
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId = "cell_type",
                  label = "Select Cell Type :",
                  choices = c("All",unique(cell_metadata$cell_type)),
                  selected = "All"),
      
      verbatimTextOutput(outputId = "marker_gene"),
      
      actionButton(inputId = "buttons",
                   label = "Generate a Plot",
                   icon = icon("image"))
    ),
    
      mainPanel(
      tableOutput(outputId = "per_gene_stats"),
      plotOutput(outputId = "scatter_plot")
    )
  )
)

# ----------------  server  --------------------------------------------
server <- function(input, output) {
  
  #compute per-gene stats
  stats_reactive <- reactive({
    req(input$cell_type)
    target_ct <- input$cell_type
    
    compute_gene_stats(
      expr_mat = expr_matrix,
      meta_df = cell_metadata,
      target_ct = target_ct)
  })
  
  #per gene stats table
  output$per_gene_stats <- renderTable({
    stats_df <- stats_reactive()
    head(stats_df[order(- stats_df$diff, - stats_df$det_in),], 20)
  })
  
  #choose marker gene
  output$marker_gene <- renderText({
  stats_df <- stats_reactive()
  top_gene <- pick_marker_gene(stats_df)
  diff_val <- stats_df$diff[stats_df$gene == top_gene]
  paste("Cell type :", input$cell_type , "| Marker :", top_gene, "| Diff =", 
        round(diff_val, 3))
})
  
  #  scatter plot 
data_plot <- eventReactive(input$buttons, {
  stats_df <- stats_reactive()
  marker_gene <- pick_marker_gene(stats_df) #get marker gene from stats table
  
  marker_expr <- expr_matrix[, marker_gene] #get per-cell expression 
  
  marker_expre_scaled <- scale_0_100(marker_expr) #scale 0 - 100
  plot_df$marker_expression <- marker_expre_scaled
  
  list(
    x = plot_df[,"UMAP_1"],
    y = plot_df[,"UMAP_2"],
    cell_type = input$cell_type,
    marker_gene = marker_gene,
    marker_expr = plot_df$marker_expression
  )
})


  output$scatter_plot <- renderPlot({
    pts <- data_plot()
    col_scale <- colorRampPalette(c("lightgrey","darkred"))(100)[
      cut(pts$marker_expr, breaks = 100, include.lowest = T)
    ]
    
    plot(x = pts$x,
         y = pts$y,
         col = col_scale,
         ylab = "UMAP_2",
         xlab = "UMAP_1",
         bg = as.factor(cell_metadata$cell_type),
         pch = 21,
         cex = 1.3,
         main = paste("UMAP: ", pts$marker_gene))
  })
}


# -------------------- Run the application ----------------------------
shinyApp(ui = ui, server = server)

