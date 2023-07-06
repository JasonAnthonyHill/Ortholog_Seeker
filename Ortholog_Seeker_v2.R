# Load required libraries
library(shiny)
library(shinyjs)
library(tidyverse)
library(BiocManager)
options(repos = BiocManager::repositories())
library(Biostrings)
library(ape)
library(ggtree)
library(ggplot2)
library(muscle)
library(magrittr)
library(arrow)
library(shinycssloaders)
library(tidytree)
library(fuzzyjoin)

#
local <- F
datapath <- if_else(local == T, "./", "https://github.com/JasonAnthonyHill/Ortholog_Seeker/raw/main/")

# Import Data
orthologs <- read_parquet(paste0(datapath, "orthologs.parquet"))
gene_list <- toupper(orthologs$gene)
if(file.exists("annotation.RDS") == FALSE) {
   download.file("https://github.com/JasonAnthonyHill/Ortholog_Seeker/raw/main/annotation.RDS", "annotation.RDS", mode = "wb")
}
spruce_gff <- readRDS("annotation.RDS")
# species_tree <- read.tree(file = "species_list_v4.nwk")


# Define UI for Shiny App
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Ortholog Seeker"),
  sidebarLayout(
    sidebarPanel(
      width = 2,
      actionButton("selectall", label="Select/Deselect all"),
      checkboxGroupInput("species",
                         "Select a species:", 
                         choices = unique(orthologs$species),
                         selected = unique(orthologs$species))
    ),
    mainPanel(tabsetPanel(
      tabPanel("Gene Tree",
               h3("Ortholog tree and sequence download"),
               br(),
               h4("Useage instructions"),
               p("Select species of interest on the sidebar and enter a gene name from one of those species in the box below, then hit GO!"),
               p("The genes present in the ortholog group will be displayed below in their inferred phylogenetic relationship, with the search gene highlighted in",
                 span("RED", style = "color:red"), "."),
               p("Use your mouse to select a subset of genes in the tree, or the entire tree. Then use the download buttons to retrieve the fasta sequences of the selected genes in either aligned or unaligned format. Also available for download are the full phylogenetic tree of the orthogroup, the GFF3 file of any identified orthologs in the spruce genome, and an extened table with more details of each ortholog in the group."),
               fluidRow(
                 column(2, textInput("gene", "Enter an Gene ID:", value = "AT2G27680")),
                 column(2, style = "margin-top: 25px;", actionButton("go", "GO!"))
               ),
               tableOutput("selectedOrthos"),
               downloadButton("download_fasta", "Download Fasta Sequences"),
               downloadButton("download_fasta_aln", "Download Aligned Fasta Sequences"),
               downloadButton("download_newick", "Download Newick Tree"),
               downloadButton("download_gff", "Download spruce GFF3"),
               downloadButton("download_filtered", "Download Filtered Ortholog List"),
               plotOutput("tree_plot", brush = "plot_brush") %>% withSpinner()
      ),
      tabPanel("Ortholog count",
               h3("Ortholog count per species"),
               plotOutput("barstack") %>% withSpinner()
      ),
      tabPanel("Ortholog table",
               h3("Ortholog Table"),
               tableOutput("table") %>% withSpinner()
      )
      )
    )
  )
)

# Define server for Shiny App
server <- function(input, output) {
  
  # Initialize with go press
  o <- observe({
    shinyjs::click("go")
    o$destroy() # destroy observer as it has no use after initial button click
  })
  
  # Select or deselect all
  observe({
    if (input$selectall > 0) {
      if (input$selectall %% 2 == 0){
        updateCheckboxGroupInput(inputId="species",
                                 choices = unique(orthologs$species),
                                 selected = unique(orthologs$species))
        
      }
      else {
        updateCheckboxGroupInput(inputId="species",
                                 choices = unique(orthologs$species),
                                 selected = c())
        
      }}
  })
  
  # Filter the Ortholog Data
  filtered_data <- eventReactive(input$go, {
    validate(
      need(toupper(trimws(input$gene)) %in% gene_list,
           "Gene ID not found in dataset")
    )
    orthologs %>%
      group_by(Ortholog_Group) %>% 
      filter(any(toupper(gene) == toupper(trimws(input$gene))), species %in% input$species) %>% 
      ungroup()
  })
  
  # Output Table of Ortholog Information
  output$table <- renderTable({
    filtered_data()
  })
  
  # Download Filtered Ortholog Information
  output$download_filtered <- downloadHandler(
    filename = function() {
      paste0("Ortholog_List_", filtered_data()$Ortholog_Group, ".tsv")
    },
    content = function(file) {
      readr::write_tsv(filtered_data(), file)
    }
  )
  
  # Load Fasta Sequences
  fasta <- eventReactive(input$go, {
    Biostrings::readAAStringSet(filepath = paste0(datapath, "Orthogroup_Sequences/", filtered_data()$Ortholog_Group %>% unique(), ".fa.gz")) %>% 
      .[grepl(stringr::str_c(filtered_data() %>% 
                      filter(species %in% input$species) %>% 
                      pull(gene),
                    collapse = "|"), .@ranges@NAMES)]
  })
  
  # Rename Fasta Sequences
  # names(fasta) <- eventReactive(input$go, {
  #   filtered_data() %>% 
  #     filter(ortholog %in% names(fasta)) %>% 
  #     pull(seqid)
  # })
  
  # Download Filtered Fasta Sequences
  output$download_fasta <- downloadHandler(
    filename = function() {
      paste0("Ortholog_Sequences_", filtered_data()$Ortholog_Group, ".fasta")
    },
    content = function(file) {
      Biostrings::writeXStringSet(fasta() %>% 
                                    .[grepl(stringr::str_c(brushedPoints(make_tree()$data, input$plot_brush) %>% 
                                                             select(species, gene, seqid) %>% 
                                                             drop_na() %>% 
                                                             pull(gene),
                                                           collapse = "|"), .@ranges@NAMES)], filepath = file)
    }
  )
  
  # Align fasta file
  aln_fasta <- eventReactive(input$go, {
    fasta() %>% 
      .[grepl(stringr::str_c(brushedPoints(make_tree()$data, input$plot_brush) %>% 
                               select(species, gene, seqid) %>% 
                               drop_na() %>% 
                               pull(gene),
                             collapse = "|"), .@ranges@NAMES)] %>% 
      muscle::muscle() %>% as("AAStringSet")
  })

  # Download Filtered and Aligned Fasta Sequences
  output$download_fasta_aln <- downloadHandler(
    filename = function() {
      paste0("Ortholog_Sequences_aln_", filtered_data()$Ortholog_Group, ".fasta")
    },
    content = function(file) {
      Biostrings::writeXStringSet(aln_fasta(), filepath = file)
    }
  )

  # Load Newick Tree
  tree <- eventReactive(input$go, {
    validate(
      need(toupper(trimws(input$gene)) %in% gene_list,
           "Gene ID not found in dataset") 
    )
    ape::read.tree(paste0(datapath, "Resolved_Gene_Trees/", 
                          filtered_data()$Ortholog_Group %>% unique(), 
                          "_tree.txt"))
  })

  #Filter Newick Tree
  tree_filtered <- eventReactive(input$go, {
    validate(
      need(toupper(trimws(input$gene)) %in% gene_list,
           "Gene ID not found in dataset")
    )
    tree <- tree()
    tree$tip.label <- gsub("gene_gene","gene",tree$tip.label)
    tree$tip.label <- gsub("gene","_gene",tree$tip.label)
    tree$tip.label <- gsub("seq_id","_seq_id",tree$tip.label)
    ape::keep.tip(tree, grep(stringr::str_c(input$species,
                          collapse = "|"),
                        tree$tip.label)
    )
  }
  )
  
  tree_tipcount <- eventReactive(input$go, {
    tree_filtered()$tip.label %>% 
      length()
  }
  )
  
  # Output filtered tree
  make_tree <- eventReactive(input$go, {
    
    treemetadata <- tibble(tiplabel = tree_filtered()$tip.label) %>% 
      regex_left_join(filtered_data(),
                      by = c("tiplabel" = "gene")) %>% 
      mutate(selected_gene_col = if_else(toupper(gene) == toupper(trimws(input$gene)), "red", "black"))
    
    p <- ggtree(tree_filtered()) %<+% treemetadata
    p +
      geom_tiplab(aes(label = paste(species, gene), color = selected_gene_col)) +
      scale_color_identity() +
      xlim(0,2) +
      labs(title = filtered_data()$Ortholog_Group)
  })
  
  output$tree_plot <- renderPlot(
    width = 1000,
    height = function() 20 * (tree_tipcount() + 1),
    {
      make_tree()
  })
  
  # Catch the selected tips in box brush
  output$selectedOrthos <- renderTable({
    brushedPoints(make_tree()$data, input$plot_brush) %>% 
      select(species, gene, seqid) %>% 
      drop_na()
  })
  
  
  # Download Filtered Newick Tree
  output$download_newick <- downloadHandler(
    filename = function() {
      paste0("Ortholog_Tree_", filtered_data()$Ortholog_Group, ".nwk")
    },
    content = function(file) {
      ape::write.tree(tree_filtered(), file = file)
      # write.tree(tree(), file = file)
    }
  )
  
  # Plot Barstack
  output$barstack <- renderPlot(filtered_data() %>%
                                  group_by(Ortholog_Group, species) %>%
                                  tally() %>%
                                  ggplot2::ggplot() +
                                  ggplot2::geom_col(aes(x = species, y = n, fill = n)) +
                                  ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                                        legend.position = "none") +
                                  ggplot2::labs(y = "Number of orthologs", 
                                       x = "species",
                                       title = filtered_data()$Ortholog_Group)
  )
  
  # Filter and download GFF
  output$download_gff <- downloadHandler(
    filename = function() {
      paste0("Picea_abies_Orthologs_", filtered_data()$Ortholog_Group, ".gff3")
    },
    content = function(file) {
      spruce_genes <- filtered_data() %>%
        filter(species == "Picea_abies") %>%
        pull(gene)
      spruce_gff_subset <- subset(spruce_gff, ID %in% spruce_genes)
      IRanges::subsetByOverlaps(spruce_gff, spruce_gff_subset) %>%
        rtracklayer::export.gff3(file)
    }
  )




}

# Run Shiny App
shinyApp(ui = ui, server = server)


