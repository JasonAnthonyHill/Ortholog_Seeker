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
      actionButton("selectall", label="Select/Deselect all"),
      checkboxGroupInput("species",
                         "Select a species:", 
                         choices = unique(orthologs$species),
                         selected = unique(orthologs$species)),
      textInput("gene", "Enter an Gene ID:", value = "AT2G27680"),
      actionButton("go", "Update"),
      br(),
      downloadButton("download_filtered", "Download Filtered Ortholog List"),
      downloadButton("download_fasta", "Download Fasta Sequences"),
      downloadButton("download_fasta_aln", "Download Aligned Fasta Sequences"),
      downloadButton("download_newick", "Download Newick Tree"),
      downloadButton("download_gff", "Download spruce GFF3")
    ),
    mainPanel(tabsetPanel(
      tabPanel("Ortholog count",
               h3("Ortholog count per species"),
               plotOutput("barstack") %>% withSpinner()
      ),
      tabPanel("Gene Tree",
               h3("Resolved gene tree"),
               plotOutput("tree_plot") %>% withSpinner()
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
      Biostrings::writeXStringSet(fasta(), filepath = file)
    }
  )
  
  # Align fasta file
  aln_fasta <- eventReactive(input$go, {
    fasta() %>% 
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
    ape::read.tree(paste0(datapath, "Resolved_Gene_Trees/", 
                          filtered_data()$Ortholog_Group %>% unique(), 
                          "_tree.txt"))
  })

  #Filter Newick Tree
  tree_filtered <- eventReactive(input$go, {
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
  output$tree_plot <- renderPlot(
    width = 1000,
    height = function() 20 * (tree_tipcount() + 1),
    {
    treefile <- tree_filtered()
    # treefile <- tree()
    ggtree::ggtree(treefile) + 
      ggtree::geom_tiplab() +
      ggtree::xlim(0,2) +
      ggplot2::labs(title = filtered_data()$Ortholog_Group)
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


