# Load required libraries
library(shiny)
library(shinyjs)
library(tidyverse)
library(Biostrings)
library(ape)
library(ggtree)
library(muscle)
library(rtracklayer)


# Import Data
orthologs <- read_tsv("Phylogenetic_Hierarchical_Orthogroups/N0.tsv") %>%
  gather("species","ortholog",4:27) %>% 
  separate_longer_delim(ortholog, delim = ", ") %>% 
  drop_na() %>% 
  separate_wider_delim(ortholog, 
                       cols_remove = F,
                       delim = " ", 
                       names = c("transcript", "gene"),
                       too_few = "align_start",
                       too_many = "drop") %>% 
  mutate(gene = if_else(grepl("gene=",gene), gene, transcript)) %>% 
  mutate(gene = str_replace_all(gene, c("gene=" = ""))) %>% 
  unite(seqid, species, gene, transcript, remove = F) %>% 
  select(-HOG, -`Gene Tree Parent Clade`) %>% 
  rename(Fasta_header = ortholog, Ortholog_Group = OG)

species_tree <- read.tree(file = "species_list_v4.nwk")
spruce_gff <- import("Picab02_codingAll.gff3.gz")

# Define UI for Shiny App
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Ortholog Seeker"),
  
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("species",
                         "Select a species:", 
                         choices = unique(orthologs$species),
                         selected = unique(orthologs$species)),
      textInput("gene", "Enter an Gene ID:", value = "AT2G27680"),
      actionButton("go", "Update"),
      br(),
      downloadButton("download_filtered", "Download Filtered Ortholog List"),
      downloadButton("download_fasta", "Download Aligned Fasta Sequences"),
      downloadButton("download_newick", "Download Newick Tree"),
      downloadButton("download_gff", "Download spruce GFF3")
    ),
    mainPanel(tabsetPanel(
      tabPanel("Ortholog count",
               h3("Ortholog count per species"),
               plotOutput("barstack")
      ),
      tabPanel("Gene Tree",
               h3("Resolved gene tree"),
               plotOutput("tree_plot")
      ),
      tabPanel("Ortholog table",
               h3("Ortholog Table"),
               tableOutput("table")
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
  
  # Filter the Ortholog Data
  filtered_data <- eventReactive(input$go, {
    orthologs %>%
      group_by(Ortholog_Group) %>% 
      filter(any(toupper(gene) == toupper(input$gene)), species %in% input$species) %>% 
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
      write_tsv(filtered_data(), file)
    }
  )
  
  # Load Fasta Sequences
  fasta <- eventReactive(input$go, {
    readAAStringSet(filepath = paste0("Orthogroup_Sequences/", filtered_data()$Ortholog_Group %>% unique(), ".fa.gz")) %>% 
      .[grepl(str_c(filtered_data() %>% 
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
  
  # Align fasta file
  aln_fasta <- eventReactive(input$go, {
    fasta() %>% 
      muscle() %>% as("AAStringSet")
  })

  # Download Filtered Fasta Sequences
  output$download_fasta <- downloadHandler(
    filename = function() {
      paste0("Ortholog_Sequences_aln_", filtered_data()$Ortholog_Group, ".fasta")
    },
    content = function(file) {
      writeXStringSet(aln_fasta(), filepath = file)
    }
  )

  # Load Newick Tree
  tree <- eventReactive(input$go, {
    read.tree(paste0("Resolved_Gene_Trees/", filtered_data()$Ortholog_Group %>% unique(), "_tree.txt"))
  })

  #Filter Newick Tree
  tree_filtered <- eventReactive(input$go, {
    tree <- tree()
    keep.tip(tree, grep(str_c(input$species,
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
    height = function() 20 * tree_tipcount(),
    {
    treefile <- tree_filtered()
    # treefile <- tree()
    ggtree(treefile) + 
      geom_tiplab() +
      xlim(0,2) +
      labs(title = filtered_data()$Ortholog_Group)
  })
  
  
  # Download Filtered Newick Tree
  output$download_newick <- downloadHandler(
    filename = function() {
      paste0("Ortholog_Tree_", filtered_data()$Ortholog_Group, ".nwk")
    },
    content = function(file) {
      write.tree(tree_filtered(), file = file)
      # write.tree(tree(), file = file)
    }
  )
  
  # Plot Barstack
  output$barstack <- renderPlot(filtered_data() %>%
                                  group_by(Ortholog_Group, species) %>%
                                  tally() %>%
                                  ggplot() +
                                  geom_col(aes(x = species, y = n, fill = n)) +
                                  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                                        legend.position = "none") +
                                  labs(y = "Number of orthologs", 
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
      subsetByOverlaps(spruce_gff, spruce_gff_subset) %>% 
        export.gff3(file)
    }
  )




}

# Run Shiny App
shinyApp(ui = ui, server = server)

