orthologs %>%
  group_by(OG) %>% 
  filter(any(ortholog == "PAC_26995852") %>%
      group_by(OG) %>% 
      filter(any(ortholog == input$ortholog), species == input$species) %>% 
      ungroup()) %>% 
  ungroup()

orthologs %>%
  group_by(OG) %>% 
  filter(any(ortholog == "rna-gnl|WGS_JAHRHJ|evm2.model.Chr12.1479")) %>% 
  ungroup()

orthologs <- read_tsv("C:/Users/hillg/Dropbox/Research/Projects/conifer/gene_fams/03_orthology_search/orthofinder/fastree_allspp/Results_Jan19/Phylogenetic_Hierarchical_Orthogroups/N0.tsv.gz") %>%
  gather("species","ortholog",4:26) %>% 
  separate_longer_delim(ortholog, delim = ", ") %>% 
  drop_na()

write_tsv(orthologs, "C:/Users/hillg/Dropbox/Research/Projects/conifer/OrthologGroup_ShinyApp/Ortholog_Seeker/orthologs.tsv")

testfasta <- readAAStringSet("fastas/OG0007253.fa")
testfasta$PAC_19649926
names(testfasta)

orthologs %>% 
  filter(ortholog %in% names(testfasta))

setdiff(names(testfasta), orthologs %>% pull(ortholog))
orthologs %>% filter(OG == "OG0007253") %>% view()

tx_gene <- read_delim("Ortholog_Seeker/trx_gene_spp_lookupTable.txt", delim = " ", col_names = c("species", "ortholog", "gene"))

fasta %>% 
  as.data.frame() %>%
  rownames_to_column(var = "name") %>% 
  filter(str_detect(name, str_c(c("nbis","PIC"), collapse = "|"))) %>%
  AAStringSet(x)
  
filtered_data <- orthologs %>% filter(OG == "OG0003567")
  readAAStringSet(filepath = paste0("Orthogroup_Sequences/", filtered_data$OG %>% unique(), ".fa.gz")) %>% 
    .[grepl(str_c(c("nbis","PIC"), collapse = "|"), .@ranges@NAMES)] %>% #names()
    writeXStringSet("test.fa")
  
orthologs <-  read_tsv("Phylogenetic_Hierarchical_Orthogroups/N0.tsv") %>%
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
    unite(seqid, species, gene, transcript, remove = F)
  
names(fasta) <- str_replace_all(names(fasta), setNames(orthologs$ortholog, orthologs$seqid))

names(fasta) <- orthologs %>% filter(ortholog %in% names(fasta)) %>% pull(seqid)

library(ape)
library(ggtree)
tree <- read.tree("Resolved_Gene_Trees/OG0012831_tree.txt")
tree$tip.label <- gsub("gene_gene","gene",tree$tip.label)
tree$tip.label <- gsub("gene","_gene",tree$tip.label)
tree$tip.label <- gsub("seq_id","_seq_id",tree$tip.label)
tree$tip.label <- gsub("seq_id","_seq_id",tree$tip.label, )
tree.tbl <- as_tibble(tree)
plot(tree)
grep("Pinus", tree$tip.label)
keep.tip(tree, grep("Pinus", tree$tip.label)) %>% plot()

filtered_data %>%
  group_by(OG, species) %>% 
  tally() %>% 
  ggplot() +
  geom_col(aes(x = species, y = n, fill = n)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "none") +
  labs(y = "Number of orthologs",
       x = "species")

tree %>% 
  ggtree() +
  geom_tiplab() +
  xlim(0,10)

##gff
library(rtracklayer)
spruce_gff <- import("Picab02_codingAll.gff3.gz")
spruce_genes <- filtered_data %>% 
  filter(species == "Picea_abies") %>% 
  pull(gene)
spruce_gff_subset <- subset(spruce_gff, ID %in% spruce_genes)
subsetByOverlaps(spruce_gff, spruce_gff_subset) %>% 
  export.gff3("../test.gff3")

test <- read_tsv("https://github.com/JasonAnthonyHill/Ortholog_Seeker/raw/main/Phylogenetic_Hierarchical_Orthogroups/N0.tsv.gz")

## Test loading times for Perquet vs. RDS
library(arrow)
start_time <- Sys.time()
# data <- readRDS("test.RDS")
# data <- read_parquet("test.parquet")
# spruce_gff <- rtracklayer::import.gff3("Picab02_codingAll.gff3.gz")
spruce_gff <- readRDS("testgff.RDS")
end_time <- Sys.time()
end_time - start_time
