#!/usr/bin/env Rscript

# Install phyloseq package before beginning, if required.
# if (!require("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")
# BiocManager::install("phyloseq")

#### Load required libraries. Run each time a new session is opened.####
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

#### CREATING THE PHYLOSEQ OBJECT (March 1st and March 2nd, 2026 by Michaela) ####

#### Loading in the data files. ####
# Metadata
fish_metadata_object <- "fish_metadata.txt"
fish_metadata <- read_delim(fish_metadata_object, delim="\t")

# OTU Table
fish_otu_object <- "feature-table.txt"
fish_otu <- read_delim(file = fish_otu_object, delim="\t", skip=1)

# Taxonomy Table
fish_taxonomy_object <- "taxonomy.tsv"
fish_taxonomy <- read_delim(fish_taxonomy_object, delim="\t")

# Rooted Tree 
fish_tree_object <- "tree.nwk"
fish_tree <- read.tree(fish_tree_object)

#### Creating a separate dataset for each body site ####
# Gill 
gill_metadata <- filter(fish_metadata, sample_type == "gill")

# Hindgut 
hindgut_metadata <- filter(fish_metadata, sample_type == "hindgut")

# Midgut 
midgut_metadata <- filter(fish_metadata, sample_type == "midgut")

# Skin 
skin_metadata <- filter(fish_metadata, sample_type == "skin")

#### Formatting the data files ####
# Converting the OTU Table into a matrix
otu_matrix <- as.matrix(fish_otu[,-1])
rownames(otu_matrix) <- fish_otu$`#OTU ID`
fish_OTU_matrix <- otu_table(otu_matrix, taxa_are_rows = TRUE) 
class(fish_OTU_matrix)

# Formatting the gill metadata 
gill_metadata_dataframe <- as.data.frame(gill_metadata[,-1])
rownames(gill_metadata_dataframe)<- gill_metadata$'#SampleID'
gill_metadata_phylo <- sample_data(gill_metadata_dataframe)
class(gill_metadata_phylo)

# Formatting the hindgut metadata 
hindgut_metadata_dataframe <- as.data.frame(hindgut_metadata[,-1])
rownames(hindgut_metadata_dataframe)<- hindgut_metadata$'#SampleID'
hindgut_metadata_phylo <- sample_data(hindgut_metadata_dataframe)
class(hindgut_metadata_phylo)

# Formatting the midgut metadata 
midgut_metadata_dataframe <- as.data.frame(midgut_metadata[,-1])
rownames(midgut_metadata_dataframe)<- midgut_metadata$'#SampleID'
midgut_metadata_phylo <- sample_data(midgut_metadata_dataframe)
class(midgut_metadata_phylo)

# Formatting the skin metadata 
skin_metadata_dataframe <- as.data.frame(skin_metadata[,-1])
rownames(skin_metadata_dataframe)<- skin_metadata$'#SampleID'
skin_metadata_phylo <- sample_data(skin_metadata_dataframe)
class(skin_metadata_phylo)

# Formatting the taxonomy file
fish_taxonomy_matrix <- fish_taxonomy %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() 
fish_taxonomy_matrix <- fish_taxonomy_matrix[,-1]
rownames(fish_taxonomy_matrix) <- fish_taxonomy$`Feature ID`
fish_taxa_table <- tax_table(fish_taxonomy_matrix)
class(fish_taxa_table)

#### Making phyloseq objects for each body site #### 
# Gill
gill_phyloseq_object <- phyloseq(fish_OTU_matrix, gill_metadata_phylo, fish_taxa_table, fish_tree)

# Hindgut
hindgut_phyloseq_object <- phyloseq(fish_OTU_matrix, hindgut_metadata_phylo, fish_taxa_table, fish_tree)

# Midgut 
midgut_phyloseq_object <- phyloseq(fish_OTU_matrix, midgut_metadata_phylo, fish_taxa_table, fish_tree)

# Skin 
skin_phyloseq_object <- phyloseq(fish_OTU_matrix, skin_metadata_phylo, fish_taxa_table, fish_tree)

#### Visualizing componenents of the phyloseq objects #### 
# Gill 
otu_table(gill_phyloseq_object)
sample_data(gill_phyloseq_object)
tax_table(gill_phyloseq_object)
phy_tree(gill_phyloseq_object)

# Hindgut
otu_table(hindgut_phyloseq_object)
sample_data(hindgut_phyloseq_object)
tax_table(hindgut_phyloseq_object)
phy_tree(hindgut_phyloseq_object)

# Midgut 
otu_table(midgut_phyloseq_object)
sample_data(midgut_phyloseq_object)
tax_table(midgut_phyloseq_object)
phy_tree(midgut_phyloseq_object)

# Skin 
otu_table(skin_phyloseq_object)
sample_data(skin_phyloseq_object)
tax_table(skin_phyloseq_object)
phy_tree(skin_phyloseq_object)



