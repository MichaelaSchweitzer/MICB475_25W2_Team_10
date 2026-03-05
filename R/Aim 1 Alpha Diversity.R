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


# Saving Objects (Can now just load these instead of running all of the above code)
save(skin_phyloseq_object, midgut_phyloseq_object, hindgut_phyloseq_object, gill_phyloseq_object, file = "all_phyloseq_objects.RData") 

load("all_phyloseq_objects.RData")


#Get skin alpha diversity metrics 
alpha_div_skin <- estimate_richness(skin_phyloseq_object, measures = "shannon")
alpha_div_hindgut <- estimate_richness(hindgut_phyloseq_object, measures = "shannon")
alpha_div_midgut <- estimate_richness(midgut_phyloseq_object, measures = "shannon")
alpha_div_gill <- estimate_richness(gill_phyloseq_object, measures = "shannon")

#Create a dataframe containing sampledata
skin_data <- sample_data(skin_phyloseq_object)
hindgut_data <- sample_data(hindgut_phyloseq_object)
midgut_data <- sample_data(midgut_phyloseq_object)
gill_data <- sample_data(gill_phyloseq_object)


#Add shannon scores to sample data 
skin_div <- cbind(alpha_div_skin, skin_data)
hindgut_div <- cbind(alpha_div_hindgut, hindgut_data)
gill_div <- cbind(alpha_div_gill, gill_data)
midgut_div <- cbind(alpha_div_midgut, midgut_data)

### Before we iterate through our phyloseq object, we have to first remove all "not applicable"'s, and convert numbers to numeric values (instead of characters)
# Replace "not applicable"'s to NA (to do this we need to remove date formatting first bc R gets confused if not)
skin_div_nodate <- select(skin_div, -"collection_date")
skin_div_nodate[skin_div_nodate == "not applicable"] <- NA

hindgut_div_nodate <- select(hindgut_div, -"collection_date")
hindgut_div_nodate[hindgut_div_nodate == "not applicable"] <- NA

midgut_div_nodate <- select(midgut_div, -"collection_date")
midgut_div_nodate[midgut_div_nodate == "not applicable"] <- NA

gill_div_nodate <- select(gill_div, -"collection_date")
gill_div_nodate[gill_div_nodate == "not applicable"] <- NA


# Make a vector of all the numeric variables
num_vars <- c("dist_to_dorsal_cm", "fl_cm", "gape_cm", "gi_cm", "host_body_mass_index",
             "host_height", "host_height_vs_max_height_tl", "mass_g", "month", "ratio_dorsal_to_tl", "ratio_gape_to_tl", 
             "ratio_gi_to_tl", "tl_cm")


# For loop to make all numeric variables numeric class
skin_div_num <- skin_div_nodate
for (n in num_vars) {
  skin_div_num[[n]] <- as.numeric(skin_div_nodate[[n]]) 
}

gill_div_num <- gill_div_nodate
for (n in num_vars) {
  gill_div_num[[n]] <- as.numeric(gill_div_nodate[[n]]) 
}

hindgut_div_num <- hindgut_div_nodate
for (n in num_vars) {
  hindgut_div_num[[n]] <- as.numeric(hindgut_div_nodate[[n]]) 
}

midgut_div_num <- midgut_div_nodate
for (n in num_vars) {
  midgut_div_num[[n]] <- as.numeric(midgut_div_nodate[[n]]) 
}

#Create a vector containing the names of our columns of interest
vars <- c("dist_to_dorsal_cm", "fl_cm", "gape_cm", "gi_cm", "host_body_mass_index",
          "host_height", "host_height_vs_max_height_tl", "mass_g", "method_gear",
          "month", "ratio_dorsal_to_tl", "ratio_gape_to_tl", "ratio_gi_to_tl", "swim_mode",
          "swim_performance", "tl_cm")

#Sanity check to prove 
for (i in 1:length(hindgut_div_num)) {
  if (class(hindgut_div_num[[i]]) == "numeric") {
    print("is numeric") }
  else {
    (print("not_numeric")) }
}


### For loop to calculate significance for alpha diversity metrics ###

#Create an empty dataframe to hold P values 
skin_results <- data.frame(variable = character(), p_value = numeric(), method = character())
hindgut_results <- data.frame(variable = character(), p_value = numeric(), method = character())
midgut_results <- data.frame(variable = character(), p_value = numeric(), method = character())
gill_results <- data.frame(variable = character(), p_value = numeric(), method = character())

#Iterate through each variable, calculating alpha diversity metrics
#Skin for loop
for (v in vars) {
  if (class(skin_div_num[[v]]) == "numeric") {
    p <- cor.test(skin_div_num$Shannon, 
                  skin_div_num[[v]],
                  method = "spearman",
                  use = "complete.obs")$p.value
    skin_results <- rbind(skin_results, data.frame(variable = v, p_value = p, method = "Spearman")) }
  else {
    p <- kruskal.test(as.formula(paste("Shannon ~", v)), data = skin_div_num)$p.value
    skin_results <- rbind(skin_results, data.frame(variable = v, p_value = p, method = "Kruskal-Wallis")) 
  }
}

#Hindgut for loop

for (v in vars) {
  if (class(hindgut_div_num[[v]]) == "numeric") {
    p <- cor.test(hindgut_div_num$Shannon, 
                  hindgut_div_num[[v]],
                  method = "spearman",
                  use = "complete.obs")$p.value
    hindgut_results <- rbind(hindgut_results, data.frame(variable = v, p_value = p, method = "Spearman")) }
  else {
    p <- kruskal.test(as.formula(paste("Shannon ~", v)), data = hindgut_div_num)$p.value
    hindgut_results <- rbind(hindgut_results, data.frame(variable = v, p_value = p, method = "Kruskal-Wallis")) 
  }
}

#Midgut for loop
for (v in vars) {
  if (class(midgut_div_num[[v]]) == "numeric") {
    p <- cor.test(midgut_div_num$Shannon, 
                  midgut_div_num[[v]],
                  method = "spearman",
                  use = "complete.obs")$p.value
    midgut_results <- rbind(midgut_results, data.frame(variable = v, p_value = p, method = "Spearman")) }
  else {
    p <- kruskal.test(as.formula(paste("Shannon ~", v)), data = midgut_div_num)$p.value
    midgut_results <- rbind(midgut_results, data.frame(variable = v, p_value = p, method = "Kruskal-Wallis")) 
  }
}

#Gill for loop
for (v in vars) {
  if (class(gill_div_num[[v]]) == "numeric") {
    p <- cor.test(gill_div_num$Shannon, 
                  gill_div_num[[v]],
                  method = "spearman",
                  use = "complete.obs")$p.value
    gill_results <- rbind(gill_results, data.frame(variable = v, p_value = p, method = "Spearman")) }
  else {
    p <- kruskal.test(as.formula(paste("Shannon ~", v)), data = gill_div_num)$p.value
    gill_results <- rbind(gill_results, data.frame(variable = v, p_value = p, method = "Kruskal-Wallis")) 
  }
}




# continuing alpha diversity analysis (march 4th, madi)

# save alpha diversity results and loop results as object 
save(skin_div_num, midgut_div_num, hindgut_div_num, gill_div_num, skin_results, midgut_results, hindgut_results, gill_results, file = "alpha_loop_results.RData") 

# load results 
load("alpha_loop_results.RData")

library(ggplot2)

# example visualization for categorical data (analyzed by Kruskal-Wallis)
kw <- ggplot(gill_div_num) + geom_boxplot(aes(y = Shannon, x = method_gear))
kw

# example visualization for continuous data (analyzed by Spearman)
sp <- ggplot(gill_div_num, aes(x = gi_cm, y = Shannon)) + geom_point() + geom_smooth()
sp






