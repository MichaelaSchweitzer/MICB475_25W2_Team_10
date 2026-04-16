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
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(ggpubr) # for plotting spearman correlation values 
library(ggsignif)
library(svglite) # for exporting as svg 

#### CREATING THE PHYLOSEQ OBJECT (March 1st and March 2nd, 2026 by Michaela, Updated March 27th, 2026 by Mirren) ####

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

# Removing non-bacterial DNA from all phyloseq objects and rarefying all to a sample size of 1000.
hindgut_filt <- subset_taxa(hindgut_phyloseq_object,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
hindgut_filt_nolow <- filter_taxa(hindgut_filt, function(x) sum(x)>5, prune = TRUE)
hindgut_final <- prune_samples(sample_sums(hindgut_filt_nolow)>100, hindgut_filt_nolow)
hindgut_rare <- rarefy_even_depth(hindgut_final, rngseed = 1, sample.size = 1000)

midgut_filt <- subset_taxa(midgut_phyloseq_object,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
midgut_filt_nolow <- filter_taxa(midgut_filt, function(x) sum(x)>5, prune = TRUE)
midgut_final <- prune_samples(sample_sums(midgut_filt_nolow)>100, midgut_filt_nolow)
midgut_rare <- rarefy_even_depth(midgut_final, rngseed = 1, sample.size = 1000)

skin_filt <- subset_taxa(skin_phyloseq_object,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
skin_filt_nolow <- filter_taxa(skin_filt, function(x) sum(x)>5, prune = TRUE)
skin_final <- prune_samples(sample_sums(skin_filt_nolow)>100, skin_filt_nolow)
skin_rare <- rarefy_even_depth(skin_final, rngseed = 1, sample.size = 1000)

gill_filt <- subset_taxa(gill_phyloseq_object,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
gill_filt_nolow <- filter_taxa(gill_filt, function(x) sum(x)>5, prune = TRUE)
gill_final <- prune_samples(sample_sums(gill_filt_nolow)>100, gill_filt_nolow)
gill_rare <- rarefy_even_depth(gill_final, rngseed = 1, sample.size = 1000)

#### Visualizing components of the non-rarefied phyloseq objects ####
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

#### Visualizing components of the rarefied phyloseq objects ####
# Gill 
otu_table(gill_rare)
sample_data(gill_rare)
tax_table(gill_rare)
phy_tree(gill_rare)

# Hindgut
otu_table(hindgut_rare)
sample_data(hindgut_rare)
tax_table(hindgut_rare)
phy_tree(hindgut_rare)

# Midgut 
otu_table(midgut_rare)
sample_data(midgut_rare)
tax_table(midgut_rare)
phy_tree(midgut_rare)

# Skin 
otu_table(skin_rare)
sample_data(skin_rare)
tax_table(skin_rare)
phy_tree(skin_rare)

# Saving Objects (Can now just load these instead of running all of the above code)
save(skin_phyloseq_object, midgut_phyloseq_object, hindgut_phyloseq_object, gill_phyloseq_object, file = "all_phyloseq_objects_non_rarefied.RData") 
save(skin_final, midgut_final, hindgut_final, gill_final, file = "all_phyloseq_objects_bacteriaonly_pruned.RData")
save(skin_rare, midgut_rare, hindgut_rare, gill_rare, file = "all_phyloseq_objects_rarefied.RData") 

#### Computing alpha diversity metrics ####

# Get alpha diversity metrics.
alpha_div_skin <- estimate_richness(skin_rare, measures = "shannon")
alpha_div_hindgut <- estimate_richness(hindgut_rare, measures = "shannon")
alpha_div_midgut <- estimate_richness(midgut_rare, measures = "shannon")
alpha_div_gill <- estimate_richness(gill_rare, measures = "shannon")

#Create a data frame containing sample data
skin_data <- sample_data(skin_rare)
hindgut_data <- sample_data(hindgut_rare)
midgut_data <- sample_data(midgut_rare)
gill_data <- sample_data(gill_rare)

#Add Shannon scores to sample data 
skin_div <- cbind(alpha_div_skin, skin_data)
hindgut_div <- cbind(alpha_div_hindgut, hindgut_data)
gill_div <- cbind(alpha_div_gill, gill_data)
midgut_div <- cbind(alpha_div_midgut, midgut_data)

### Before we iterate through our phyloseq object, we have to first remove all "not applicable"'s, and convert numbers to numeric values (instead of characters)
# Replace "not applicable"'s to NA (to do this we need to remove date formatting first because R gets confused if not)
skin_div_nodate <- select(skin_div, -"collection_date")
skin_div_nodate[skin_div_nodate == "not applicable"] <- NA

hindgut_div_nodate <- select(hindgut_div, -"collection_date")
hindgut_div_nodate[hindgut_div_nodate == "not applicable"] <- NA

midgut_div_nodate <- select(midgut_div, -"collection_date")
midgut_div_nodate[midgut_div_nodate == "not applicable"] <- NA

gill_div_nodate <- select(gill_div, -"collection_date")
gill_div_nodate[gill_div_nodate == "not applicable"] <- NA

# Make a vector of all the numerical variables
num_vars <- c("dist_to_dorsal_cm", "fl_cm", "gape_cm", "gi_cm", "host_body_mass_index",
             "host_height", "host_height_vs_max_height_tl", "mass_g", "month", "ratio_dorsal_to_tl", "ratio_gape_to_tl", 
             "ratio_gi_to_tl", "tl_cm")

# For loop to make all numerical variables numeric class
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

#### For loop to calculate significance for alpha diversity metrics ####

#Create an empty data frame to hold P values 
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

# Saving the alpha diversity and loop results as an object 
save(skin_div_num, midgut_div_num, hindgut_div_num, gill_div_num, skin_results, midgut_results, hindgut_results, gill_results, file = "alpha_loop_results.RData")

write.csv(skin_results, "skin_alpha_results.csv", row.names = FALSE)
write.csv(hindgut_results, "hindgut_alpha_results.csv", row.names = FALSE)
write.csv(midgut_results, "midgut_alpha_results.csv", row.names = FALSE)
write.csv(gill_results, "gill_alpha_results.csv", row.names = FALSE)

#### P-value table figure ####
# Creating a table to look at p-values for each variable across each organ 
skin_wo_method <- skin_results[,-3] %>%
  setNames(c("variable", "Skin"))
gill_wo_method <- gill_results[,-3] %>%
  setNames(c("variable", "Gill"))
midgut_wo_method <- midgut_results[,-3] %>%
  setNames(c("variable", "Midgut")) 
hindgut_wo_method <- hindgut_results[,-3] %>%
  setNames(c("variable", "Hindgut"))

# Merging the tables for each organ into a single table 
merged_table <- list(skin_wo_method, gill_wo_method, midgut_wo_method, hindgut_wo_method) %>%
  reduce(full_join, by = "variable") %>%
  

# Saving the merged table 
write.csv(merged_table, "alpha_div_pvalues.csv", row.names = FALSE)

#### Within variable alpha diversity plots ####
# Skin, dist_to_dorsal_cm
skin_dist_dorsal <- ggplot(skin_div_num, aes(x = dist_to_dorsal_cm, y = Shannon)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  xlab("Distance to Dorsal (cm)") + 
  ylab("Shannon Index") +
  stat_cor(method = "spearman", label.x = 30) 
skin_dist_dorsal 

# Gill, dist_to_dorsal_cm
gill_dist_dorsal <- ggplot(gill_div_num, aes(x = dist_to_dorsal_cm, y = Shannon)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  xlab("Distance to Dorsal (cm)") + 
  ylab("Shannon Index") +
  stat_cor(method = "spearman", label.x = 30) 
gill_dist_dorsal 

# Gill, fl_cm
gill_fl_cm <- ggplot(gill_div_num, aes(x = fl_cm, y = Shannon)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  xlab("Fork Length (cm)") + 
  ylab("Shannon Index") + 
  stat_cor(method = "spearman", label.x = 75) 
gill_fl_cm 

# Gill, gape_cm
gill_gape_cm <- ggplot(gill_div_num, aes(x = gape_cm, y = Shannon)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  xlab("Gape Length (cm)") + 
  ylab("Shannon Index") + 
  stat_cor(method = "spearman", label.x = 7) 
gill_gape_cm 

# Gill, gi_cm 
gill_gi_cm <- ggplot(gill_div_num, aes(x = gi_cm, y = Shannon)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  xlab("Length of GI Tract (cm)") +
  ylab("Shannon Index") +
  stat_cor(method = "spearman", label.x = 45) 
gill_gi_cm 

# Gill, host_height 
gill_host_height <- ggplot(gill_div_num, aes(x = host_height, y = Shannon)) +
  geom_point() + 
  geom_smooth(method = "lm") + 
  xlab("Host Height (cm)") + 
  ylab("Shannon Index") +
  stat_cor(method = "spearman", label.x = 80) 
gill_host_height

# Gill, mass_g 
gill_mass_g <- ggplot(gill_div_num, aes(x = mass_g, y = Shannon)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  xlab("Mass (g)") +
  ylab("Shannon Index") +
  stat_cor(method = "spearman", label.x = 2500) 
gill_mass_g

# Gill, method_gear 
lm_gill_method_gear <- lm(log(Shannon) ~ `method_gear`, data=gill_div_num)
anova_gill_log <- aov(lm_gill_method_gear)
summary(anova_gill_log)
TukeyHSD(anova_gill_log)

gill_method_gear <- ggplot(gill_div_num, aes(x = method_gear, y = Shannon)) + 
  geom_boxplot() +
  xlab("Method Gear") +
  ylab("Shannon Index") +
  geom_signif(comparisons = list(c("rod and reel","net")),
              y_position = 5,
              annotations = "*") 
gill_method_gear

ggsave(gill_method_gear, file = "Gill_Method_Gear_Alpha.svg")

# Gill, tl_cm 
gill_tail_cm <- ggplot(gill_div_num, aes(x = tl_cm, y = Shannon)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  xlab("Tail Length (cm)") +
  ylab("Shannon Index") +
  stat_cor(method = "spearman", label.x = 80) 
gill_tail_cm

# Midgut, method_gear 
lm_midgut_method_gear <- lm(log(Shannon) ~ `method_gear`, data=midgut_div_num)
anova_midgut_log <- aov(lm_midgut_method_gear)
summary(anova_midgut_log)
TukeyHSD(anova_midgut_log)

midgut_method_gear <- ggplot(midgut_div_num, aes(x = method_gear, y = Shannon)) + 
  geom_boxplot() +
  xlab("Method Gear") +
  ylab("Shannon Index") +
  geom_signif(comparisons = list(c("rod and reel","benthic trawl")),
              y_position = 5,
              annotations = "*") 
midgut_method_gear 

ggsave(midgut_method_gear, file = "Midgut_Method_Gear_Alpha.svg")

# Midgut, month 
midgut_div_num$month <- factor(midgut_div_num$month)
midgut_month <- ggplot(midgut_div_num) + 
  geom_boxplot(aes(y = Shannon, x = month)) +
  xlab("Month") +
  ylab("Shannon Index")
midgut_month

# Midgut, ratio_dorsal_to_tl
midgut_dorsal_tl <- ggplot(midgut_div_num, aes(x = ratio_dorsal_to_tl, y = Shannon)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  xlab("Ratio of dorsal length to tail length") +
  ylab("Shannon Index") +
  stat_cor(method = "spearman", label.x = 0.5) 
midgut_dorsal_tl

# Midgut, gape_to_tl
midgut_gape_tl <- ggplot(midgut_div_num, aes(x = ratio_gape_to_tl, y = Shannon)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  xlab("Ratio of gape length to tail length") +
  ylab("Shannon Index") +
  stat_cor(method = "spearman", label.x = 0.2) 
midgut_gape_tl

# Hindgut, method_gear 
lm_hindgut_method_gear <- lm(log(Shannon) ~ `method_gear`, data=hindgut_div_num)
anova_hindgut_log <- aov(lm_hindgut_method_gear)
summary(anova_hindgut_log)
TukeyHSD(anova_hindgut_log)

hindgut_method_gear <- ggplot(hindgut_div_num, aes(x = method_gear, y = Shannon)) + 
  geom_boxplot() +
  xlab("Method Gear") +
  ylab("Shannon Index") 
hindgut_method_gear 

ggsave(hindgut_method_gear, file = "Hindgut_Method_Gear_Alpha.svg")

# Hindgut, month 
hindgut_div_num$month <- factor(hindgut_div_num$month)
hindgut_month <- ggplot(hindgut_div_num) + 
  geom_boxplot(aes(y = Shannon, x = month)) +
  xlab("Month") +
  ylab("Shannon Index")
hindgut_month

# Hindgut, gape_to_tl
hindgut_gape_tl <- ggplot(hindgut_div_num, aes(x = ratio_gape_to_tl, y = Shannon)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  xlab("Ratio of gape length to tail length") +
  ylab("Shannon Index") +
  stat_cor(method = "spearman", label.x = 0.2) 
hindgut_gape_tl

# Skin, method_gear
lm_skin_method_gear <- lm(log(Shannon) ~ `method_gear`, data=skin_div_num)
anova_skin_log <- aov(lm_skin_method_gear)
summary(anova_skin_log)
TukeyHSD(anova_skin_log)

skin_method_gear <- ggplot(skin_div_num, aes(x = method_gear, y = Shannon)) + 
  geom_boxplot() +
  xlab("Method Gear") +
  ylab("Shannon Index") 
skin_method_gear 

ggsave(skin_method_gear, file = "Skin_Method_Gear_Alpha.svg")


