#!/usr/bin/env Rscript

# Install phyloseq package before beginning, if required.
# if (!require("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")
# BiocManager::install("phyloseq")

# Install ggpubr
install.packages("ggpubr")

#### Load required libraries. Run each time a new session is opened.####
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)
library(ggplot2)
library(ggpubr)

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

# Saving formatted metadata
write.csv(gill_metadata_phylo, "gillmetadata.csv", row.names = FALSE)
write.csv(hindgut_metadata_phylo, "hindgutmetadata.csv", row.names = FALSE)
write.csv(midgut_metadata_phylo, "midgutmetadata.csv", row.names = FALSE)
write.csv(skin_metadata_phylo, "skinmetadata.csv", row.names = FALSE)

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

write.csv(skin_results, "skin_alpha_results.csv", row.names = FALSE)
write.csv(hindgut_results, "hindgut_alpha_results.csv", row.names = FALSE)
write.csv(midgut_results, "midgut_alpha_results.csv", row.names = FALSE)
write.csv(gill_results, "gill_alpha_results.csv", row.names = FALSE)

# load results 
load("alpha_loop_results.RData")

# example visualization for categorical data (analyzed by Kruskal-Wallis)
kw <- ggplot(gill_div_num) + geom_boxplot(aes(y = Shannon, x = method_gear))

# example visualization for continuous data (analyzed by Spearman)
sp <- ggplot(gill_div_num, aes(x = gi_cm, y = Shannon)) + geom_point() + geom_smooth(method = lm)

# Making plots for alpha diversity (March 22nd, Michaela)

# Plot to count significant variables
skin_wo_method <- skin_results[,-3] %>%
                  setNames(c("variable", "Skin"))
gill_wo_method <- gill_results[,-3] %>%
                  setNames(c("variable", "Gill"))
midgut_wo_method <- midgut_results[,-3] %>%
                    setNames(c("variable", "Midgut")) 
hindgut_wo_method <- hindgut_results[,-3] %>%
                     setNames(c("variable", "Hindgut"))

merged_table <- list(skin_wo_method, gill_wo_method, midgut_wo_method, hindgut_wo_method) %>%
                reduce(full_join, by = "variable")

write.csv(merged_table, "alpha_div_pvalues.csv", row.names = FALSE)

# Bubble plot to show relative p-values 
bubble_plot_data <- merged_table %>%
                      pivot_longer(cols = c("Skin", "Gill", "Midgut", "Hindgut"), names_to = "dataset", values_to = "p_value") %>%
                      mutate(significant = p_value < 0.05, minusLogP = -log10(p_value))

bubble_plot <- ggplot(bubble_plot_data, aes(x = dataset, y = variable, size = minusLogP, color = minusLogP)) +
               geom_point(alpha = 0.7) +
               scale_size(range = c(3, 12)) +                # adjust bubble sizes
               scale_color_gradient(low = "lightblue", high = "darkblue") +  # gradient
               theme_minimal() +
               labs(
               x = "Tissue",
               y = "Variable",
               size = "-log10(p-value)",
               color = "-log10(p-value)"
               )
bubble_plot

ggsave(file = "bubble_plot.png",
              , plot = bubble_plot
              , height=10, width=6)

# Within variable alpha diversity plots
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
              stat_cor(method = "spearman", label.x = 80) 
gill_fl_cm 

# Gill, gape_cm
gill_gape_cm <- ggplot(gill_div_num, aes(x = gape_cm, y = Shannon)) + 
                geom_point() + 
                geom_smooth(method = "lm") + 
                xlab("Gape Length (cm)") + 
                ylab("Shannon Index") + 
                stat_cor(method = "spearman", label.x = 10) 
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
               stat_cor(method = "spearman", label.x = 20000) 
gill_mass_g

# Gill, method_gear 
gill_method_gear <- ggplot(gill_div_num) + 
                    geom_boxplot(aes(y = Shannon, x = method_gear)) +
                    xlab("Method Gear") +
                    ylab("Shannon Index")
gill_method_gear

# Gill, tl_cm 
gill_tail_cm <- ggplot(gill_div_num, aes(x = tl_cm, y = Shannon)) + 
                geom_point() + 
                geom_smooth(method = "lm") + 
                xlab("Tail Length (cm)") +
                ylab("Shannon Index") +
                stat_cor(method = "spearman", label.x = 80) 
gill_tail_cm

# Midgut, method_gear 
midgut_method_gear <- ggplot(midgut_div_num) + 
                      geom_boxplot(aes(y = Shannon, x = method_gear)) +
                      xlab("Method Gear") +
                      ylab("Shannon Index")
midgut_method_gear 

# Midgut, month 
midgut_div_num$month <- factor(midgut_div_num$month)
midgut_month <- ggplot(midgut_div_num) + 
                geom_boxplot(aes(y = Shannon, x = month)) +
                xlab("Month") +
                ylab("Shannon Index")
midgut_month

# Hindgut, host BMI 
hindgut_host_bmi <- ggplot(hindgut_div_num, aes(x = host_body_mass_index, y = Shannon)) + 
                    geom_point() + 
                    geom_smooth(method = "lm") + 
                    xlab("Host Body Mass Index") +
                    ylab("Shannon Index") + 
                    stat_cor(method = "spearman", label.x = 1.7) 
hindgut_host_bmi

# Hindgut, method_gear 
hindgut_method_gear <- ggplot(hindgut_div_num) + 
                       geom_boxplot(aes(y = Shannon, x = method_gear)) +
                       xlab("Method Gear") +
                       ylab("Shannon Index")
hindgut_method_gear 

# Hindgut, month 
hindgut_div_num$month <- factor(hindgut_div_num$month)
hindgut_month <- ggplot(hindgut_div_num) + 
                geom_boxplot(aes(y = Shannon, x = month)) +
                xlab("Month") +
                ylab("Shannon Index")
hindgut_month

# Hindgut, swim_mode
hindgut_swim_mode <- ggplot(hindgut_div_num) + 
                     geom_boxplot(aes(y = Shannon, x = swim_mode)) +
                     xlab("Swim Mode") +
                     ylab("Shannon Index")
hindgut_swim_mode

# Saving plots
ggsave(file = "gilldisttodorsal.png",
       , plot = gill_dist_dorsal
       , height=6, width=10)

ggsave(file = "gillflcm.png",
       , plot = gill_fl_cm
       , height=6, width=10)

ggsave(file = "gillgapecm.png",
       , plot = gill_gape_cm
       , height=6, width=10)

ggsave(file = "gillgicm.png",
       , plot = gill_gi_cm
       , height=6, width=10)

ggsave(file = "gillhostheight.png",
       , plot = gill_host_height
       , height=6, width=10)

ggsave(file = "gillmassg.png",
       , plot = gill_mass_g
       , height=6, width=10)

ggsave(file = "gillmethodgear.png",
       , plot = gill_method_gear
       , height=6, width=10)

ggsave(file = "gilltlcm.png",
       , plot = gill_tail_cm
       , height=6, width=10)

ggsave(file = "midgutmethodgear.png",
       , plot = midgut_method_gear
       , height=6, width=10)

ggsave(file = "midgutmonth.png",
       , plot = midgut_month
       , height=6, width=10)

ggsave(file = "hindgutbmi.png",
       , plot = hindgut_host_bmi
       , height=6, width=10)

ggsave(file = "hindgutmethodgear.png",
       , plot = hindgut_method_gear
       , height=6, width=10)

ggsave(file = "hindgutmonth.png",
       , plot = hindgut_month
       , height=6, width=10)

ggsave(file = "hindgutswimmode.png",
       , plot = hindgut_swim_mode
       , height=6, width=10)

