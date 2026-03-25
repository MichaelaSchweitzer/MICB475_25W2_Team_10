#!/usr/bin/env Rscript

#### Install phyloseq package before beginning, if required. ####
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("phyloseq")

#### Load required libraries ####
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

#### Load filtered phyloseq objects from Aim 1 ####
load("all_phyloseq_objects_rodonly.RData")

#### Create data frames from sample_data ####
skin_div_num <- data.frame(sample_data(skin_rodonly_phyloseq_object))
hindgut_div_num <- data.frame(sample_data(hindgut_rodonly_phyloseq_object))
midgut_div_num <- data.frame(sample_data(midgut_rodonly_phyloseq_object))
gill_div_num <- data.frame(sample_data(gill_rodonly_phyloseq_object))

#### Keep sample IDs as an explicit column ####
skin_div_num$SampleID <- rownames(skin_div_num)
hindgut_div_num$SampleID <- rownames(hindgut_div_num)
midgut_div_num$SampleID <- rownames(midgut_div_num)
gill_div_num$SampleID <- rownames(gill_div_num)

#### Remove collection_date if present ####
if ("collection_date" %in% colnames(skin_div_num)) {
  skin_div_num <- skin_div_num %>% select(-collection_date)
}
if ("collection_date" %in% colnames(hindgut_div_num)) {
  hindgut_div_num <- hindgut_div_num %>% select(-collection_date)
}
if ("collection_date" %in% colnames(midgut_div_num)) {
  midgut_div_num <- midgut_div_num %>% select(-collection_date)
}
if ("collection_date" %in% colnames(gill_div_num)) {
  gill_div_num <- gill_div_num %>% select(-collection_date)
}

#### Replace "not applicable" with NA BEFORE numeric conversion ####
skin_div_num[skin_div_num == "not applicable"] <- NA
hindgut_div_num[hindgut_div_num == "not applicable"] <- NA
midgut_div_num[midgut_div_num == "not applicable"] <- NA
gill_div_num[gill_div_num == "not applicable"] <- NA

#### Make a vector of all numeric variables ####
num_vars <- c(
  "dist_to_dorsal_cm", "fl_cm", "gape_cm", "gi_cm", "host_body_mass_index",
  "host_height", "host_height_vs_max_height_tl", "mass_g", "month",
  "ratio_dorsal_to_tl", "ratio_gape_to_tl", "ratio_gi_to_tl", "tl_cm"
)

#### Create a vector containing the variables of interest for Aim 2 ####
# method_gear removed because all samples are now rod and reel only
aim2_vars <- c(
  "dist_to_dorsal_cm", "fl_cm", "gape_cm", "gi_cm", "host_body_mass_index",
  "host_height", "host_height_vs_max_height_tl", "mass_g",
  "month", "ratio_dorsal_to_tl", "ratio_gape_to_tl", "ratio_gi_to_tl",
  "swim_mode", "swim_performance", "tl_cm"
)

#### Convert numeric variables ####
for (n in num_vars) {
  if (n %in% colnames(skin_div_num)) {
    skin_div_num[[n]] <- as.numeric(as.character(skin_div_num[[n]]))
  }
}
for (n in num_vars) {
  if (n %in% colnames(hindgut_div_num)) {
    hindgut_div_num[[n]] <- as.numeric(as.character(hindgut_div_num[[n]]))
  }
}
for (n in num_vars) {
  if (n %in% colnames(midgut_div_num)) {
    midgut_div_num[[n]] <- as.numeric(as.character(midgut_div_num[[n]]))
  }
}
for (n in num_vars) {
  if (n %in% colnames(gill_div_num)) {
    gill_div_num[[n]] <- as.numeric(as.character(gill_div_num[[n]]))
  }
}

#### Create empty results data frames ####
skin_beta_results <- data.frame(
  variable = character(), p_value = numeric(), r2 = numeric(),
  f_statistic = numeric(), method = character(),
  stringsAsFactors = FALSE
)

hindgut_beta_results <- data.frame(
  variable = character(), p_value = numeric(), r2 = numeric(),
  f_statistic = numeric(), method = character(),
  stringsAsFactors = FALSE
)

midgut_beta_results <- data.frame(
  variable = character(), p_value = numeric(), r2 = numeric(),
  f_statistic = numeric(), method = character(),
  stringsAsFactors = FALSE
)

gill_beta_results <- data.frame(
  variable = character(), p_value = numeric(), r2 = numeric(),
  f_statistic = numeric(), method = character(),
  stringsAsFactors = FALSE
)

#### Helper function for Aim 2 beta diversity testing ####
run_beta_loop <- function(metadata_df, physeq_obj, vars, num_vars) {
  
  results <- data.frame(
    variable = character(), p_value = numeric(),
    r2 = numeric(), f_statistic = numeric(),
    method = character(), stringsAsFactors = FALSE
  )
  
  for (v in vars) {
    
    # Skip if variable is missing from metadata
    if (!(v %in% colnames(metadata_df))) next
    
    # Keep only samples with data for this variable
    meta_sub <- metadata_df[!is.na(metadata_df[[v]]), , drop = FALSE]
    
    # Skip if too few samples
    if (nrow(meta_sub) < 4) next
    
    # Prune phyloseq object to those samples
    ps_sub <- prune_samples(meta_sub$SampleID, physeq_obj)
    ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
    
    # Skip if too few samples or taxa remain
    if (nsamples(ps_sub) < 4 || ntaxa(ps_sub) < 2) next
    
    # Reorder metadata to match phyloseq sample order
    meta_sub <- meta_sub[match(sample_names(ps_sub), meta_sub$SampleID), , drop = FALSE]
    
    # For categorical variables, make sure there are at least 2 groups
    if (!(v %in% num_vars)) {
      meta_sub[[v]] <- factor(meta_sub[[v]])
      meta_sub[[v]] <- droplevels(meta_sub[[v]])
      
      if (nlevels(meta_sub[[v]]) < 2) next
      if (any(table(meta_sub[[v]]) < 2)) next
    }
    
    # For numeric variables, make sure there is variation
    if (v %in% num_vars) {
      if (length(unique(meta_sub[[v]])) < 2) next
    }
    
    # Calculate weighted UniFrac distance
    dist_mat <- phyloseq::distance(ps_sub, method = "wunifrac")
    
    # Run PERMANOVA
    adonis_res <- adonis2(dist_mat ~ meta_sub[[v]], permutations = 999)
    
    # Save results
    results <- rbind(
      results,
      data.frame(
        variable = v,
        p_value = adonis_res$`Pr(>F)`[1],
        r2 = adonis_res$R2[1],
        f_statistic = adonis_res$F[1],
        method = "Weighted UniFrac PERMANOVA",
        stringsAsFactors = FALSE
      )
    )
  }
  
  if (nrow(results) > 0) {
    results$p_adj_bh <- p.adjust(results$p_value, method = "BH")
    results <- results %>% arrange(p_value)
  }
  
  return(results)
}

#### Optional cleaning step for gill metadata ####
# Count number of NAs per sample (only for Aim 2 variables)
gill_div_num$na_count <- rowSums(is.na(gill_div_num[, aim2_vars]))

# Keep samples with fewer than 50% missing values
threshold <- length(aim2_vars) * 0.5
gill_div_num_clean <- gill_div_num[gill_div_num$na_count <= threshold, ]

# Drop helper column
gill_div_num_clean$na_count <- NULL

# Prune phyloseq object to match cleaned metadata
gill_rodonly_phyloseq_object_clean <- prune_samples(
  gill_div_num_clean$SampleID,
  gill_rodonly_phyloseq_object
)

#### Run loop for each body site ####
skin_beta_results <- run_beta_loop(
  skin_div_num, skin_rodonly_phyloseq_object, aim2_vars, num_vars
)

hindgut_beta_results <- run_beta_loop(
  hindgut_div_num, hindgut_rodonly_phyloseq_object, aim2_vars, num_vars
)

midgut_beta_results <- run_beta_loop(
  midgut_div_num, midgut_rodonly_phyloseq_object, aim2_vars, num_vars
)

gill_beta_results <- run_beta_loop(
  gill_div_num_clean, gill_rodonly_phyloseq_object_clean, aim2_vars, num_vars
)

#### View results ####
skin_beta_results
hindgut_beta_results
midgut_beta_results
gill_beta_results

#### Save results ####
save(
  skin_beta_results, hindgut_beta_results, midgut_beta_results, gill_beta_results,
  file = "aim2_weighted_unifrac_results_rodonly.RData"
)

write.csv(skin_beta_results, "skin_weighted_unifrac_results_rodonly.csv", row.names = FALSE)
write.csv(hindgut_beta_results, "hindgut_weighted_unifrac_results_rodonly.csv", row.names = FALSE)
write.csv(midgut_beta_results, "midgut_weighted_unifrac_results_rodonly.csv", row.names = FALSE)
write.csv(gill_beta_results, "gill_weighted_unifrac_results_rodonly.csv", row.names = FALSE)

# Making a merged table with adjusted p-values for beta diversity
skin_padj <- skin_beta_results %>%
  select(variable, p_adj_bh) %>%
  setNames(c("variable", "Skin"))

gill_padj <- gill_beta_results %>%
  select(variable, p_adj_bh) %>%
  setNames(c("variable", "Gill"))

midgut_padj <- midgut_beta_results %>%
  select(variable, p_adj_bh) %>%
  setNames(c("variable", "Midgut"))

hindgut_padj <- hindgut_beta_results %>%
  select(variable, p_adj_bh) %>%
  setNames(c("variable", "Hindgut"))

beta_div_padjvalues <- list(skin_padj, gill_padj, midgut_padj, hindgut_padj) %>%
  reduce(full_join, by = "variable")

beta_div_padjvalues
write.csv(beta_div_padjvalues, "beta_div_padjvalues.csv", row.names = FALSE)

####THIS SECTION ONWARDS IS NOT UPDATED JUST FYI #######

# Plotting Data (March 22nd, Michaela)

# Making a merged table with just Adonis R^2 values 
skin_R2 <- skin_beta_results%>%
           select(-c(2, 4, 5, 6)) %>%
           setNames(c("variable", "Skin"))
gill_R2 <- gill_beta_results%>%
           select(-c(2, 4, 5, 6)) %>%
           setNames(c("variable", "Gill"))
midgut_R2 <- midgut_beta_results%>%
             select(-c(2, 4, 5, 6)) %>%
             setNames(c("variable", "Midgut"))
hindgut_R2 <- hindgut_beta_results%>%
              select(-c(2, 4, 5, 6)) %>%
              setNames(c("variable", "Hindgut"))

merged_table_beta <- list(skin_R2, gill_R2, midgut_R2, hindgut_R2) %>%
  reduce(full_join, by = "variable")

merged_table_beta_long <- merged_table_beta %>% 
                          pivot_longer(
                          cols = c(`Skin`, `Gill`, `Midgut`, `Hindgut`),
                          names_to = "Tissue",
                          values_to = "R2_Value") 
merged_table_beta_long

# Plotting the graph 
R2_values <- ggplot(merged_table_beta_long, aes(x = variable, y = R2_Value, fill = Tissue)) + 
             geom_bar(stat = "identity", position = "dodge") +
             theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
             ylab("Adonis R^2") +
             xlab("Variable")
R2_values

# Saving the graph 
ggsave(file = "betadiv.png", 
      , plot = R2_values, 
      , height=8, width=20)
