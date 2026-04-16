library(phyloseq)
library(ape)
library(tidyverse)

# Load data
load("all_phyloseq_objects_rarefied.RData")

# Metadata
midgut_meta <- data.frame(sample_data(midgut_rare))
midgut_meta$SampleID <- rownames(midgut_meta)

# Clean method_gear
midgut_meta <- midgut_meta %>%
  filter(!is.na(method_gear),
         method_gear != "not applicable")

# Prune phyloseq
midgut_ps <- prune_samples(midgut_meta$SampleID, midgut_rare)
midgut_ps <- prune_taxa(taxa_sums(midgut_ps) > 0, midgut_ps)

# Reorder metadata
midgut_meta <- midgut_meta[match(sample_names(midgut_ps), midgut_meta$SampleID), ]

# Distance + PCoA
dist_mat <- phyloseq::distance(midgut_ps, method = "wunifrac")
pcoa_res <- ape::pcoa(as.matrix(dist_mat))

# Coordinates
pcoa_df <- data.frame(
  SampleID = rownames(pcoa_res$vectors),
  Axis1 = pcoa_res$vectors[,1],
  Axis2 = pcoa_res$vectors[,2]
) %>%
  left_join(midgut_meta, by = "SampleID")

# Variance explained
var_exp <- round(pcoa_res$values$Relative_eig[1:2] * 100, 1)

# Custom soft palette 
custom_colors <- c(
  "benthic trawl" = "#E69F00",
  "gillnet" = "#56B4E9",
  "midwater trawl" = "#009E73",
  "net" = "#F0E442",
  "rod and reel" = "#0072B2",
  "spear" = "#CC79A7"
)

# Plot
ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = method_gear)) +
  geom_point(size = 3, alpha = 0.85) +
  
  scale_color_manual(values = custom_colors) +
  
  labs(
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    color = "Method gear"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

