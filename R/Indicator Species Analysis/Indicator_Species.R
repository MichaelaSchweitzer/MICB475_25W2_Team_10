#### Indicator Species Analysis ####

# Load required libraries
library(phyloseq)
library(tidyverse)
library(indicspecies)

# Load necessary objects
load("all_phyloseq_objects_bacteriaonly_pruned.RData")

# Midgut analysis 
mpt_genus_midgut <- tax_glom(midgut_final, "Genus", NArm = FALSE)
mpt_genus_RA_midgut <- transform_sample_counts(mpt_genus_midgut, fun=function(x) x/sum(x))

isa_midgut <- multipatt(t(otu_table(mpt_genus_RA_midgut)), 
                        cluster = sample_data(mpt_genus_RA_midgut)$method_gear)
summary(isa_midgut)

taxtable_midgut <- tax_table(midgut_final) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ASV")

indic_res_midgut <- isa_midgut$sign %>%
  rownames_to_column(var = "ASV") %>%
  left_join(taxtable_midgut) %>%
  filter(p.value < 0.05) %>%
  pivot_longer(cols = starts_with("s."),
               names_to = "method_gear",
               values_to = "indicator") %>%
  filter(indicator == 1) %>%
  mutate(method_gear = gsub("s\\.", "", method_gear)) 

View(indic_res_midgut)


# Data visualization with dot plot 
indic_res_midgut_clean <- indic_res_midgut %>%
  filter(!is.na(Genus),
         Genus != "Incertae_Sedis",
         !grepl("Incertae", Genus)) %>%
  mutate(Genus = gsub("g__", "", Genus))

indic_res_midgut_clean %>%
  ggplot(aes(x = stat, 
             y = reorder(Genus, as.numeric(as.factor(method_gear)) * 100 + stat),
             color = method_gear)) +
  geom_point(size = 3) +                        # increase point size (default is 1.5)
  scale_color_manual(values = c(
    "#E41A1C",  # red
    "#377EB8",  # blue
    "#4DAF4A",  # green
    "#984EA3",  # purple
    "#A65628",  # brown
    "#F781BF"   # pink
  )) +
  theme_bw() +
  theme(
    axis.text.y = element_text(face = "italic", size = 9),
    axis.text.x = element_text(size = 9),
    legend.position = "right",
    legend.key.size = unit(1.5, "lines"),        # increase legend key size (default is 1)
    legend.text = element_text(size = 12),       # optional: larger legend text
    legend.title = element_text(size = 14),      # optional: larger legend title
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(x = "Indicator Value (stat)", 
       y = NULL,
       color = "Gear Type",
       title = "Midgut Indicator Taxa by Gear Type")

ggsave("midgut_ISA_dot.png", 
       width = 8, 
       height = 12)
