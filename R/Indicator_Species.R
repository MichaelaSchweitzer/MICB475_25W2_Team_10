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

# Hindgut analysis 
mpt_genus_hindgut <- tax_glom(hindgut_final, "Genus", NArm = FALSE)
mpt_genus_RA_hindgut <- transform_sample_counts(mpt_genus_hindgut, fun=function(x) x/sum(x))

isa_hindgut <- multipatt(t(otu_table(mpt_genus_RA_hindgut)), 
                         cluster = sample_data(mpt_genus_RA_hindgut)$method_gear)
summary(isa_hindgut)

taxtable_hindgut <- tax_table(hindgut_final) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ASV")

indic_res_hindgut <- isa_hindgut$sign %>%
  rownames_to_column(var = "ASV") %>%
  left_join(taxtable_hindgut) %>%
  filter(p.value < 0.05) %>%
  pivot_longer(cols = starts_with("s."),
               names_to = "method_gear",
               values_to = "indicator") %>%
  filter(indicator == 1) %>%
  mutate(method_gear = gsub("s\\.", "", method_gear)) 

View(indic_res_hindgut)

# Gill analysis 
mpt_genus_gill <- tax_glom(gill_final, "Genus", NArm = FALSE)
mpt_genus_RA_gill <- transform_sample_counts(mpt_genus_gill, fun=function(x) x/sum(x))

isa_gill <- multipatt(t(otu_table(mpt_genus_RA_gill)), 
                         cluster = sample_data(mpt_genus_RA_gill)$method_gear)
summary(isa_gill)

taxtable_gill <- tax_table(gill_final) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ASV")

indic_res_gill <- isa_gill$sign %>%
  rownames_to_column(var = "ASV") %>%
  left_join(taxtable_gill) %>%
  filter(p.value < 0.05) %>%
  pivot_longer(cols = starts_with("s."),
               names_to = "method_gear",
               values_to = "indicator") %>%
  filter(indicator == 1) %>%
  mutate(method_gear = gsub("s\\.", "", method_gear)) 

View(indic_res_gill)

# Skin analysis 
mpt_genus_skin <- tax_glom(skin_final, "Genus", NArm = FALSE)
mpt_genus_RA_skin <- transform_sample_counts(mpt_genus_skin, fun=function(x) x/sum(x))

isa_skin <- multipatt(t(otu_table(mpt_genus_RA_skin)), 
                      cluster = sample_data(mpt_genus_RA_skin)$method_gear)
summary(isa_skin)

taxtable_skin <- tax_table(skin_final) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ASV")

indic_res_skin <- isa_skin$sign %>%
  rownames_to_column(var = "ASV") %>%
  left_join(taxtable_skin) %>%
  filter(p.value < 0.05) %>%
  pivot_longer(cols = starts_with("s."),
               names_to = "method_gear",
               values_to = "indicator") %>%
  filter(indicator == 1) %>%
  mutate(method_gear = gsub("s\\.", "", method_gear)) 

View(indic_res_skin)

# Data visualization with dot plot 
indic_res_midgut_clean <- indic_res_midgut %>%
  filter(!is.na(Genus),
         Genus != "Incertae_Sedis",
         !grepl("Incertae", Genus)) %>%
  mutate(Genus = gsub("g__", "", Genus))

indic_res_midgut_clean %>%
  ggplot(aes(x = stat, 
             y = reorder(Genus, as.numeric(as.factor(method_gear)) * 100 + stat),
             color = method_gear, )) +
  geom_point() +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(
    axis.text.y = element_text(face = "italic", size = 9),
    axis.text.x = element_text(size = 9),
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(x = "Indicator Value (stat)", 
       y = NULL,
       color = "Gear Type",
       title = "Midgut Indicator Taxa by Gear Type")

ggsave("midgut_ISA_dot.png", 
       width = 10, 
       height = 8)

indic_res_hindgut_clean <- indic_res_hindgut %>%
  filter(!is.na(Genus),
         Genus != "Incertae_Sedis",
         !grepl("Incertae", Genus)) %>%
  mutate(Genus = gsub("g__", "", Genus))

indic_res_hindgut_clean %>%
  ggplot(aes(x = stat, 
             y = reorder(Genus, as.numeric(as.factor(method_gear)) * 100 + stat),
             color = method_gear, )) +
  geom_point() +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(
    axis.text.y = element_text(face = "italic", size = 9),
    axis.text.x = element_text(size = 9),
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(x = "Indicator Value (stat)", 
       y = NULL,
       color = "Gear Type",
       title = "Hindgut Indicator Taxa by Gear Type")

ggsave("hindgut_ISA_dot.png", 
       width = 10, 
       height = 8)

indic_res_gill_clean <- indic_res_gill %>%
  filter(!is.na(Genus),
         Genus != "Incertae_Sedis",
         !grepl("Incertae", Genus)) %>%
  mutate(Genus = gsub("g__", "", Genus))

indic_res_gill_clean %>%
  ggplot(aes(x = stat, 
             y = reorder(Genus, as.numeric(as.factor(method_gear)) * 100 + stat),
             color = method_gear, )) +
  geom_point() +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(
    axis.text.y = element_text(face = "italic", size = 9),
    axis.text.x = element_text(size = 9),
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(x = "Indicator Value (stat)", 
       y = NULL,
       color = "Gear Type",
       title = "Gill Indicator Taxa by Gear Type")

ggsave("gill_ISA_dot.png", 
       width = 10, 
       height = 8)

indic_res_skin_clean <- indic_res_skin %>%
  filter(!is.na(Genus),
         Genus != "Incertae_Sedis",
         !grepl("Incertae", Genus)) %>%
  mutate(Genus = gsub("g__", "", Genus))

indic_res_skin_clean %>%
  ggplot(aes(x = stat, 
             y = reorder(Genus, as.numeric(as.factor(method_gear)) * 100 + stat),
             color = method_gear, )) +
  geom_point() +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(
    axis.text.y = element_text(face = "italic", size = 9),
    axis.text.x = element_text(size = 9),
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(x = "Indicator Value (stat)", 
       y = NULL,
       color = "Gear Type",
       title = "Skin Indicator Taxa by Gear Type")

ggsave("skin_ISA_dot.png", 
       width = 10, 
       height = 8)
