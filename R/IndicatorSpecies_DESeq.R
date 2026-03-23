# Aim 3 - Indicator species and DESeq

#### Load required libraries ####
library(phyloseq)
library(tidyverse)
library(DESeq2)
library(indicspecies)
library(ggplot2)

#### Load objects from Aim 1 ####
load("all_phyloseq_objects.RData")
load("alpha_loop_results.RData")


# Check levels of your significant variables
levels(factor(sample_data(midgut_phyloseq_object)$method_gear))
levels(factor(sample_data(gill_phyloseq_object)$method_gear))

# Remove "not applicable" samples from ALL variables in gill
gill_phyloseq_object_clean <- subset_samples(
  gill_phyloseq_object,
  method_gear != "not applicable" & dist_to_dorsal_cm != "not applicable" & fl_cm != "not applicable" & gape_cm != "not applicable" & gi_cm != "not applicable" & host_height != "not applicable" & mass_g != "not applicable" & tl_cm != "not applicable")


# Remove any taxa that are now all zeros after sample removal
gill_phyloseq_object_clean <- prune_taxa(
  taxa_sums(gill_phyloseq_object_clean) > 0, 
  gill_phyloseq_object_clean
)

# Check how many samples kept vs removed
nsamples(gill_phyloseq_object)       # original
nsamples(gill_phyloseq_object_clean) # after cleaning

# Verify levels in each body site 
levels(factor(sample_data(midgut_phyloseq_object)$method_gear))
levels(factor(sample_data(hindgut_phyloseq_object)$method_gear))
levels(factor(sample_data(gill_phyloseq_object_clean)$method_gear))

#### INDICATOR SPECIES ANALYSES ####

# Midgut 
isa_midgut <- multipatt(t(otu_table(midgut_phyloseq_object)), 
                        cluster = sample_data(midgut_phyloseq_object)$method_gear)
summary(isa_midgut)

taxtable_midgut <- tax_table(midgut_phyloseq_object) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ASV")

indic_res_midgut <- isa_midgut$sign %>%
  rownames_to_column(var = "ASV") %>%
  left_join(taxtable_midgut) %>%
  filter(p.value < 0.05)

View(indic_res_midgut)

# Hindgut
isa_hindgut <- multipatt(t(otu_table(hindgut_phyloseq_object)), 
                         cluster = sample_data(hindgut_phyloseq_object)$method_gear)
summary(isa_hindgut)

taxtable_hindgut <- tax_table(hindgut_phyloseq_object) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ASV")

indic_res_hindgut <- isa_hindgut$sign %>%
  rownames_to_column(var = "ASV") %>%
  left_join(taxtable_hindgut) %>%
  filter(p.value < 0.05)

View(indic_res_hindgut)

# Gill
isa_gill <- multipatt(t(otu_table(gill_phyloseq_object_clean)), 
                      cluster = sample_data(gill_phyloseq_object_clean)$method_gear)
summary(isa_gill)

taxtable_gill <- tax_table(gill_phyloseq_object_clean) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ASV")

indic_res_gill <- isa_gill$sign %>%
  rownames_to_column(var = "ASV") %>%
  left_join(taxtable_gill) %>%
  filter(p.value < 0.05)

View(indic_res_gill)

#### DESeq ANALYSES ####

# Midgut
midgut_plus1 <- transform_sample_counts(midgut_phyloseq_object, function(x) x + 1)
midgut_deseq2 <- phyloseq_to_deseq2(midgut_plus1, ~ method_gear)
midgut_deseq <- DESeq(midgut_deseq2)

deseq_res_midgut <- results(midgut_deseq, tidy = TRUE)
View(deseq_res_midgut)

sigASVs_midgut <- deseq_res_midgut %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row)

sigASVs_vec_midgut <- sigASVs_midgut %>%
  pull(ASV)

midgut_filt <- prune_taxa(sigASVs_vec_midgut, midgut_phyloseq_object)

merged_results_midgut <- tax_table(midgut_phyloseq_object) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_midgut) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

method_gear_midgut <- ggplot(merged_results_midgut) +
  geom_bar(aes(x = Genus, y = log2FoldChange), stat = "identity") +
  geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, 
                    ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Midgut - Differentially Abundant ASVs by Gear Type",
       x = "Genus", y = "Log2 Fold Change")
method_gear_midgut

ggsave(file = "midgutdeseqmethodgear.png", 
       plot = method_gear_midgut,
       height = 8, width = 20)

midgut_deseq2_month <- phyloseq_to_deseq2(midgut_plus1, ~ month)
midgut_deseq_month <- DESeq(midgut_deseq2_month)

deseq_res_midgut_month <- results(midgut_deseq_month, tidy = TRUE)
View(deseq_res_midgut_month)

sigASVs_midgut_month <- deseq_res_midgut_month %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row)

sigASVs_vec_midgut_month <- sigASVs_midgut_month %>%
  pull(ASV)

midgut_filt_month <- prune_taxa(sigASVs_vec_midgut_month, midgut_phyloseq_object)

merged_results_midgut_month <- tax_table(midgut_phyloseq_object) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_midgut_month) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

month_midgut <- ggplot(merged_results_midgut_month) +
  geom_bar(aes(x = Genus, y = log2FoldChange), stat = "identity") +
  geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, 
                    ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Midgut - Differentially Abundant ASVs by Month",
       x = "Genus", y = "Log2 Fold Change")
month_midgut

ggsave(file = "midgutdeseqmonth.png", 
       plot = month_midgut,
       height = 8, width = 20)


# Hindgut
hindgut_plus1 <- transform_sample_counts(hindgut_phyloseq_object, function(x) x + 1)
hindgut_deseq2 <- phyloseq_to_deseq2(hindgut_plus1, ~ method_gear)
hindgut_deseq <- DESeq(hindgut_deseq2)

deseq_res_hindgut <- results(hindgut_deseq, tidy = TRUE)
View(deseq_res_hindgut)

sigASVs_hindgut <- deseq_res_hindgut %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row)

sigASVs_vec_hindgut <- sigASVs_hindgut %>%
  pull(ASV)

hindgut_filt <- prune_taxa(sigASVs_vec_hindgut, hindgut_phyloseq_object)

merged_results_hindgut <- tax_table(hindgut_phyloseq_object) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_hindgut) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

method_gear_hindgut_plot <- ggplot(merged_results_hindgut) +
  geom_bar(aes(x = Genus, y = log2FoldChange), stat = "identity") +
  geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, 
                    ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Hindgut - Differentially Abundant ASVs by Gear Type",
       x = "Genus", y = "Log2 Fold Change")
method_gear_hindgut_plot

ggsave(file = "hindgutdeseqmethodgear.png", 
       plot = method_gear_hindgut_plot,
       height = 8, width = 20)

hindgut_deseq2_bmi <- phyloseq_to_deseq2(hindgut_plus1, ~ host_body_mass_index)
hindgut_deseq_bmi <- DESeq(hindgut_deseq2_bmi)

deseq_res_hindgut_bmi <- results(hindgut_deseq_bmi, tidy = TRUE)
View(deseq_res_hindgut_bmi)

sigASVs_hindgut_bmi <- deseq_res_hindgut_bmi %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row)

sigASVs_vec_hindgut_bmi <- sigASVs_hindgut_bmi %>%
  pull(ASV)

hindgut_filt_bmi <- prune_taxa(sigASVs_vec_hindgut_bmi, hindgut_phyloseq_object)

merged_results_hindgut_bmi <- tax_table(hindgut_phyloseq_object) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_hindgut_bmi) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

bmi_hindgut_plot <- ggplot(merged_results_hindgut_bmi) +
  geom_bar(aes(x = Genus, y = log2FoldChange), stat = "identity") +
  geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, 
                    ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Hindgut - Differentially Abundant ASVs by Host Body Mass Index",
       x = "Genus", y = "Log2 Fold Change")
bmi_hindgut_plot

ggsave(file = "hindgutdeseqbmi.png", 
       plot = bmi_hindgut_plot,
       height = 8, width = 20)

hindgut_deseq2_month <- phyloseq_to_deseq2(hindgut_plus1, ~ month)
hindgut_deseq_month <- DESeq(hindgut_deseq2_month)

deseq_res_hindgut_month <- results(hindgut_deseq_month, tidy = TRUE)
View(deseq_res_hindgut_month)

sigASVs_hindgut_month <- deseq_res_hindgut_month %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row)

sigASVs_vec_hindgut_month <- sigASVs_hindgut_month %>%
  pull(ASV)

hindgut_filt_month <- prune_taxa(sigASVs_vec_hindgut_month, hindgut_phyloseq_object)

merged_results_hindgut_month <- tax_table(hindgut_phyloseq_object) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_hindgut_month) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

month_hindgut_plot <- ggplot(merged_results_hindgut_month) +
  geom_bar(aes(x = Genus, y = log2FoldChange), stat = "identity") +
  geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, 
                    ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Hindgut - Differentially Abundant ASVs by Month",
       x = "Genus", y = "Log2 Fold Change")
month_hindgut_plot

ggsave(file = "hindgutdeseqmonth.png", 
       plot = month_hindgut_plot,
       height = 8, width = 20)

hindgut_deseq2_swim <- phyloseq_to_deseq2(hindgut_plus1, ~ swim_mode)
hindgut_deseq_swim <- DESeq(hindgut_deseq2_swim)

deseq_res_hindgut_swim <- results(hindgut_deseq_swim, tidy = TRUE)
View(deseq_res_hindgut_swim)

sigASVs_hindgut_swim <- deseq_res_hindgut_swim %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row)

sigASVs_vec_hindgut_swim <- sigASVs_hindgut_swim %>%
  pull(ASV)

hindgut_filt_swim <- prune_taxa(sigASVs_vec_hindgut_swim, hindgut_phyloseq_object)

merged_results_hindgut_swim <- tax_table(hindgut_phyloseq_object) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_hindgut_swim) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

swim_hindgut_plot <- ggplot(merged_results_hindgut_swim) +
  geom_bar(aes(x = Genus, y = log2FoldChange), stat = "identity") +
  geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, 
                    ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Hindgut - Differentially Abundant ASVs by Swim Mode",
       x = "Genus", y = "Log2 Fold Change")
swim_hindgut_plot

ggsave(file = "hindgutdeseqswim.png", 
       plot = swim_hindgut_plot,
       height = 8, width = 20)

# Gill
gill_plus1 <- transform_sample_counts(gill_phyloseq_object_clean, function(x) x + 1)
gill_deseq2 <- phyloseq_to_deseq2(gill_plus1, ~ method_gear)
gill_deseq <- DESeq(gill_deseq2)

deseq_res_gill <- results(gill_deseq, tidy = TRUE)
View(deseq_res_gill)

sigASVs_gill <- deseq_res_gill %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row)

sigASVs_vec_gill <- sigASVs_gill %>%
  pull(ASV)

gill_filt <- prune_taxa(sigASVs_vec_gill, gill_phyloseq_object_clean)

merged_results_gill <- tax_table(gill_phyloseq_object_clean) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_gill) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

gill_method_gear <- ggplot(merged_results_gill) +
  geom_bar(aes(x = Genus, y = log2FoldChange), stat = "identity") +
  geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, 
                    ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Gill - Differentially Abundant ASVs by Gear Type",
       x = "Genus", y = "Log2 Fold Change")
gill_method_gear

ggsave(file = "gilldeseqmethodgear.png", 
       plot = gill_method_gear,
       height = 8, width = 20)

gill_deseq2_distdorsal <- phyloseq_to_deseq2(gill_plus1, ~ dist_to_dorsal_cm)
gill_deseq_distdorsal <- DESeq(gill_deseq2_distdorsal)

deseq_res_gill_distdorsal <- results(gill_deseq_distdorsal, tidy = TRUE)
View(deseq_res_gill_distdorsal)

sigASVs_gill_distdorsal <- deseq_res_gill_distdorsal %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row)

sigASVs_vec_gill_distdorsal <- sigASVs_gill_distdorsal %>%
  pull(ASV)

gill_filt_distdorsal <- prune_taxa(sigASVs_vec_gill_distdorsal, gill_phyloseq_object_clean)

merged_results_gill_distdorsal <- tax_table(gill_phyloseq_object_clean) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_gill_distdorsal) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

gill_distdorsal <- ggplot(merged_results_gill_distdorsal) +
  geom_bar(aes(x = Genus, y = log2FoldChange), stat = "identity") +
  geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, 
                    ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Gill - Differentially Abundant ASVs by Distance to Dorsal",
       x = "Genus", y = "Log2 Fold Change")
gill_distdorsal

ggsave(file = "gilldeseqmethoddistdorsal.png", 
       plot = gill_distdorsal,
       height = 8, width = 20)

gill_deseq2_fl <- phyloseq_to_deseq2(gill_plus1, ~ fl_cm)
gill_deseq_fl <- DESeq(gill_deseq2_fl)

deseq_res_gill_fl <- results(gill_deseq_fl, tidy = TRUE)
View(deseq_res_gill_fl)

sigASVs_gill_fl <- deseq_res_gill_fl %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row)

sigASVs_vec_gill_fl <- sigASVs_gill_fl %>%
  pull(ASV)

gill_filt_fl <- prune_taxa(sigASVs_vec_gill_fl, gill_phyloseq_object_clean)

merged_results_gill_fl <- tax_table(gill_phyloseq_object_clean) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_gill_fl) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

gill_fl <- ggplot(merged_results_gill_fl) +
  geom_bar(aes(x = Genus, y = log2FoldChange), stat = "identity") +
  geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, 
                    ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Gill - Differentially Abundant ASVs by Fork Length",
       x = "Genus", y = "Log2 Fold Change")
gill_fl

ggsave(file = "gilldeseqfl.png", 
       plot = gill_fl,
       height = 8, width = 20)

gill_deseq2_gape <- phyloseq_to_deseq2(gill_plus1, ~ gape_cm)
gill_deseq_gape <- DESeq(gill_deseq2_gape)

deseq_res_gill_gape <- results(gill_deseq_gape, tidy = TRUE)
View(deseq_res_gill_gape)

sigASVs_gill_gape <- deseq_res_gill_gape %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row)

sigASVs_vec_gill_gape <- sigASVs_gill_gape %>%
  pull(ASV)

gill_filt_gape <- prune_taxa(sigASVs_vec_gill_gape, gill_phyloseq_object_clean)

merged_results_gill_gape <- tax_table(gill_phyloseq_object_clean) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_gill_gape) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

gill_gape <- ggplot(merged_results_gill_gape) +
  geom_bar(aes(x = Genus, y = log2FoldChange), stat = "identity") +
  geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, 
                    ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Gill - Differentially Abundant ASVs by Gape Length",
       x = "Genus", y = "Log2 Fold Change")
gill_gape

ggsave(file = "gilldeseqgape.png", 
       plot = gill_gape,
       height = 8, width = 20)

gill_deseq2_gi <- phyloseq_to_deseq2(gill_plus1, ~ gi_cm)
gill_deseq_gi <- DESeq(gill_deseq2_gi)

deseq_res_gill_gi <- results(gill_deseq_gi, tidy = TRUE)
View(deseq_res_gill_gi)

sigASVs_gill_gi <- deseq_res_gill_gi %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row)

sigASVs_vec_gill_gi <- sigASVs_gill_gi %>%
  pull(ASV)

gill_filt_gi <- prune_taxa(sigASVs_vec_gill_gi, gill_phyloseq_object_clean)

merged_results_gill_gi <- tax_table(gill_phyloseq_object_clean) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_gill_gi) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

gill_gi <- ggplot(merged_results_gill_gi) +
  geom_bar(aes(x = Genus, y = log2FoldChange), stat = "identity") +
  geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, 
                    ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Gill - Differentially Abundant ASVs by GI Tract Length",
       x = "Genus", y = "Log2 Fold Change")
gill_gi

ggsave(file = "gilldeseqgi.png", 
       plot = gill_gi,
       height = 8, width = 20)

gill_deseq2_height <- phyloseq_to_deseq2(gill_plus1, ~ host_height)
gill_deseq_height <- DESeq(gill_deseq2_height)

deseq_res_gill_height <- results(gill_deseq_height, tidy = TRUE)
View(deseq_res_gill_height)

sigASVs_gill_height <- deseq_res_gill_height %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row)

sigASVs_vec_gill_height <- sigASVs_gill_height %>%
  pull(ASV)

gill_filt_height <- prune_taxa(sigASVs_vec_gill_height, gill_phyloseq_object_clean)

merged_results_gill_height <- tax_table(gill_phyloseq_object_clean) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_gill_height) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

gill_height <- ggplot(merged_results_gill_height) +
  geom_bar(aes(x = Genus, y = log2FoldChange), stat = "identity") +
  geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, 
                    ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Gill - Differentially Abundant ASVs by GI Tract Length",
       x = "Genus", y = "Log2 Fold Change")
gill_height

ggsave(file = "gilldeseqheight.png", 
       plot = gill_height,
       height = 8, width = 20)

gill_deseq2_mass <- phyloseq_to_deseq2(gill_plus1, ~ mass_g)
gill_deseq_mass <- DESeq(gill_deseq2_mass)

deseq_res_gill_mass <- results(gill_deseq_mass, tidy = TRUE)
View(deseq_res_gill_mass)

sigASVs_gill_mass <- deseq_res_gill_mass %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row)

sigASVs_vec_gill_mass <- sigASVs_gill_mass %>%
  pull(ASV)

gill_filt_mass <- prune_taxa(sigASVs_vec_gill_mass, gill_phyloseq_object_clean)

merged_results_gill_mass <- tax_table(gill_phyloseq_object_clean) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_gill_mass) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

gill_mass <- ggplot(merged_results_gill_mass) +
  geom_bar(aes(x = Genus, y = log2FoldChange), stat = "identity") +
  geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, 
                    ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Gill - Differentially Abundant ASVs by Mass",
       x = "Genus", y = "Log2 Fold Change")
gill_mass

ggsave(file = "gilldeseqmass.png", 
       plot = gill_mass,
       height = 8, width = 20)

gill_deseq2_tail <- phyloseq_to_deseq2(gill_plus1, ~ tl_cm)
gill_deseq_tail <- DESeq(gill_deseq2_tail)

deseq_res_gill_tail <- results(gill_deseq_tail, tidy = TRUE)
View(deseq_res_gill_tail)

sigASVs_gill_tail <- deseq_res_gill_tail %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row)

sigASVs_vec_gill_tail <- sigASVs_gill_tail %>%
  pull(ASV)

gill_filt_tail <- prune_taxa(sigASVs_vec_gill_tail, gill_phyloseq_object_clean)

merged_results_gill_tail <- tax_table(gill_phyloseq_object_clean) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_gill_tail) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

gill_tail <- ggplot(merged_results_gill_tail) +
  geom_bar(aes(x = Genus, y = log2FoldChange), stat = "identity") +
  geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, 
                    ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Gill - Differentially Abundant ASVs by Tail Length",
       x = "Genus", y = "Log2 Fold Change")
gill_tail

ggsave(file = "gilldeseqtail.png", 
       plot = gill_tail,
       height = 8, width = 20)

resultsNames(midgut_deseq)
resultsNames(hindgut_deseq)
resultsNames(gill_deseq)
