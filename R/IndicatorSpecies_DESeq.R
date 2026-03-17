# Aim 3 - Indicator species and DESeq

#### Load required libraries ####
library(phyloseq)
library(tidyverse)
library(DESeq2)
library(indicspecies)

#### Load objects from Aim 1 ####
load("all_phyloseq_objects.RData")
load("alpha_loop_results.RData")


# Check levels of your significant variables
levels(factor(sample_data(midgut_phyloseq_object)$method_gear))
levels(factor(sample_data(gill_phyloseq_object)$method_gear))

# Remove "not applicable" samples from BOTH variables in gill
gill_phyloseq_object_clean <- subset_samples(
  gill_phyloseq_object,
  method_gear != "not applicable" & swim_performance != "not applicable"
)

# Remove any taxa that are now all zeros after sample removal
gill_phyloseq_object_clean <- prune_taxa(
  taxa_sums(gill_phyloseq_object_clean) > 0, 
  gill_phyloseq_object_clean
)

# Verify "not applicable" is gone
levels(factor(sample_data(gill_phyloseq_object_clean)$method_gear))

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

ggplot(merged_results_midgut) +
  geom_bar(aes(x = Genus, y = log2FoldChange), stat = "identity") +
  geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, 
                    ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Midgut - Differentially Abundant ASVs by Gear Type",
       x = "Genus", y = "Log2 Fold Change")

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

ggplot(merged_results_hindgut) +
  geom_bar(aes(x = Genus, y = log2FoldChange), stat = "identity") +
  geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, 
                    ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Hindgut - Differentially Abundant ASVs by Gear Type",
       x = "Genus", y = "Log2 Fold Change")

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

ggplot(merged_results_gill) +
  geom_bar(aes(x = Genus, y = log2FoldChange), stat = "identity") +
  geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, 
                    ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Gill - Differentially Abundant ASVs by Gear Type",
       x = "Genus", y = "Log2 Fold Change")

resultsNames(midgut_deseq)
resultsNames(hindgut_deseq)
resultsNames(gill_deseq)
