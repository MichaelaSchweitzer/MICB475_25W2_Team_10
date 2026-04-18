# Aim 3 - Indicator species and DESeq

#### Load required libraries ####
library(phyloseq)
library(tidyverse)
library(DESeq2)
library(indicspecies)
library(ggplot2)

#### Load objects from Aim 1 ####
load("all_phyloseq_objects_bacteriaonly_pruned.RData")
load("alpha_loop_results.RData")

# Check levels of the method_gear variable across the three tissues where there is significance
levels(factor(sample_data(hindgut_final)$method_gear))
levels(factor(sample_data(midgut_final)$method_gear))
levels(factor(sample_data(gill_final)$method_gear))

# Remove "not applicable" samples from ALL variables in gill
gill_phyloseq_object_clean <- subset_samples(
  gill_final,
  method_gear != "not applicable")

# Remove any taxa that are now all zeros after sample removal
gill_phyloseq_object_clean <- prune_taxa(
  taxa_sums(gill_phyloseq_object_clean) > 0, 
  gill_phyloseq_object_clean
)

# Check how many samples kept vs removed
nsamples(gill_final)       # original
nsamples(gill_phyloseq_object_clean) # after cleaning

# Verify levels in each body site 
levels(factor(sample_data(midgut_final)$method_gear))
levels(factor(sample_data(hindgut_final)$method_gear))
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
midgut_plus1 <- transform_sample_counts(midgut_final, function(x) x + 1)
midgut_deseq2 <- phyloseq_to_deseq2(midgut_plus1, ~ method_gear)
midgut_deseq <- DESeq(midgut_deseq2)

deseq_res_midgut <- results(midgut_deseq, tidy = TRUE)
View(deseq_res_midgut)

sigASVs_midgut <- deseq_res_midgut %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row)

sigASVs_vec_midgut <- sigASVs_midgut %>%
  pull(ASV)

midgut_filt <- prune_taxa(sigASVs_vec_midgut, midgut_final)

merged_results_midgut <- tax_table(midgut_final) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_midgut) %>%
  arrange(log2FoldChange) %>%
  filter(Genus != "g__Incertae_Sedis") %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

method_gear_midgut <- ggplot(merged_results_midgut) +
  geom_bar(aes(x = Genus, y = log2FoldChange), stat = "identity") +
  geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, 
                    ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Midgut",
       x = "Genus", y = "Log2 Fold Change") + 
  theme(panel.grid = element_blank())
method_gear_midgut

ggsave(file = "midgutdeseqmethodgear.svg", 
       plot = method_gear_midgut,
       height = 8, width = 20)


# Volcano plot
midgut_volcano_deseq <- deseq_res_midgut %>%
                        mutate(significant = ifelse(padj<0.01 & abs(log2FoldChange)>2, "Significant", "Not significant")) %>%
                        ggplot(aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
                        geom_point() +
                        geom_vline(aes(xintercept = 0)) +
                        scale_color_manual(name = "Legend", values = c("Significant" = "darkblue", "Not significant" = "lightblue")) +
                        theme(panel.grid = element_blank())
midgut_volcano_deseq

ggsave(filename="midgut_vol_plot.svg", midgut_volcano_deseq)

# Hindgut
hindgut_plus1 <- transform_sample_counts(hindgut_final, function(x) x + 1)
hindgut_deseq2 <- phyloseq_to_deseq2(hindgut_plus1, ~ method_gear)
hindgut_deseq <- DESeq(hindgut_deseq2)

deseq_res_hindgut <- results(hindgut_deseq, tidy = TRUE)
View(deseq_res_hindgut)

sigASVs_hindgut <- deseq_res_hindgut %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row)

sigASVs_vec_hindgut <- sigASVs_hindgut %>%
  pull(ASV)

hindgut_filt <- prune_taxa(sigASVs_vec_hindgut, hindgut_final)

merged_results_hindgut <- tax_table(hindgut_final) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_hindgut) %>%
  arrange(log2FoldChange) %>%
  filter(Genus != "g__Incertae_Sedis") %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

method_gear_hindgut_plot <- ggplot(merged_results_hindgut) +
  geom_bar(aes(x = Genus, y = log2FoldChange), stat = "identity") +
  geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, 
                    ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Hindgut",
       x = "Genus", y = "Log2 Fold Change") + 
  theme(panel.grid = element_blank())
method_gear_hindgut_plot

ggsave(file = "hindgutdeseqmethodgear.svg", 
       plot = method_gear_hindgut_plot,
       height = 8, width = 20)

# Volcano plot
hindgut_volcano_deseq <- deseq_res_hindgut %>%
  mutate(significant = ifelse(padj<0.01 & abs(log2FoldChange)>2, "Significant", "Not significant")) %>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
  geom_point() +
  geom_vline(aes(xintercept = 0)) +
  scale_color_manual(name = "Legend", values = c("Significant" = "darkblue", "Not significant" = "lightblue")) +
  theme(panel.grid = element_blank())
hindgut_volcano_deseq

ggsave(filename="hindgut_vol_plot.svg", hindgut_volcano_deseq)

# Gill
gill_plus1 <- transform_sample_counts(gill_phyloseq_object_clean, function(x) x + 1)
gill_deseq2 <- phyloseq_to_deseq2(gill_plus1, ~ method_gear)
gill_deseq <- DESeq(gill_deseq2)

deseq_res_gill <- results(gill_deseq, tidy = TRUE)
View(deseq_res_gill)

sigASVs_gill <- deseq_res_gill %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row)

sigASVs_vec_gill <- sigASVs_gill %>%
  pull(ASV)

gill_filt <- prune_taxa(sigASVs_vec_gill, gill_phyloseq_object_clean)

merged_results_gill <- tax_table(gill_phyloseq_object_clean) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_gill) %>%
  arrange(log2FoldChange) %>%
  filter(Genus != "g__Incertae_Sedis") %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

gill_method_gear <- ggplot(merged_results_gill) +
  geom_bar(aes(x = Genus, y = log2FoldChange), stat = "identity") +
  geom_errorbar(aes(x = Genus, ymin = log2FoldChange - lfcSE, 
                    ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Gill",
       x = "Genus", y = "Log2 Fold Change") + 
  theme(panel.grid = element_blank())
gill_method_gear

ggsave(file = "gilldeseqmethodgear.svg", 
       plot = gill_method_gear,
       height = 8, width = 20)

# Volcano plot
gill_volcano_deseq <- deseq_res_gill %>%
  mutate(significant = ifelse(padj<0.01 & abs(log2FoldChange)>2, "Significant", "Not significant")) %>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
  geom_point() +
  geom_vline(aes(xintercept = 0)) +
  scale_color_manual(name = "Legend", values = c("Significant" = "darkblue", "Not significant" = "lightblue")) +
  theme(panel.grid = element_blank())
gill_volcano_deseq


ggsave(filename="gill_vol_plot.svg", gill_volcano_deseq)

resultsNames(midgut_deseq)
resultsNames(hindgut_deseq)
resultsNames(gill_deseq)

# Making Table 3 (April 5th, 2026 by Michaela) 
# Hindgut summary
hindgut_deseq_summary_table <- data.frame(number_of_ASVs = nrow(deseq_res_hindgut),
                                          number_of_differentially_abundant_ASVs = nrow(subset(deseq_res_hindgut, padj < 0.01 & abs(log2FoldChange) > 2)),
                                          number_of_increased_ASVs = nrow(subset(deseq_res_hindgut, padj < 0.01 & abs(log2FoldChange) > 2 & log2FoldChange > 0)),
                                          number_of_decreased_ASVs = nrow(subset(deseq_res_hindgut, padj < 0.01 & abs(log2FoldChange) > 2 & log2FoldChange < 0))) 
hindgut_deseq_summary_table <- mutate(hindgut_deseq_summary_table, percent_changed_ASVs = ((number_of_differentially_abundant_ASVs) / (number_of_ASVs)) * 100)
hindgut_deseq_summary_table <- mutate(hindgut_deseq_summary_table, percent_increased_ASVs = ((number_of_increased_ASVs) / (number_of_differentially_abundant_ASVs)) * 100)
hindgut_deseq_summary_table <- mutate(hindgut_deseq_summary_table, percent_decreased_ASVs = ((number_of_decreased_ASVs) / (number_of_differentially_abundant_ASVs)) * 100)

# Midgut summary 
midgut_deseq_summary_table <- data.frame(number_of_ASVs = nrow(deseq_res_midgut),
                                          number_of_differentially_abundant_ASVs = nrow(subset(deseq_res_midgut, padj < 0.01 & abs(log2FoldChange) > 2)),
                                          number_of_increased_ASVs = nrow(subset(deseq_res_midgut, padj < 0.01 & abs(log2FoldChange) > 2 & log2FoldChange > 0)),
                                          number_of_decreased_ASVs = nrow(subset(deseq_res_midgut, padj < 0.01 & abs(log2FoldChange) > 2 & log2FoldChange < 0))) 
midgut_deseq_summary_table <- mutate(midgut_deseq_summary_table, percent_changed_ASVs = ((number_of_differentially_abundant_ASVs) / (number_of_ASVs)) * 100)
midgut_deseq_summary_table <- mutate(midgut_deseq_summary_table, percent_increased_ASVs = ((number_of_increased_ASVs) / (number_of_differentially_abundant_ASVs)) * 100)
midgut_deseq_summary_table <- mutate(midgut_deseq_summary_table, percent_decreased_ASVs = ((number_of_decreased_ASVs) / (number_of_differentially_abundant_ASVs)) * 100)

# Gill summary 
gill_deseq_summary_table <- data.frame(number_of_ASVs = nrow(deseq_res_gill),
                                         number_of_differentially_abundant_ASVs = nrow(subset(deseq_res_gill, padj < 0.01 & abs(log2FoldChange) > 2)),
                                         number_of_increased_ASVs = nrow(subset(deseq_res_gill, padj < 0.01 & abs(log2FoldChange) > 2 & log2FoldChange > 0)),
                                         number_of_decreased_ASVs = nrow(subset(deseq_res_gill, padj < 0.01 & abs(log2FoldChange) > 2 & log2FoldChange < 0))) 
gill_deseq_summary_table <- mutate(gill_deseq_summary_table, percent_changed_ASVs = ((number_of_differentially_abundant_ASVs) / (number_of_ASVs)) * 100)
gill_deseq_summary_table <- mutate(gill_deseq_summary_table, percent_increased_ASVs = ((number_of_increased_ASVs) / (number_of_differentially_abundant_ASVs)) * 100)
gill_deseq_summary_table <- mutate(gill_deseq_summary_table, percent_decreased_ASVs = ((number_of_decreased_ASVs) / (number_of_differentially_abundant_ASVs)) * 100)

# Combining summary tables into 1 table 
deseq_results_summary <- rbind(hindgut_deseq_summary_table, midgut_deseq_summary_table, gill_deseq_summary_table)
deseq_results_summary$Tissue <- c("Hindgut", "Midgut", "Gill")
deseq_results_summary <- select(deseq_results_summary, Tissue, everything())
deseq_results_summary <- rename (deseq_results_summary, 
                                 number_of_ASVs = "Total Number of ASVs",
                                 number_of_differentially_abundant_ASVs = "Number of Differentially Abundant ASVs",
                                 number_of_increased_ASVs = "Number of ASVs that Increased in Abundance",
                                 number_of_decreased_ASVs = "Number of ASVs that Decreased in Abundance", 
                                 percent_changed_ASVs = "% of Total ASVs that are Differentially Abundant",
                                 percent_increased_ASVs = "% of Differentially Abundant ASVs that Increased in Abundance", 
                                 percent_decreased_ASVs = "% of Differentially Abundant ASVs that Decreased in Abundance")

# Saving 
write.csv(deseq_results_summary, "deseqresults.csv", row.names = FALSE)