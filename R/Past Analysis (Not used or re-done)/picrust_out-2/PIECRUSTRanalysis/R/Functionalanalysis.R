---
# title: "Group 10 - ggpicrust2"
# output: html_document
# Adapted from https://github.com/armetcal/MICB_305/blob/main/R%20Scripts/Module%2015%20-%20ggpicrust2.Rmd
---
#Phoebe McNair-Luxon - March 23rd 2025
#In this R script, we will run the ggpicrust2 pipeline on our KO term data. Specifically, we'll took at ___

#ggpicrust2 has several dependencies that aren't automatically installed, so keep an eye on the code output as it'll throw an error and may or may not break the code. Errors about missing packages mean you need to install them first.
#ggpicrust2 is a fantastic tool when it works properly. Unfortunately, it's very buggy. This script will show how to theoretically run the entire pipeline in one command, but as of summer 2025, the package needs updating. We will therefore run the pipeline piece by piece instead in order to make debugging easier.

 #### {r Load data and libraries}

install.packages("devtools") #ran
devtools::install_github("cafferychen777/ggpicrust2") #ran, updated all packages
install.packages('MicrobiomeStat') #run
BiocManager::install("KEGGREST") #run
install.packages("GGally") #run

library(tidyverse)
library(phyloseq)
library(ggpicrust2)

# Load data
# Metadata - I'm loading in from the previously constructed phyloseq object as the data is already formatted!

meta_gill_fromphyloseq = gill_phyloseq_object %>% .@sam_data %>% data.frame() %>%
  rownames_to_column('sample_name')

meta_hindgut_fromphyloseq = hindgut_phyloseq_object %>% .@sam_data %>% data.frame() %>%
  rownames_to_column('sample_name')

meta_skin_fromphyloseq = skin_phyloseq_object %>% .@sam_data %>% data.frame() %>%
  rownames_to_column('sample_name')

meta_midgut_fromphyloseq = midgut_phyloseq_object %>% .@sam_data %>% data.frame() %>%
  rownames_to_column('sample_name')

KO_data = read.delim('KO_metagenome_out/pred_metagenome_unstrat/pred_metagenome_unstrat.tsv')

KO_edit = read.delim('KO_metagenome_out/pred_metagenome_unstrat/pred_metagenome_unstrat.tsv', row.names = 1)

####HINDGUT AND MIDGUT ->
###From video KO linear differential abundance analysis
meta_hindgut_midgut = rbind(meta_hindgut_fromphyloseq, meta_midgut_fromphyloseq)

#remove X from abundance column names
colnames(KO_edit) <- gsub("^X", "", colnames(KO_edit))
intersect(colnames(KO_edit), meta_hindgut_midgut$sample_name)

#filter KO data to only hindgut and midgut ->
KO_edit_hindgut_midgut <- KO_edit[, meta_hindgut_midgut$sample_name]

#check matches (continue if TRUE)
all(colnames(KO_edit_hindgut_midgut) == meta_hindgut_midgut$sample_name)

# run linear differential abundance
daa_results_hindgut_midgut = pathway_daa(abundance = KO_edit_hindgut_midgut, #can run this on all the different inputs (there should be 3)
                    metadata = meta_hindgut_midgut,
                    group = "sample_type", #what we are interested in
                    daa_method = "LinDA") #Linear differential abundance

pathway_annotation(pathway = "KO",
                   daa_results_df = daa_results_hindgut_midgut,
                   ko_to_kegg = TRUE) #transform individual KO terms in to Kegg pathways, can only set this to true if using KO terms

#Save KEGG pathway results as RDF
saveRDS(daa_annotated_results_df,'KO_metagenome_out/pred_metagenome_unstrat/daa_results_hindgut_midgut.rds')

#plot KEGG pathways dimension and aesthetics attempt 1
p_midgut_hindgut = pathway_errorbar(abundance = ko_filt,
                     daa_results_df = daa_annotated_results_df,
                     Group = lpalmvstongue$body.site,
                     p_values_threshold = 5e-11,
                     order = "pathway_class",
                     ko_to_kegg = T,
                     p_value_bar = TRUE,
                     x_lab = "pathway_name")
p_midgut_hindgut

#plot KEGG pathways dimension and aesthetics attempt 2 (debugged)
source('ggpicrust2_errorbar_function_fixed.R')

daa_annotated_results_df %>% filter(p_adjust< 5e-11, abs(log2FoldChange) > 5) %>% nrow()

peb = pathway_errorbar_fixed(abundance = ko_filt,
                             daa_results_df = daa_annotated_results_df,
                             Group = lpalmvstongue$body.site,
                             wrap_label = T, wraplength=60,
                             fc_cutoff = 5, order_by_log = F,
                             p_values_threshold = 5e-11,
                             order = "pathway_class",
                             ko_to_kegg = TRUE,
                             p_value_bar = TRUE,
                             x_lab = "pathway_name")
peb
ggsave("../Results/Plots/Module 15 - KEGG Error Bar Fixed.png",
       plot = peb, width = 18, height = 10)


####HINDGUT AND GILL ->
###From video KO linear differential abundance analysis
meta_hindgut_midgut = rbind(meta_hindgut_fromphyloseq, meta_gill_fromphyloseq)

#remove X from abundance column names
colnames(KO_edit_h_g) <- gsub("^X", "", colnames(KO_edit))
intersect(colnames(KO_edit_h_g), meta_hindgut_gill$sample_name)

#filter KO data to only hindgut and gill ->
KO_edit_hindgut_midgut <- KO_edit_h_g[, meta_hindgut_gill$sample_name]

#check matches (continue if TRUE)
all(colnames(KO_edit_hindgut_gill) == meta_hindgut_gill$sample_name)

# run linear differential abundance
daa_results_hindgut_gill = pathway_daa(abundance = KO_edit_hindgut_gill, #can run this on all the different inputs (there should be 3)
                                         metadata = meta_hindgut_gill,
                                         group = "sample_type", #what we are interested in
                                         daa_method = "LinDA") #Linear differential abundance

pathway_annotation(pathway = "KO",
                   daa_results_df = daa_results_hindgut_gill,
                   ko_to_kegg = TRUE) #transform individual KO terms in to Kegg pathways, can only set this to true if using KO terms




###NOT MY CODE FROM HERE ON ->
# ATTEMPT 1 - Whole Pipeline

# The ggpicrust2 library comes with a main function that will perform the entire pipeline in a single command.
# The ggpicrust2 command is designed for relatively straightforward analyses. For example, the below code only works when the variable of interest is numeric or has only two levels. If you need to do a more complex analysis, you can run the individual steps of the pipeline separately. See the ggpicrust2 documentation for more details.

# ```{r run ggpicrust2 on KO terms}
# Filter to only compare two body sites
lpalmvstongue = meta %>% filter(body.site %in% c("left palm", "tongue"))
ko_filt = ko %>% select(`function.`,all_of(lpalmvstongue$sample_name))

# This is throwing an error - see if you can use debugonce() to find it.
debugonce(ggpicrust2)
ggpicrust2(data = ko_filt,
           metadata = lpalmvstongue,
           group = "body.site",
           pathway = "KO",
           daa_method = "LinDA",
           ko_to_kegg = TRUE,
           order = "pathway_class",
           p_values_bar = TRUE,
           x_lab = "pathway_name")
```

# ATTEMPT 2 - Run each step separately

As you can see above, the error bar function throws an error. Unfortunately, the ggpicrust2 command doesn't let us skip that step.

```{r Differential analysis}
# First column should be row names when running pipeline manually, so I'm just reloading it as such
ko = read.delim('../Datasets/pred_metagenome_unstrat.tsv',row.names = 1)

# Filter to only compare two body sites
lpalmvstongue = meta %>% filter(body.site %in% c("left palm", "tongue"))
ko_filt = ko %>% select(all_of(lpalmvstongue$sample_name))

# Perform pathway differential abundance analysis (DAA) using LinDA method
daa_results_df = pathway_daa(abundance = ko_filt,
                             metadata = lpalmvstongue,
                             group = "body.site",
                             daa_method = "LinDA",
                             select = NULL, reference = NULL)

# Annotate pathway results using KO to KEGG conversion
# This will take a while to run as it downloads KEGG data - save as an .rds file if you want to avoid re-running it
# daa_annotated_results_df = pathway_annotation(pathway = "KO",
#                                               daa_results_df = daa_results_df,
#                                               ko_to_kegg = TRUE)
# saveRDS(daa_annotated_results_df,'../Datasets/module15_daa_annotated_results_df.rds')
daa_annotated_results_df = readRDS('../Datasets/module15_daa_annotated_results_df.rds')

# Generate pathway error bar plot
# This function is very buggy and may not work properly
p = pathway_errorbar(abundance = ko_filt,
                     daa_results_df = daa_annotated_results_df,
                     Group = lpalmvstongue$body.site,
                     p_values_threshold = 5e-11,
                     order = "pathway_class",
                     ko_to_kegg = T,
                     p_value_bar = TRUE,
                     x_lab = "pathway_name")
p
```

The plot is kind of messed up. The P values aren't visible, for example. If you have only one significant pathway, it might also throw an error.

Try using this function instead, which I have edited to work a bit better.

Feel free to use this function, edit the function to make it even better, or make your own simpler plot to visualize your data. The differential abundance plot options would be perfectly acceptable.

```{r fixed function}
source('ggpicrust2_errorbar_function_fixed.R')

daa_annotated_results_df %>% filter(p_adjust< 5e-11, abs(log2FoldChange) > 5) %>% nrow()

peb = pathway_errorbar_fixed(abundance = ko_filt,
                     daa_results_df = daa_annotated_results_df,
                     Group = lpalmvstongue$body.site,
                     wrap_label = T, wraplength=60,
                     fc_cutoff = 5, order_by_log = F,
                     p_values_threshold = 5e-11,
                     order = "pathway_class",
                     ko_to_kegg = TRUE,
                     p_value_bar = TRUE,
                     x_lab = "pathway_name")
peb
ggsave("../Results/Plots/Module 15 - KEGG Error Bar Fixed.png",
       plot = peb, width = 18, height = 10)
```

Let's see which features are essentially the same (they'll be highly correlated). Keep this in mind when writing up your results.

```{r identical features}
sig_features = daa_annotated_results_df %>%
  filter(p_adjust < 5e-11, abs(log2FoldChange) > 5) %>%
  pull('feature') %>% unique()

ko_relab = read.delim('../Datasets/pred_metagenome_unstrat.tsv',row.names = 1) %>%
  apply(2,function(x) x/sum(x)) %>% as.data.frame()
colSums(ko_relab) # Should be all 1

stats = ko_relab %>% t() %>%  # switch rows and columns
  as.data.frame() %>% # Lets us use select()
  select(all_of(sig_features)) %>%
  cor(method = 'pearson') #Pearson tests between every feature pair

# Define the color palette
color_palette <- colorRampPalette(c("blue", "white", "red"))(40)

# Define breaks for the color scale
breaks <- seq(-1, 1, length.out = 41)  # 40 colors + 1 for the endpoint

# Create the clustered heatmap with centered color at zero
# install.packages("pheatmap")
library('pheatmap')
pheatmap(stats,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = color_palette,
         breaks = breaks,
         main = '',
         fontsize_row = 10,
         fontsize_col = 10,
         filename = "../Results/Plots/Module 15 - KEGG Heatmap.png",
         height= 9, width = 9)
```

See how similar your sample groups are. Unfortunately, the size of the text can't be increased, even with theme() commands.

                                                   ```{r PCA}
                                                   # Generate pathway PCA plot
                                                   pathway_pca(abundance = ko_filt,
                                                               metadata = lpalmvstongue,
                                                               group = "body.site")
