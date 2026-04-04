install.packages("UpSetR")
install.packages("ComplexHeatmap")
library(grid)
library(UpSetR)
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

load("all_phyloseq_objects_bacteriaonly_pruned.RData")


#Get the relative abundance of each body site
gill_RA <- transform_sample_counts(gill_final, fun=function(x) x/sum(x))
hindgut_RA <- transform_sample_counts(hindgut_final, fun=function(x) x/sum(x))
midgut_RA <- transform_sample_counts(midgut_final, fun=function(x) x/sum(x))
skin_RA <- transform_sample_counts(skin_final, fun=function(x) x/sum(x))


###All the options (just for my reference) 
# benthic trawl 
# gillnet 
# midwater trawl 
# net 
# rod and reel 
# spear


### Break it up by method_gear ###
# Gill 

gill_bt <- subset_samples(gill_RA, method_gear=="benthic trawl")
gill_gnet <- subset_samples(gill_RA, method_gear=="gillnet")
gill_mt <- subset_samples(gill_RA, method_gear=="midwater trawl")
gill_net <- subset_samples(gill_RA, method_gear=="net")
gill_rr <- subset_samples(gill_RA, method_gear=="rod and reel")
gill_sp <- subset_samples(gill_RA, method_gear=="spear")

#Checking the number of samples in each category
nsamples(gill_bt)
nsamples(gill_gnet)
nsamples(gill_mt)
nsamples(gill_net)
nsamples(gill_rr)
nsamples(gill_sp)

# Hindgutgill_mt# Hindgut
hindgut_bt <- subset_samples(hindgut_RA, method_gear=="benthic trawl")
hindgut_gnet <- subset_samples(hindgut_RA, method_gear=="gillnet")
hindgut_mt <- subset_samples(hindgut_RA, method_gear=="midwater trawl")
hindgut_net <- subset_samples(hindgut_RA, method_gear=="net")
hindgut_rr <- subset_samples(hindgut_RA, method_gear=="rod and reel")
hindgut_sp <- subset_samples(hindgut_RA, method_gear=="spear")

#Checking the number of samples in each category
nsamples(hindgut_bt)
nsamples(hindgut_gnet)
nsamples(hindgut_mt)
nsamples(hindgut_net)
nsamples(hindgut_rr)
nsamples(hindgut_sp)


# Midgut 
midgut_bt <- subset_samples(midgut_RA, method_gear=="benthic trawl")
midgut_gnet <- subset_samples(midgut_RA, method_gear=="gillnet")
midgut_mt <- subset_samples(midgut_RA, method_gear=="midwater trawl")
midgut_net <- subset_samples(midgut_RA, method_gear=="net")
midgut_rr <- subset_samples(midgut_RA, method_gear=="rod and reel")
midgut_sp <- subset_samples(midgut_RA, method_gear=="spear")

#Checking the number of samples in each category
nsamples(midgut_bt)
nsamples(midgut_gnet)
nsamples(midgut_mt)
nsamples(midgut_net)
nsamples(midgut_rr)
nsamples(midgut_sp)

# Skin (I know we aren't planning on using it but thought I'd have it here in case we change our mind)
skin_bt <- subset_samples(skin_RA, method_gear=="benthic trawl")
skin_gnet <- subset_samples(skin_RA, method_gear=="gillnet")
skin_mt <- subset_samples(skin_RA, method_gear=="midwater trawl")
skin_net <- subset_samples(skin_RA, method_gear=="net")
skin_rr <- subset_samples(skin_RA, method_gear=="rod and reel")
skin_sp <- subset_samples(skin_RA, method_gear=="spear")

#Checking the number of samples in each category
nsamples(skin_bt)
nsamples(skin_gnet)
nsamples(skin_mt)
nsamples(skin_net)
nsamples(skin_rr)
nsamples(skin_sp)

### Get the core members ###

# Gill
gill_bt_core <- core_members(gill_bt, detection=0, prevalence=0.5)
gill_gnet_core <- core_members(gill_gnet, detection=0, prevalence=0.5)
gill_mt_core <- core_members(gill_mt, detection=0, prevalence=0.5)
gill_net_core <- core_members(gill_net, detection=0, prevalence=0.5)
gill_rr_core <- core_members(gill_rr, detection=0, prevalence=0.5)
gill_sp_core <- core_members(gill_sp, detection=0, prevalence=0.5)

# Hindgut
hindgut_bt_core <- core_members(hindgut_bt, detection=0, prevalence=0.5)
hindgut_gnet_core <- core_members(hindgut_gnet, detection=0, prevalence=0.5)
hindgut_mt_core <- core_members(hindgut_mt, detection=0, prevalence=0.5)
hindgut_net_core <- core_members(hindgut_net, detection=0, prevalence=0.5)
hindgut_rr_core <- core_members(hindgut_rr, detection=0, prevalence=0.5)
hindgut_sp_core <- core_members(hindgut_sp, detection=0, prevalence=0.5)

# Midgut
midgut_bt_core <- core_members(midgut_bt, detection=0, prevalence=0.5)
midgut_gnet_core <- core_members(midgut_gnet, detection=0, prevalence=0.5)
midgut_mt_core <- core_members(midgut_mt, detection=0, prevalence=0.5)
midgut_net_core <- core_members(midgut_net, detection=0, prevalence=0.5)
midgut_rr_core <- core_members(midgut_rr, detection=0, prevalence=0.5)
midgut_sp_core <- core_members(midgut_sp, detection=0, prevalence=0.5)

# Skin 
skin_bt_core <- core_members(skin_bt, detection=0, prevalence=0.5)
skin_gnet_core <- core_members(skin_gnet, detection=0, prevalence=0.5)
skin_mt_core <- core_members(skin_mt, detection=0, prevalence=0.5)
skin_net_core <- core_members(skin_net, detection=0, prevalence=0.5)
skin_rr_core <- core_members(skin_rr, detection=0, prevalence=0.5)
skin_sp_core <- core_members(skin_sp, detection=0, prevalence=0.5)

### Visualizing some of these with a table (Can delete) ###
prune_taxa(midgut_gnet_core,midgut_final) %>%
  tax_table()

prune_taxa(gill_bt_core, gill_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`method_gear`, scales ="free")

prune_taxa(midgut_gnet_core, midgut_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`method_gear`, scales ="free")

### RELAXING PARAMETERS ### 
# This is needed since even with prevalence of 50% we lose so much and our Venn Diagram won't be very informative 

gill_bt_core_relax <- core_members(gill_bt, detection=0.001, prevalence=0.1)
gill_gnet_core_relax <- core_members(gill_gnet, detection=0.001, prevalence=0.1)
gill_mt_core_relax <- core_members(gill_mt, detection=0.001, prevalence=0.1)
gill_net_core_relax <- core_members(gill_net, detection=0.001, prevalence=0.1)
gill_rr_core_relax <- core_members(gill_rr, detection=0.001, prevalence=0.1)
gill_sp_core_relax <- core_members(gill_sp, detection=0.001, prevalence=0.1)

# Hindgut
hindgut_bt_core_relax <- core_members(hindgut_bt, detection=0.001, prevalence=0.1)
hindgut_gnet_core_relax <- core_members(hindgut_gnet, detection=0.001, prevalence=0.1)
hindgut_mt_core_relax <- core_members(hindgut_mt, detection=0.001, prevalence=0.1)
hindgut_net_core_relax <- core_members(hindgut_net, detection=0.001, prevalence=0.1)
hindgut_rr_core_relax <- core_members(hindgut_rr, detection=0.001, prevalence=0.1)
hindgut_sp_core_relax <- core_members(hindgut_sp, detection=0.001, prevalence=0.1)

# Midgut
midgut_bt_core_relax <- core_members(midgut_bt, detection=0.001, prevalence=0.1)
midgut_gnet_core_relax <- core_members(midgut_gnet, detection=0.001, prevalence=0.1)
midgut_mt_core_relax <- core_members(midgut_mt, detection=0.001, prevalence=0.1)
midgut_net_core_relax <- core_members(midgut_net, detection=0.001, prevalence=0.1)
midgut_rr_core_relax <- core_members(midgut_rr, detection=0.001, prevalence=0.1)
midgut_sp_core_relax <- core_members(midgut_sp, detection=0.001, prevalence=0.1)

# Skin 
skin_bt_core_relax <- core_members(skin_bt, detection=0.001, prevalence=0.1)
skin_gnet_core_relax <- core_members(skin_gnet, detection=0.001, prevalence=0.1)
skin_mt_core_relax <- core_members(skin_mt, detection=0.001, prevalence=0.1)
skin_net_core_relax <- core_members(skin_net, detection=0.001, prevalence=0.1)
skin_rr_core_relax <- core_members(skin_rr, detection=0.001, prevalence=0.1)
skin_sp_core_relax <- core_members(skin_sp, detection=0.001, prevalence=0.1)

## Make this all into a list ##

gill_core_list <- list("Benthic Trawl" = gill_bt_core_relax, "Gill Net" = gill_gnet_core_relax, 
                       "Midwater Trawl" = gill_mt_core_relax, "Net" = gill_net_core_relax,
                       "Rod and Reel" = gill_rr_core_relax, "Spear" = gill_sp_core_relax)

midgut_core_list <- list("Benthic Trawl" = midgut_bt_core_relax, "Gill Net" = midgut_gnet_core_relax, 
                         "Midwater Trawl" = midgut_mt_core_relax, "Net" = midgut_net_core_relax,
                         "Rod and Reel" = midgut_rr_core_relax, "Spear" = midgut_sp_core_relax)

hindgut_core_list <- list("Benthic Trawl" = hindgut_bt_core_relax, "Gill Net" = hindgut_gnet_core_relax, 
                          "Midwater Trawl" = hindgut_mt_core_relax, "Net" = hindgut_net_core_relax,
                          "Rod and Reel" = hindgut_rr_core_relax, "Spear" = hindgut_sp_core_relax)

skin_core_list <- list("Benthic Trawl" = skin_bt_core_relax, "Gill Net" = skin_gnet_core_relax, 
                       "Midwater Trawl" = skin_mt_core_relax, "Net" = skin_net_core_relax,
                       "Rod and Reel" = skin_rr_core_relax, "Spear" = skin_sp_core_relax)

### Installing SF package MAKE SURE TO TYPE "no" WHEN IT PROMPTS YOU
# install.packages("sf")
# library(sf)


skin_venn <- ggVennDiagram(x = skin_core_list)
skin_venn

### The Venn diagram looks really bad (!) since theres 6 groups. Going to make an UpSet plot instead which is 
# much cleaner and shows the same thing
 
# For now I set "nintersects" to 64 (the max possible number) so that all intersects are shown.
# this can definitely be changed

### If you want to remake these images, just uncomment (ctrl + shift + C)


# # Gill upset plot
# 
# png(filename = "gill_upset.png", width = 1800, height = 1800, res = 300)
# 
# 
# upset(fromList(gill_core_list), 
#       nsets = 6,
#       order.by = "freq", 
#       keep.order = TRUE,
#       sets = c("Benthic Trawl", "Net", "Spear", "Gill Net", "Midwater Trawl", "Rod and Reel"),
#       nintersects = 64)
# 
# dev.off()
# # Gill upset plot for powerpoint
# png(filename = "gill_upset_ppt.png", width = 2500, height = 1800, res = 300)
# 
# 
# upset(fromList(gill_core_list), 
#       nsets = 6,
#       order.by = "freq", 
#       keep.order = TRUE,
#       sets = c("Benthic Trawl", "Net", "Spear", "Gill Net", "Midwater Trawl", "Rod and Reel"),
#       nintersects = 64)
# 
# dev.off()
# 
# # Hindgut Upset Plot
# png(filename = "hindgut_upset.png", width = 1800, height = 1800, res = 300)
# 
# upset(fromList(hindgut_core_list), 
#       nsets = 6,
#       order.by = "freq", 
#       keep.order = FALSE,
#       sets = c("Benthic Trawl", "Net", "Spear", "Gill Net", "Midwater Trawl", "Rod and Reel"),
#       nintersects = 64)
# 
# dev.off()
# 
# # Hindgut upset plot for powerpoint
# png(filename = "hindgut_upset_ppt.png", width = 2500, height = 1800, res = 300)
# 
# upset(fromList(hindgut_core_list), 
#       nsets = 6,
#       order.by = "freq", 
#       keep.order = FALSE,
#       sets = c("Benthic Trawl", "Net", "Spear", "Gill Net", "Midwater Trawl", "Rod and Reel"),
#       nintersects = 64)
# 
# dev.off()
# 
# 
# # Midgut Upset Plot
# png(filename = "midgut_upset.png", width = 1800, height = 1800, res = 300)
# 
# upset(fromList(midgut_core_list), 
#       nsets = 6,
#       order.by = "freq", 
#       keep.order = FALSE,
#       sets = c("Benthic Trawl", "Net", "Spear", "Gill Net", "Midwater Trawl", "Rod and Reel"),
#       nintersects = 64)
# dev.off()
# 
# #Midgut upset plot for powerpoint
# png(filename = "midgut_upset_ppt.png", width = 2500, height = 1800, res = 300)
# 
# 
# upset(fromList(midgut_core_list), 
#       nsets = 6,
#       order.by = "freq", 
#       keep.order = FALSE,
#       sets = c("Benthic Trawl", "Net", "Spear", "Gill Net", "Midwater Trawl", "Rod and Reel"),
#       nintersects = 64)
# 
# dev.off()
# 
# 
# # Skin Upset plot
# png(filename = "skin_upset.png", width = 1800, height = 1800, res = 300)
# 
# upset(fromList(skin_core_list), 
#       nsets = 6,
#       order.by = "freq", 
#       keep.order = FALSE,
#       sets = c("Benthic Trawl", "Net", "Spear", "Gill Net", "Midwater Trawl", "Rod and Reel"),
#       nintersects = 64)
# dev.off()
# 
# #Skin upset plot for ppt
# png(filename = "skin_upset_ppt.png", width = 2500, height = 1800, res = 300)
# 
# 
# upset(fromList(skin_core_list), 
#       nsets = 6,
#       order.by = "freq", 
#       keep.order = FALSE,
#       sets = c("Benthic Trawl", "Net", "Spear", "Gill Net", "Midwater Trawl", "Rod and Reel"),
#       nintersects = 64)
# 
# dev.off()
# 
# 
# 
# 










