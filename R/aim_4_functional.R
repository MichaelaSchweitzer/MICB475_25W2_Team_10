#!/usr/bin/env Rscript

# PICRUSt2 analysis (aim 4)

# march 22nd, 2026, madi + phoebe 
# picrust analysis (R part, following module 19 tutorial) 
# initial generation of picrust object is done on the server 
# see PICRUSt-code.txt file 

library(ggpicrust2)
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)
library(ALDEx2)

#import metadata 
fish_metadata <- read_delim("fish_metadata.txt", delim="\t")

# import picrust2 objects 
arc_tree <- read.tree("picrust_out/arc.tre")
bac_tree <- read.tree("picrust_out/bac.tre")
arc <- read_delim('picrust_out/arc_marker_predicted_and_nsti.tsv.gz', delim="\t")
bac <- read_delim('picrust_out/bac_EC_predicted.tsv.gz', delim="\t")
comb <- read_delim('picrust_out/combined_marker_predicted_and_nsti.tsv.gz', delim="\t")

da <- ggpicrust2(data = bac, metadata = fish_metadata)

