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
colnames(fish_metadata)[1] <- "sample_name"

#trying to fix metadata so that it imports correctly 
#make sample IDs rownames instead of separate column 
metadata_formatted <- fish_metadata %>% 
  column_to_rownames("sample_name")

# import picrust2 objects 
x <- read_delim('picrust_out-2/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz', delim="\t")



da <- ggpicrust2(data = x, metadata = metadata_formatted, group = "sample_type")

# i think problem is different # of samples between abundance table and metadata
#i don't know why not all samples are in the abundance data though 
