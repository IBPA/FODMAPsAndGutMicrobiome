#!/usr/bin/env Rscript

# Libraries
library(phyloseq)
library(biomformat)
library(microbiome)
library(stringr)
library(reshape2)

source("lib/abundance_management.R")

src_dirs <- c("./preprocessing/s1/microbiome_data/", 
              "./preprocessing/s2/microbiome_data/", 
              "./preprocessing/s3/microbiome_data/",
              "./preprocessing/s6/microbiome_data/")
ps <- integrate_data_taxa(src_dirs)
out_filename <- file.path("./integrated_analysis/integrated_all/DADA2.rds")
saveRDS(ps, out_filename)
print(sprintf("Saved to %s", out_filename))
