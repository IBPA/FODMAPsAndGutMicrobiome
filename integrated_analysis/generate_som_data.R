#!/usr/bin/env Rscript
#' Generate supplementary online data relating to pre-diet microbiome (16S rRNA)

# Libraries
library(phyloseq)
library(biomformat)
library(ggplot2)
library(ggrepel)
library(microbiome)
library(stringr)
library(reshape2)
library(xlsx)

# Shared code
source("lib/util.R")
source("lib/config_management.R")
source("lib/metadata_management.R")
source("lib/abundance_management.R")
source("lib/theme_util.R")

# 1) Load metadata from all IBS patients
df_metadata <- load_metadata()
df_metadata <- df_metadata[df_metadata$technology == "16S rRNA" &
                           !is.na(df_metadata$Response) &
                           !is.na(df_metadata$pre_diet_microbiome_sample_id),]
rownames(df_metadata) <- df_metadata$pre_diet_microbiome_sample_id
df_metadata <- df_metadata[, c("reference","Response")]

# 2) Load PICRUSt Data
pathway_abundances <- load_pathway_abundances(INTEGRATED_IBS_DIR, df_metadata)
colnames(pathway_abundances$df_otu) <- pathway_abundances$friendly_names[colnames(pathway_abundances$df_otu)]

# 3) Load Taxa Data
taxa_abundances <- load_taxa_abundances(INTEGRATED_IBS_DIR, df_metadata)

# 4) Save in excel
filename <- file.path("results", "Data.xlsx")
write.xlsx(df_metadata, file=filename, sheetName="Metadata", row.names=TRUE)
write.xlsx(pathway_abundances$df_otu, file=filename, sheetName="KEGG Pathways (L3)", append=TRUE, row.names=TRUE)
write.xlsx(taxa_abundances$df_otu, file=filename, sheetName="Taxa (Genus)", append=TRUE, row.names=TRUE)
