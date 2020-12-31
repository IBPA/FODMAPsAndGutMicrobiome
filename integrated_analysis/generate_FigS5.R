#!/usr/bin/env Rscript
#' Generate Fig S5: A 3-genus microbiome biomarker for identifying the response
#'                  of IBS patients to low-FODMAP diet.

# Libraries
library(phyloseq)
library(biomformat)
library(ggplot2)
library(ggrepel)
library(microbiome)
library(stringr)
library(reshape2)
library(ggpubr)

# Shared code
source("lib/util.R")
source("lib/config_management.R")
source("lib/metadata_management.R")
source("lib/abundance_management.R")
source("lib/theme_util.R")

# A) Load data
df_metadata <- load_metadata()
df_metadata <- df_metadata[df_metadata$technology == "16S rRNA" &
                             !is.na(df_metadata$Response) &
                             !is.na(df_metadata$pre_diet_microbiome_sample_id),]
rownames(df_metadata) <- df_metadata$pre_diet_microbiome_sample_id

# B) Load Abundances
taxa_abundances <- load_taxa_abundances(INTEGRATED_IBS_DIR, df_metadata)
df_comb <- cbind(taxa_abundances$df_otu[rownames(df_metadata),], df_metadata)

# C) Create labels
response_freq <- plyr::count(df_comb$Response)
get_label <- function(response){
  return(sprintf("%s\n(n=%d)", response, 
                 response_freq[response_freq$x == response, c("freq")]))
}
response_labels <- get_label(levels(df_metadata$Response))
names(response_labels) <- levels(df_metadata$Response)

# D) Prepare for p-value comparison
my_comparisons <- list( c("High", "No"),
                        c("High", "Low"),
                        c("Low", "No"))

# E) Create the biomarker
df_comb$Genus_Biomarker <- df_comb$`Ruminococcaceae UCG-002` +
                           df_comb$`Ruminococcus 1` +
                           df_comb$Anaerostipes

# F) Plot and save
gPlot_Genus_Biomarker <- ggplot(df_comb, aes(x=Response, 
                                             y=Genus_Biomarker, 
                                             fill=Response))+
  stat_compare_means(comparisons = my_comparisons, size=4, method = "t.test",
                     method.args = list(alternative = "greater"))+
  scale_x_discrete("Low-FODMAP diet response", labels = response_labels)+
  scale_y_continuous("Sum of CLR-transformed abundances for\n genuses Ruminococcaceae UCG-002,\n Ruminococcus 1 and Anaerostipes")+
  geom_boxplot(outlier.shape = 21)+
  scale_fill_manual(values=c("No"=unname(colors_assigned["D"]),
                             "Low"=unname(colors_assigned["B"]),
                             "High"=unname(colors_assigned["A"])))+
  my_base_theme %+replace% theme(legend.position = "none")
print(gPlot_Genus_Biomarker)

lib.util.ggsave(file.path("results", "FigS5.png"),
                gPlot_Genus_Biomarker, dpi = 400)
