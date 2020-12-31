#!/usr/bin/env Rscript
#' Generate Fig S6: Differential pathway analysis of pre-treatment microbiome 
#'                  for individual studies with 16S rRNA data

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
my_base_theme <- get_base_theme(0.8)

# A) Load data
df_metadata <- load_metadata()
df_metadata <- df_metadata[df_metadata$technology == "16S rRNA" &
                             !is.na(df_metadata$Response) &
                             !is.na(df_metadata$pre_diet_microbiome_sample_id),]
rownames(df_metadata) <- df_metadata$pre_diet_microbiome_sample_id

# B) Load Abundances
pathway_abundances <- load_pathway_abundances(INTEGRATED_IBS_DIR, df_metadata)
df_comb <- cbind(pathway_abundances$df_otu[rownames(df_metadata),], df_metadata)

# Prepare for p-value comparison
my_comparisons <- list( c("High", "No"),
                        c("High", "Low"),
                        c("Low", "No"))

# C) Methane metabolism
gPlot_methane <- ggplot(df_comb, aes(x=Response, 
                                     y=Methane.metabolism, 
                                     fill=Response))+
  facet_grid(cols = vars(reference))+
  stat_compare_means(comparisons = my_comparisons, size=2,
                     label.y = c(2.7, 2.8, 2.88),
                     method.args = list(alternative = "greater"))+
  scale_x_discrete("Low-FODMAP diet response")+
  scale_y_continuous("CLR-Transformed abundance of\nmicrobiome methane metabolism")+
  expand_limits(y=2.9)+
  geom_boxplot(outlier.shape = 21)+
  scale_fill_manual(values=c("No"=unname(colors_assigned["D"]),
                             "Low"=unname(colors_assigned["B"]),
                             "High"=unname(colors_assigned["A"])))+
  my_base_theme %+replace% theme(legend.position = "none")
print(gPlot_methane)

# D) Fatty acid metabolism
gPlot_fatty_acid <- ggplot(df_comb, aes(x=Response, 
                                        y=Fatty.acid.metabolism, 
                                        fill=Response))+
  facet_grid(cols = vars(reference))+
  stat_compare_means(comparisons = my_comparisons, size=2,
                     label.y = c(1, 1.12, 1.2),
                     method.args = list(alternative = "greater"))+
  scale_x_discrete("Low-FODMAP diet response")+
  scale_y_continuous("CLR-Transformed abundance of\nmicrobiome fatty acid metabolism")+
  geom_boxplot(outlier.shape = 21)+
  scale_fill_manual(values=c("No"=unname(colors_assigned["D"]),
                             "Low"=unname(colors_assigned["B"]),
                             "High"=unname(colors_assigned["A"])))+
  my_base_theme %+replace% theme(legend.position = "none")
print(gPlot_fatty_acid)

# E) Plot all
gPlot <- cowplot::plot_grid(gPlot_methane, gPlot_fatty_acid,
                            ncol = 1, nrow = 2, labels = c("A", "B"),
                            align = "hv", axis = "tblr")
print(gPlot)

lib.util.ggsave(file.path("results", "FigS6.png"),
              gPlot, dpi = 400, width = 5.5, height = 4)
