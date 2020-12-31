#!/usr/bin/env Rscript
#' Generate Fig S3: Classification performance based on fatty acid metabolism 
#'                  pathway enrichment of 16S rRNA data.
#' 1) Predict High vs. No response
#' 2) Predict High vs. Low response

# Library
library(phyloseq)
library(biomformat)
library(ggplot2)
library(ggrepel)
library(microbiome)
library(stringr)
library(reshape2)
library(nsprcomp)
library(PRROC)
library(ggpubr)

# Shared code
source("lib/util.R")
source("lib/config_management.R")
source("lib/metadata_management.R")
source("lib/abundance_management.R")
source("lib/classification.R")
source("lib/theme_util.R")
my_base_theme <- get_base_theme(0.9)

# A) Load data
df_metadata <- load_metadata()
df_metadata <- load_metadata()
df_metadata <- df_metadata[df_metadata$technology == "16S rRNA" &
                             !is.na(df_metadata$Response) &
                             !is.na(df_metadata$pre_diet_microbiome_sample_id),]
rownames(df_metadata) <- df_metadata$pre_diet_microbiome_sample_id
pathway_abundances <- load_pathway_abundances(INTEGRATED_IBS_DIR, df_metadata)
df_comb <- cbind(pathway_abundances$df_otu[rownames(df_metadata),], df_metadata)

# B) Predict High vs. No response
RF_MaxNodes <- 2
df_res <- s_classify_binary_RF_loo_cr(df_comb, x_cols = c("Fatty.acid.metabolism") , 
                                      label_col = "Response",
                                      label_group = list(true=c("High"), false=c("No")))
gPlot_roc_High_No <- get_roc_plot(df_res, color_code = "D", add_star = FALSE)
gPlot_pr_High_No <- get_pr_plot(df_res, color_code = "D", add_star = FALSE)

# C) Predict High vs. Low response
df_res <- s_classify_binary_RF_loo_cr(df_comb, x_cols = c("Fatty.acid.metabolism") , 
                                      label_col = "Response",
                                      label_group = list(true=c("High"), false=c("Low")))
gPlot_roc_High_Low <- get_roc_plot(df_res, color_code = "B", add_star = FALSE)
gPlot_pr_High_Low <- get_pr_plot(df_res, color_code = "B", add_star = FALSE)

# D) Plot and save
gPlot <- cowplot::plot_grid(gPlot_roc_High_No + ggtitle("High vs. No response"), 
                            gPlot_roc_High_Low + ggtitle("High vs. Low response"), 
                            gPlot_pr_High_No, 
                            gPlot_pr_High_Low,
                            nrow = 2, ncol = 2, labels = c("A", "C", "B", "D"))
print(gPlot)
lib.util.load_fonts()
lib.util.ggsave(file.path("results", "FigS3.png"),
                gPlot, width = 5, height = 5, dpi=400)
