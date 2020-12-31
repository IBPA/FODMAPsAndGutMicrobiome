#!/usr/bin/env Rscript
#' Generate Fig2: Pre-diet microbial differential abundances for IBS patients 
#'                with high versus no response to the low-FODMAP diet
#' 1) Get distribution IBS-SSS scores before and after low-FODMAP diet.
#' 2) Perform microbiome differential abundance analysis for high vs. no response
#'    IBS patients and plot.
#' 3) Plot above together.

# Libraries
library(phyloseq)
library(biomformat)
library(ggplot2)
library(ggrepel)
library(microbiome)
library(stringr)
library(reshape2)

# Shared code
source("lib/util.R")
source("lib/config_management.R")
source("lib/metadata_management.R")
source("lib/abundance_management.R")
source("lib/theme_util.R")

###############################################################################
############## 1. Plot distribution of pre/post IBS-SSS scores ################
###############################################################################
# 1.A) Load metadata from all IBS patients
df_metadata <- load_metadata()
df_metadata <- df_metadata[!is.na(df_metadata$Response),]

# 1.B) Scatter-plot pre/post IBS-SSS scores
gPlot_severities <- ggplot(df_metadata, aes(x=pre_diet_IBS_SSS, y=post_diet_IBS_SSS, 
                                              shape=reference, fill=Response ))+
  geom_point(size=2)+
  scale_x_continuous("IBS Symptom Severity\n(Before Low-FODMAP Diet)", limits = c(0,500), breaks=c(0, 250, 500))+
  scale_y_continuous("IBS Symptom Severity\n(After Low-FODMAP Diet)", limits = c(0,500), breaks=c(0, 250, 500))+
  scale_shape_manual("Study", values = c(21, 22, 23, 24, 25))+
  geom_abline(intercept = -No_Response_Max_Improvement, slope = 1, linetype = "dashed")+
  geom_abline(intercept = -Low_Response_Max_Improvement, slope = 1, linetype = "dashed")+
  scale_fill_manual(values=c("No"=unname(colors_assigned["D"]),
                             "Low"=unname(colors_assigned["B"]),
                             "High"=unname(colors_assigned["A"])))+
  ggtitle(sprintf("(No. of patients = %d)", nrow(df_metadata)))+
  guides(shape = guide_legend(order = 1),
         fill = guide_legend(override.aes=list(shape=21), order = 2))+
  coord_equal()+
  my_base_theme %+replace% theme(legend.position="bottom",
                                  legend.direction = "vertical",
                                  legend.justification = "left")
print(gPlot_severities)

###############################################################################
############## 2. Differential Abundant Analysis of KEGG Pathways #############
###############################################################################
# 2.A) Load pathway abundances (16S rRNA only)
df_metadata <- df_metadata[df_metadata$technology == "16S rRNA" &
                           !is.na(df_metadata$Response) &
                           !is.na(df_metadata$pre_diet_microbiome_sample_id),]
rownames(df_metadata) <- df_metadata$pre_diet_microbiome_sample_id
pathway_abundances <- load_pathway_abundances(INTEGRATED_IBS_DIR, df_metadata)

# 2.B) Identify differentially abuntant pathways
data_High_No <- get_selected_labels(pathway_abundances$df_otu, df_metadata, selected_labels = c("High", "No"))
df_diff_abundances_High_No <- get_diff_abundances(df_otu = data_High_No$df_otu, labels = data_High_No$labels)

# 2.C) Box-plot top differentially abundant pathways
gPlot_Pathways <- plot_diff_abundances(df_otu = data_High_No$df_otu, 
                                       labels = data_High_No$labels, 
                                       diff_abundance_info = df_diff_abundances_High_No,
                                       friendly_names = pathway_abundances$friendly_names,
                                       top_n = 5, 
                                       feature_group_name = "KEGG Pathways")
print(gPlot_Pathways)

###############################################################################
############### 3. Differential Abundant Analysis of Genus Taxa ###############
###############################################################################
# 3.A) Load taxa abundances
taxa_abundances <- load_taxa_abundances(INTEGRATED_IBS_DIR, df_metadata)

# 3.B) Identify differentially abuntant genus taxa
data_High_No <- get_selected_labels(taxa_abundances$df_otu, df_metadata, selected_labels = c("High", "No"))
df_diff_abundances_High_No <- get_diff_abundances(df_otu = data_High_No$df_otu, labels = data_High_No$labels)

# 3.C) Box-plot top differentially abundant taxa
gPlot_Taxa <- plot_diff_abundances(df_otu = data_High_No$df_otu,
                                   labels = data_High_No$labels, 
                                   diff_abundance_info = df_diff_abundances_High_No[!is.na(df_diff_abundances_High_No$p.value),],
                                   friendly_names = taxa_abundances$friendly_names,
                                   top_n = 5, 
                                   feature_group_name = "Genus Taxa")
print(gPlot_Taxa)

###############################################################################
############### 4. Combine plots ##############################################
###############################################################################
gPlot_B_C <- cowplot::plot_grid(gPlot_Pathways, NULL, gPlot_Taxa, ncol = 1,
                                rel_heights = c(1, 0.1, 1), rel_widths = c(1, 1, 1), labels = c("B", "","C"),
                                align = "h", axis = "lr")

gPLot_comb <- cowplot::plot_grid(gPlot_severities, NULL, gPlot_B_C, ncol = 3,
                                 rel_heights = c(1,1,1), rel_widths = c(2, 0.2, 5), labels = c("A"),
                                 align = "h", axis = "c")

lib.util.ggsave(file.path("results", "Fig2.png"),
                gPLot_comb, width = 10, height = 5, dpi=400)
lib.util.ggsave(file.path("results", "Fig2.pdf"),
                gPLot_comb, width = 10, height = 5, dpi=400)

