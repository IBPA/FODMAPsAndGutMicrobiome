#!/usr/bin/env Rscript
#' Generate Fig 3: Prediction of response to low-FODMAP diet given pre-diet 
#'                 microbiome data
#' 1) Predict using pathway abundances (16S rRNA data)
#' 2) Predict using genus abundances (16S rRNA data)
#' 3) Predict using GA-map probes
#' 4) Visualize performance using ROC,PR and recursive feature elimination plots

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

###############################################################################
######## 1. Predict response using pathway abundances (16S rRNA data) #########
###############################################################################
# 1.A) Load data
df_metadata <- load_metadata()
df_metadata <- df_metadata[df_metadata$technology == "16S rRNA" &
                             !is.na(df_metadata$Response) &
                             !is.na(df_metadata$pre_diet_microbiome_sample_id),]
rownames(df_metadata) <- df_metadata$pre_diet_microbiome_sample_id
pathway_abundances <- load_pathway_abundances(INTEGRATED_IBS_DIR, df_metadata)

# 1.B) Plot classification performance for  prediction High vs. Low/No response
df_comb <- cbind(pathway_abundances$df_otu[rownames(df_metadata),], df_metadata)
gPlot_ML_Pathways <- plot_classification_perf(df_comb, 
                                              microbiome_cols = colnames(pathway_abundances$df_otu),
                                              binary_label_groups = list(true=c("High"), false=c("Low", "No")),
                                              friendly_names = pathway_abundances$friendly_names,
                                              entity_label = "pathways",
                                              label_col = "Response")
print(gPlot_ML_Pathways$main)

###############################################################################
######### 2. Predict response using genus abundances (16S rRNA data) ##########
###############################################################################
# 2.A) Load data
df_metadata <- load_metadata()
df_metadata <- df_metadata[df_metadata$technology == "16S rRNA" &
                             !is.na(df_metadata$Response) &
                             !is.na(df_metadata$pre_diet_microbiome_sample_id),]
rownames(df_metadata) <- df_metadata$pre_diet_microbiome_sample_id
taxa_abundances <- load_taxa_abundances(INTEGRATED_IBS_DIR, df_metadata)

# 2.B) Plot classification performance for  prediction High vs. Low/No response
df_comb <- cbind(taxa_abundances$df_otu[rownames(df_metadata),], df_metadata)
gPlot_ML_Taxa <- plot_classification_perf(df_comb, 
                                          microbiome_cols = colnames(taxa_abundances$df_otu),
                                          binary_label_groups = list(true=c("High"), false=c("Low", "No")),
                                          friendly_names = taxa_abundances$friendly_names,
                                          entity_label = "pathways",
                                          label_col = "Response")
print(gPlot_ML_Taxa$main)

###############################################################################
################## 3. Predict response using GA-map probes ####################
###############################################################################
# 3.A) Load data
df_metadata <- load_metadata()
df_metadata <- df_metadata[df_metadata$technology == "GA-map" &
                          !is.na(df_metadata$Response),]
rownames(df_metadata) <- df_metadata$host_id
GA_map_abundances <- load_GA_map_abundances(INTEGRATED_GA_MAP_DIR, df_metadata)

# 3.B) Plot classification performance for  prediction High vs. Low/No response
df_comb <- cbind(GA_map_abundances$df_otu[rownames(df_metadata),], df_metadata)
gPlot_ML_GA_map <- plot_classification_perf(df_comb, 
                                            microbiome_cols = colnames(GA_map_abundances$df_otu),
                                            binary_label_groups = list(true=c("High"), false=c("Low", "No")),
                                            friendly_names = GA_map_abundances$friendly_names,
                                            entity_label = "probes",
                                            label_col = "Response")

print(gPlot_ML_GA_map$main)


###############################################################################
########################## 4. Combine Plots and Save  #########################
###############################################################################
align_top <- cowplot::align_plots(gPlot_ML_Pathways$top_left, gPlot_ML_Pathways$top_right,
                                  gPlot_ML_Taxa$top_left, gPlot_ML_Taxa$top_right,
                                  gPlot_ML_GA_map$top_left, gPlot_ML_GA_map$top_right,
                                  align = "h", axis = "tb" )

align_bottom <- cowplot::align_plots(gPlot_ML_Pathways$bottom,
                                     gPlot_ML_Taxa$bottom,
                                     gPlot_ML_GA_map$bottom,
                                     align = "h", axis = "tb")


gPlotA_top <- cowplot::plot_grid(align_top[[1]], align_top[[2]], nrow=1, ncol=2, align = "h")
gPlotA <- cowplot::plot_grid(gPlotA_top, align_bottom[[1]], nrow=2, ncol=1,
                             rel_heights = c(1,2.2), labels = c("", "D"))

gPlotB_top <- cowplot::plot_grid(align_top[[3]], align_top[[4]], nrow=1, ncol=2, align = "h")
gPlotB <- cowplot::plot_grid(gPlotB_top, align_bottom[[2]], nrow=2, ncol=1,
                             rel_heights = c(1,2.2), labels = c("", "E"))

gPlotC_top <- cowplot::plot_grid(align_top[[5]], align_top[[6]], nrow=1, ncol=2, align = "h")
gPlotC <- cowplot::plot_grid(gPlotC_top, align_bottom[[3]], nrow=2, ncol=1,
                             rel_heights = c(1,2.2), labels = c("", "F"))

gPlot <- cowplot::plot_grid(gPlotA, gPlotB, gPlotC, nrow = 1, ncol = 3, labels = "AUTO")
lib.util.ggsave(file.path("results", "Fig3.png"),
              gPlot, width = 10, height = 6, dpi=400)

lib.util.load_fonts()
lib.util.ggsave(file.path("results", "Fig3.pdf"),
              gPlot, width = 10, height = 6, dpi=400)


