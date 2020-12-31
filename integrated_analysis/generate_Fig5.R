#!/usr/bin/env Rscript
#' Generate Fig 5: Prediction of diet (low-FODMAP vs. other) given microbiome data
#' 1) Predict using pathway abundances (16S rRNA data)
#' 2) Predict using genus abundances (16S rRNA data)

# Libraries
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
########## 1. Predict diet using pathway abundances (16S rRNA data) ###########
###############################################################################
# 1.A) Load data
df_metadata <- load_metadata_diet_all_16S()
pathway_abundances <- load_pathway_abundances(INTEGRATED_ALL_16S_DIR, df_metadata)
df_metadata <- df_metadata[rownames(pathway_abundances$df_otu),]

# 1.B) Plot classification performance for  prediction (True vs. False low-FODMAP diet)
df_comb <- cbind(pathway_abundances$df_otu[rownames(df_metadata),], df_metadata)
gPlot_ML_Pathways <- plot_classification_perf(df_comb,
                                              microbiome_cols = colnames(pathway_abundances$df_otu),
                                              binary_label_groups = list(true=c("TRUE"), false=c("FALSE")),
                                              friendly_names = pathway_abundances$friendly_names,
                                              entity_label = "pathways",
                                              label_col = "low_fodmap")
print(gPlot_ML_Pathways$main)

###############################################################################
########### 2. Predict diet using genus abundances (16S rRNA data) ############
###############################################################################
# 2.A) Load data
df_metadata <- load_metadata_diet_all_16S()
taxa_abundances <- load_taxa_abundances(INTEGRATED_ALL_16S_DIR, df_metadata)
df_metadata <- df_metadata[rownames(taxa_abundances$df_otu),]

# 2.B) Plot classification performance for  prediction (True vs. False low-FODMAP diet)
df_comb <- cbind(taxa_abundances$df_otu[rownames(df_metadata),], df_metadata)
gPlot_ML_Taxa <- plot_classification_perf(df_comb,
                                              microbiome_cols = colnames(taxa_abundances$df_otu),
                                              binary_label_groups = list(true=c("TRUE"), false=c("FALSE")),
                                              friendly_names = taxa_abundances$friendly_names,
                                              entity_label = "taxa",
                                              label_col = "low_fodmap")
print(gPlot_ML_Taxa$main)

###############################################################################
########################## 3. Combine Plots and Save  #########################
###############################################################################
align_top <- cowplot::align_plots(gPlot_ML_Pathways$top_left, gPlot_ML_Pathways$top_right,
                                  gPlot_ML_Taxa$top_left, gPlot_ML_Taxa$top_right,
                                  align = "h", axis = "tb" )

align_bottom <- cowplot::align_plots(gPlot_ML_Pathways$bottom,
                                     gPlot_ML_Taxa$bottom,
                                     align = "h", axis = "tb")

gPlotB_top <- cowplot::plot_grid(align_top[[1]], align_top[[2]], nrow=1, ncol=2, align = "h")
gPlotB <- cowplot::plot_grid(gPlotB_top, align_bottom[[1]], nrow=2, ncol=1,
                             rel_heights = c(1,2.2), labels = c("", "D"))

gPlotC_top <- cowplot::plot_grid(align_top[[3]], align_top[[4]], nrow=1, ncol=2, align = "h")
gPlotC <- cowplot::plot_grid(gPlotC_top, align_bottom[[2]], nrow=2, ncol=1,
                             rel_heights = c(1,2.2), labels = c("", "E"))

gPlot_B_C <- cowplot::plot_grid(gPlotB, gPlotC, nrow = 1, ncol = 2, labels = c("B", "C"))

gPlotA <- cowplot::ggdraw()+
  cowplot::draw_image(file.path("./fig/Fig5A.png"))

gPlot_all <- cowplot::plot_grid(gPlotA, gPlot_B_C, nrow = 2, ncol = 1, labels = c("A", ""),
                                rel_heights = c(2, 7))

lib.util.ggsave(file.path("results", "Fig5.png"),
              gPlot_all, width = 6, height = 7, dpi=400)

lib.util.load_fonts()
lib.util.ggsave(file.path("results", "Fig5.pdf"),
              gPlot_all, width = 6, height = 7, dpi=400)
