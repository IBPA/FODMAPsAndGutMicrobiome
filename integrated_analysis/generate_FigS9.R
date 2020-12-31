#!/usr/bin/env Rscript
#' Generate Fig S9: Improvement in IBS-SSS score amongst studies with 16S rRNA 
#'                  microbiome data. 

# Libraries
library(ggplot2)
library(ggrepel)
library(stringr)
library(reshape2)
library(ggpubr)

# Shared code
source("lib/util.R")
source("lib/config_management.R")
source("lib/metadata_management.R")
source("lib/theme_util.R")

df_metadata <- load_metadata()
df_metadata <- df_metadata[df_metadata$technology == "16S rRNA" &
                             !is.na(df_metadata$Response) &
                             !is.na(df_metadata$pre_diet_microbiome_sample_id),]
# Prepare for p-value comparison
my_comparisons <- list(c("(Harvie et al., 2017)", "(Schumann et al., 2018)"),
                       c("(Harvie et al., 2017)", "(McIntosh et al., 2017)"),
                       c("(McIntosh et al., 2017)", "(Schumann et al., 2018)"))

gPlot <- ggplot(df_metadata, aes(x=reference, y=pre_diet_IBS_SSS-post_diet_IBS_SSS))+
                stat_compare_means(comparisons = my_comparisons)+
                geom_boxplot(outlier.shape = 21, fill=colors_assigned["A"])+
                xlab("Study")+
                scale_y_continuous("Improvement in IBS-SSS after following a low-FODMAP diet")+
                my_base_theme %+replace% theme(legend.position = "none")
print(gPlot)

lib.util.ggsave(file.path("results", "FigS9.png"),
                gPlot, dpi=400)
lib.util.ggsave(file.path("results", "FigS9.pdf"),
                gPlot, dpi=400)
