#!/usr/bin/env Rscript
#' Generate Fig S4: Cumulative percentage of ΔIBS-SSS for 152 IBS patients 
#'                  following a low-FODMAP diet.

library(ggplot2)
library(ggrepel)
library(stringr)
library(reshape2)

# Shared code
source("lib/util.R")
source("lib/config_management.R")
source("lib/metadata_management.R")
source("lib/theme_util.R")

# 1.A) Load data
df_meta <- load_metadata()

# 1.B) Identify cumulative percentage thresholds
t1 <- 124 # Reported mean of ΔIBS-SSS - STD for placebo treatment
t2 <- 162 # Reported mean of ΔIBS-SSS + STD for placebo treatment

df_meta$impr_IBS_SSS <- df_meta$pre_diet_IBS_SSS - df_meta$post_diet_IBS_SSS
df_meta <- df_meta[!is.na(df_meta$impr_IBS_SSS),]
y_t1 <- sum(df_meta$impr_IBS_SSS < t1)/nrow(df_meta)
y_t2 <- sum(df_meta$impr_IBS_SSS < t2)/nrow(df_meta)
df_lines = data.frame(x=c(t1, t2, -Inf, -Inf), 
                      y=c(-Inf, -Inf, y_t1, y_t2), 
                      xend=c(t1, t2, t1, t2), 
                      yend=c(y_t1, y_t2, y_t1, y_t2))

# 1.C) Plot and save
gPlot_CP <- ggplot(df_meta, aes(x=impr_IBS_SSS))+
  stat_ecdf(geom="area", fill=colors_assigned["A"])+
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend), data = df_lines,
               linetype="dashed")+
  scale_y_continuous("Cumulative percentage", labels=scales::percent, breaks = c(0, 0.2, 0.4, y_t1, y_t2, 0.80, 1.0))+
  scale_x_continuous(expression(bold(Delta[IBS-SSS])), breaks = c(-100, 0, 100, t1, t2, 200, 300))+
  my_base_theme

print(gPlot_CP)
lib.util.ggsave(file.path("results", "FigS4.png"), gPlot_CP, dpi = 400)