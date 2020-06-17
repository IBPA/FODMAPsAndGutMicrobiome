#' Generate Fig2
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
source("lib/config_management.R")
source("lib/metadata_management.R")
source("lib/theme_util.R")

###############################################################################
############## 1. Plot distribution of pre/post IBS-SSS scores ################
###############################################################################
# 1.A) Load metadata from all IBS patients
df_metadata <- load_metadata()

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
# 2.A) Load pathway abundances

# 2.B) Identify differentially abuntant pathways

# 2.C) Box-plot top differentially abundant pathways

###############################################################################
############### 3. Differential Abundant Analysis of Genus Taxa ###############
###############################################################################
# 3.A) Load taxa abundances

# 3.B) Identify differentially abuntant genus taxa

# 3.C) Box-plot top differentially abundant taxa

###############################################################################
############### 4. Combine plots ##############################################
###############################################################################
