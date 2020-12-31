#!/usr/bin/env Rscript
#' Generate Fig 4: Prediction of response to low-FODMAP diet given pre-diet microbial
#'                 abundances for methane and fatty acid metabolism pathways
#'
#' 1) Differential abundances
#' 2) Predict High vs. No response
#' 3) Predict High vs. Low response

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
my_base_theme <- get_base_theme(0.8)

###############################################################################
######################### 1. Differential Abundances ##########################
###############################################################################
# 1.A) Load data
df_metadata <- load_metadata()
df_metadata <- df_metadata[df_metadata$technology == "16S rRNA" &
                             !is.na(df_metadata$Response) &
                             !is.na(df_metadata$pre_diet_microbiome_sample_id),]
rownames(df_metadata) <- df_metadata$pre_diet_microbiome_sample_id
pathway_abundances <- load_pathway_abundances(INTEGRATED_IBS_DIR, df_metadata)
df_comb <- cbind(pathway_abundances$df_otu[rownames(df_metadata),], df_metadata)

# 1.B) Create labels
response_freq <- plyr::count(df_comb$Response)
get_label <- function(response){
  return(sprintf("%s\n(n=%d)", response, response_freq[response_freq$x == response, c("freq")]))
}
response_labels <- get_label(levels(df_metadata$Response))
names(response_labels) <- levels(df_metadata$Response)

# 1.C) Prepare for p-value comparison
my_comparisons <- list( c("High", "No"),
                        c("High", "Low"),
                        c("Low", "No"))

# 1.D) Methane metabolism
gPlot_methane <- ggplot(df_comb, aes(x=Response, 
                                     y=Methane.metabolism, 
                                     fill=Response))+
  stat_compare_means(comparisons = my_comparisons, size=2,
                     label.y = c(2.7, 2.8, 2.88),
                     method.args = list(alternative = "greater"))+
  scale_x_discrete("Low-FODMAP diet response", labels = response_labels)+
  scale_y_continuous("CLR-Transformed abundance of\nmicrobiome methane metabolism")+
  expand_limits(y=2.9)+
  geom_boxplot(outlier.shape = 21)+
  scale_fill_manual(values=c("No"=unname(colors_assigned["D"]),
                             "Low"=unname(colors_assigned["B"]),
                             "High"=unname(colors_assigned["A"])))+
  my_base_theme %+replace% theme(legend.position = "none")
print(gPlot_methane)

# 1.E) Fatty Acid metabolism
gPlot_fatty_acid <- ggplot(df_comb, aes(x=Response, 
                                        y=Fatty.acid.metabolism, 
                                        fill=Response))+
  stat_compare_means(comparisons = my_comparisons, size=2,
                     label.y = c(1, 1.12, 1.2),
                     method.args = list(alternative = "greater"))+
  scale_x_discrete("Low-FODMAP diet response", labels = response_labels)+
  scale_y_continuous("CLR-Transformed abundance of\nmicrobiome fatty acid metabolism")+
  geom_boxplot(outlier.shape = 21)+
  scale_fill_manual(values=c("No"=unname(colors_assigned["D"]),
                             "Low"=unname(colors_assigned["B"]),
                             "High"=unname(colors_assigned["A"])))+
  my_base_theme %+replace% theme(legend.position = "none")
print(gPlot_fatty_acid)

###############################################################################
######################## 2. Predict High vs. No response ######################
###############################################################################
selected_pathways <- c("Methane.metabolism", "Fatty.acid.metabolism")
df_res <- s_classify_binary_RF_loo_cr(df_comb, x_cols = selected_pathways, label_col = "Response",
                                      label_group = list(true=c("High"), false=c("No")))
gPlot_High_vs_No_roc <- get_roc_plot(df_res, color_code = "D") + 
  ggtitle("Predict High vs. No\n response given methane \nand fatty acid abundances")
print(gPlot_High_vs_No_roc)

gPlot_High_vs_No_pr <- get_pr_plot(df_res, color_code = "D"); print(gPlot_High_vs_No_pr)
message(sprintf("F1-Score: %.3f", get_F1_score(df_res)))

###############################################################################
###################### 3. Predict High vs. Low/No response ####################
###############################################################################
df_res <- s_classify_binary_RF_loo_cr(df_comb, x_cols = selected_pathways, label_col = "Response",
                                      label_group = list(true=c("High"), false=c("Low", "No")))
gPlot_High_vs_LowNo_roc <- get_roc_plot(df_res, color_code = "E")+
  ggtitle("Predict High vs. Low or No\n response given methane \n and fatty acid abundances")
print(gPlot_High_vs_LowNo_roc)
gPlot_High_vs_LowNo_pr <- get_pr_plot(df_res, color_code = "E"); print(gPlot_High_vs_LowNo_pr)

###############################################################################
######################### 4. Combine Plots and Save ###########################
###############################################################################
gPlot <- cowplot::plot_grid(gPlot_methane, gPlot_High_vs_No_roc , gPlot_High_vs_LowNo_roc,
                            gPlot_fatty_acid, gPlot_High_vs_No_pr, gPlot_High_vs_LowNo_pr,
                            ncol = 3, nrow = 2, labels = c("A", "C", "E", "B", "D", "F"),
                            align = "hv", axis = "tblr")
print(gPlot)
lib.util.ggsave(file.path("results", "Fig4.png"),
              gPlot, dpi = 400, width = 5.5, height = 4)
lib.util.ggsave(file.path("results", "Fig4.pdf"),
              gPlot, dpi = 400, width = 5.5, height = 4)

