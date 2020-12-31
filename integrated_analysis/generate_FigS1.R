#!/usr/bin/env Rscript
#' Generate Fig S1: GA-map analysis: differential abundance of methane biomarker,
#'                  and low-FODMAP response prediction in (Bennet et al., 2018)
#' 1) Differential abundance of methane metabolism biomarker (3 taxa)
#' 2) Differential abundance of methane metabolism biomarker (1 taxa)
#' 3) Predict response using all 3 methane metabolism biomarkers

# Library
library(ggplot2)
library(ggrepel)
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

###############################################################################
##### 1. Differential abundance of methane metabolism biomarker (3 taxa) ######
###############################################################################

# 1.A) Load data
df_metadata <- load_metadata()
df_metadata <- df_metadata[df_metadata$reference == "(Bennet et al., 2018)" &
                             !is.na(df_metadata$Response),]
rownames(df_metadata) <- df_metadata$host_id
GA_map_abundances <- load_GA_map_abundances(INTEGRATED_GA_MAP_DIR, df_metadata)
df_comb <- cbind(GA_map_abundances$df_otu[rownames(df_metadata),], df_metadata)

# 1.B) Create Methane metabolism biomarker using GA-map probes
# We use probes that detect fecal taxa with reported negative correlations to  
# breath methane according to (Parthasarathy et al., 2016), Table 3:
#  1) AG0581: Dorea [g] (specified in the manuscript)
#  2) AG1267: Pseudomonas [g] (missed to specify in the manuscript)
#  3) IG0005: Proteobacteria [ph] (missed to specify in the manuscript)
df_comb$Methane.metabolism <-  -(df_comb$AG0581 + df_comb$AG1267 + df_comb$IG0005)

# 1.C) Create labels
response_freq <- plyr::count(df_comb$Response)
get_label <- function(response){
  return(sprintf("%s\n(n=%d)", response, response_freq[response_freq$x == response, c("freq")]))
}
response_labels <- get_label(levels(df_metadata$Response))
names(response_labels) <- levels(df_metadata$Response)

# 1.D) Prepare for p-value comparison
my_comparisons <- list( c("High", "No"),
                        c("High", "Low"),
                        c("Low", "No"))

# 1.E) Plot
gPlot_DA_3 <- ggplot(df_comb, aes(x=Response, y=Methane.metabolism, fill=Response))+
                    geom_boxplot()+
                    scale_fill_manual(values=c("No"=unname(colors_assigned["D"]),
                                               "Low"=unname(colors_assigned["B"]),
                                               "High"=unname(colors_assigned["A"])))+
                    stat_compare_means(comparisons = my_comparisons, size=2,
                                       label.y = c(3.1, 3.8, 4.2),
                                       method.args = list(alternative = "greater"))+
                    scale_x_discrete("Low-FODMAP diet response", labels = response_labels)+
                    scale_y_continuous("Microbiome methane metabolism")+
                    geom_boxplot(outlier.shape = 21)+
                    expand_limits(y=4.5)+
                    my_base_theme %+replace% theme(legend.position = "none")
print(gPlot_DA_3)

###############################################################################
##### 2. Differential abundance of methane metabolism biomarker (1 taxa) ######
###############################################################################
# 2) Use Dorea genus only as biomarker and plot 
df_comb$Methane.metabolism <-  -(df_comb$AG0581)
gPlot_DA_1 <- ggplot(df_comb, aes(x=Response, y=Methane.metabolism, fill=Response))+
  geom_boxplot()+
  scale_fill_manual(values=c("No"=unname(colors_assigned["D"]),
                             "Low"=unname(colors_assigned["B"]),
                             "High"=unname(colors_assigned["A"])))+
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3.1, 3.8, 4.2), size=2,
                     method.args = list(alternative = "greater"))+
  scale_x_discrete("Low-FODMAP diet response", labels = response_labels)+
  scale_y_continuous("Microbiome methane metabolism \n (Dorea genus only)")+
  geom_boxplot(outlier.shape = 21)+
  expand_limits(y=4.5)+
  my_base_theme %+replace% theme(legend.position = "none")
print(gPlot_DA_1)

###############################################################################
######## 3. Predict response using all 3 methane metabolism biomarkers ########
###############################################################################

# 3.A) Classify LOO CV (High vs. No)
methane_probes <- c("AG0581", "AG1267", "IG0005")
df_res <- s_classify_binary_RF_loo_cr(df_comb, x_cols = methane_probes, label_col = "Response",
                                        label_group = list(true=c("High"), false=c("No")))
gPlot_roc <- get_roc_plot(df_res, color_code = "D", add_star = FALSE)
gPlot_pr <- get_pr_plot(df_res, color_code = "D", add_star = FALSE)

###############################################################################
######################## 4. Combine and save figure ###########################
###############################################################################
gPlot <- cowplot::plot_grid(gPlot_DA_3, gPlot_DA_1, gPlot_roc, gPlot_pr,
                            nrow=2, ncol=2, labels = c("A", "B", "C", "D"),
                            align = "hv", axis = "tblr")

lib.util.ggsave(file.path("results", "FigS1.png"),
                gPlot, dpi = 400, width = 5, height = 5)
