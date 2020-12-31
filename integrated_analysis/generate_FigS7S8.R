#!/usr/bin/env Rscript
# Linear discriminant analysis (LDA) (Figure S7 for SOM)

# Libraries
library(phyloseq)
library(biomformat)
library(ggplot2)
library(ggrepel)
library(microbiome)
library(stringr)
library(reshape2)
library(MASS)

# Shared code
source("lib/util.R")
source("lib/config_management.R")
source("lib/metadata_management.R")
source("lib/abundance_management.R")
source("lib/theme_util.R")


# Load metadata
df_metadata <- load_metadata()
df_metadata <- df_metadata[df_metadata$technology == "16S rRNA" &
                             !is.na(df_metadata$Response) &
                             !is.na(df_metadata$pre_diet_microbiome_sample_id),]
rownames(df_metadata) <- df_metadata$pre_diet_microbiome_sample_id

###############################################################################
####################### A) LDA for Pathway Abundances #########################
###############################################################################

# A.1) Load pathway abundances
pathway_abundances <- load_pathway_abundances(INTEGRATED_IBS_DIR, df_metadata)
data_High_No <- get_selected_labels(pathway_abundances$df_otu, df_metadata, selected_labels = c("High", "No"))
df_diff_abundances_High_No <- get_diff_abundances(df_otu = data_High_No$df_otu, labels = data_High_No$labels)

# A.2 Select otus based on p-value threshold alpha
top_n <- 10
data.table::setorderv(df_diff_abundances_High_No, cols = c("q.value", "p.value", "name"), order = 1)
df_selected_otus <- df_diff_abundances_High_No[1:top_n, c("name", "p.value")]
df_otu <- data_High_No$df_otu[,df_selected_otus$name]

# A.3) LDA
lda_res = MASS::lda(x=df_otu, grouping = data_High_No$labels)
df_otu_weights <- data.frame(weight=lda_res$scaling[,1],
                             pvalue=df_selected_otus[df_selected_otus$name == names(lda_res$scaling[,1]), c("p.value")],
                             name=pathway_abundances$friendly_names[names(lda_res$scaling[,1])])
df_otu_weights$name = str_wrap(df_otu_weights$name, width = 45)

# A.4) Plot weights
gPlot_Weights <- ggplot(df_otu_weights, aes(x=reorder(name, -pvalue) , y=abs(weight), fill=factor(sign(weight), levels = c(1, -1))))+
                        geom_bar(stat="identity") + 
                        scale_y_continuous("| LDA Coefficient |")+
                        scale_fill_manual("", values = c("1" = unname(colors_assigned["B"]),
                                                        "-1" = unname(colors_assigned["D"])),
                                              labels = c("1"="Positive", "-1" = "Negative")) +
                        xlab("Feature name") +
                        coord_flip() +
                        my_base_theme %+replace% theme(legend.position = "top")

gPlot_Pvalues <- ggplot(df_otu_weights, aes(x=reorder(name, -pvalue) , y=pvalue))+
                        geom_bar(stat="identity", fill=unname(colors_assigned["A"])) +
                        scale_y_continuous("p-value")+
                        coord_flip() +
                        my_base_theme %+replace% theme(axis.title.y=element_blank(),
                                                       axis.text.y=element_blank(),
                                                       axis.ticks.y=element_blank())
gPlot_comb <- cowplot::plot_grid(gPlot_Weights, gPlot_Pvalues, ncol = 2,
                                 rel_heights = c(1,1), rel_widths = c(4, 1),
                                 align = "h", axis = "tb")

lib.util.ggsave(file.path("results", "FigS7.png"),
                gPlot_comb, width = 10, height = 5, dpi=400)
lib.util.ggsave(file.path("results", "FigS7.pdf"),
                gPlot_comb, width = 10, height = 5, dpi=400)



###############################################################################
######################## B) LDA for Taxa Abundances ###########################
###############################################################################
# B.1) Load Taxa abundances
taxa_abundances <- load_taxa_abundances(INTEGRATED_IBS_DIR, df_metadata)

# B.2) Identify differentially abuntant genus taxa
data_High_No <- get_selected_labels(taxa_abundances$df_otu, df_metadata, selected_labels = c("High", "No"))
df_diff_abundances_High_No <- get_diff_abundances(df_otu = data_High_No$df_otu, labels = data_High_No$labels)
df_diff_abundances_High_No <- df_diff_abundances_High_No[!is.nan(df_diff_abundances_High_No$p.value),]

# A.2 Select otus based on p-value threshold alpha
top_n <- 10
data.table::setorderv(df_diff_abundances_High_No, cols = c("q.value", "p.value", "name"), order = 1)
df_selected_otus <- df_diff_abundances_High_No[1:top_n, c("name", "p.value")]
df_otu <- data_High_No$df_otu[,df_selected_otus$name]

# A.3) LDA
lda_res = MASS::lda(x=df_otu, grouping = data_High_No$labels)
df_otu_weights <- data.frame(weight=lda_res$scaling[,1],
                             pvalue=df_selected_otus[df_selected_otus$name == names(lda_res$scaling[,1]), c("p.value")],
                             name=taxa_abundances$friendly_names[names(lda_res$scaling[,1])])
df_otu_weights$name = str_wrap(df_otu_weights$name, width = 45)

# A.4) Plot weights
gPlot_Weights <- ggplot(df_otu_weights, aes(x=reorder(name, -pvalue) , y=abs(weight), fill=factor(sign(weight), levels = c(1, -1))))+
  geom_bar(stat="identity") + 
  scale_y_continuous("| LDA Coefficient |")+
  scale_fill_manual("", values = c("1" = unname(colors_assigned["B"]),
                                   "-1" = unname(colors_assigned["D"])),
                    labels = c("1"="Positive", "-1" = "Negative")) +
  xlab("Feature name") +
  coord_flip() +
  my_base_theme %+replace% theme(legend.position = "top")

gPlot_Pvalues <- ggplot(df_otu_weights, aes(x=reorder(name, -pvalue) , y=pvalue))+
  geom_bar(stat="identity", fill=unname(colors_assigned["A"])) +
  scale_y_continuous("p-value", breaks = seq(0, 0.1, 0.02))+
  coord_flip() +
  my_base_theme %+replace% theme(axis.title.y=element_blank(),
                                 axis.text.y=element_blank(),
                                 axis.ticks.y=element_blank())
gPlot_comb <- cowplot::plot_grid(gPlot_Weights, gPlot_Pvalues, ncol = 2,
                                 rel_heights = c(1,1), rel_widths = c(4, 1),
                                 align = "h", axis = "tb")
print(gPlot_comb)     
lib.util.ggsave(file.path("results", "FigS8.png"),
                gPlot_comb, width = 10, height = 5, dpi=400)
lib.util.ggsave(file.path("results", "FigS8.pdf"),
                gPlot_comb, width = 10, height = 5, dpi=400)
