# Beta diversity across three studies
# ******** Not doing this *********

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

df_metadata <- load_metadata()
df_metadata <- df_metadata[df_metadata$technology == "16S rRNA" &
                             !is.na(df_metadata$Response) &
                             !is.na(df_metadata$pre_diet_microbiome_sample_id),]
rownames(df_metadata) <- df_metadata$pre_diet_microbiome_sample_id
pathway_abundances <- load_pathway_abundances(INTEGRATED_IBS_DIR, df_metadata)

pca_res <- prcomp(pathway_abundances$df_otu)
df_plot <- data.frame(pca_res$x[,c("PC1", "PC2")])
df_plot <- cbind(df_plot, df_metadata$reference)
colnames(df_plot) <- c("PC1", "PC2", "Study")

ggplot(df_plot, aes(x=PC1, y=PC2, colour=Study))+
  geom_point()

#plot_landscape(taxa_abundances$ps, "t-SNE", distance = "euclidean", perplexity=2)
#plot_ordination(taxa_abundances$ps, ps.ord, type="samples")
