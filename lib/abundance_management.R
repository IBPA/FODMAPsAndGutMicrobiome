#' Utility functions for management of microbiome pathway abundances

multi_study_filter <- function(ps, num_min_shared = 1){
  df_metadata <- meta(ps)
  if( !"study_name" %in% colnames(df_metadata) || 
      length(unique(df_metadata$study_name)) <= 1 ){
    return(ps) # Filter is not applicable, return
  }
  
  lib.util.log("multiple studies, filter out non-shared OTUs")
  otu_intersect <- taxa(ps)
  for(study_name_curr in unique(df_metadata$study_name)){
    sample_names_curr <- df_metadata[df_metadata$study_name == study_name_curr, c("name")]
    ps_curr <- prune_samples(df_metadata$study_name == study_name_curr, ps)
    otu_keep <- filter_taxa(ps_curr, function(x) sum(x) >= num_min_shared, FALSE)
    otu_curr <- names(otu_keep[otu_keep == TRUE])
    otu_intersect <- intersect(otu_intersect, otu_curr)
  }
  
  return(prune_taxa(otu_intersect, ps))
}

apply_clr <- function(ps.prop){
  otu_t <- otu_table(ps.prop)
  df_otu <- otu_t@.Data
  
  # Ensure samples are in rows
  if (otu_t@taxa_are_rows == TRUE) {
    df_otu <- t(df_otu)
  }
  
  # Trasnform
  df_otu_clr <- compositions::clr(df_otu)
  
  # replace zeros with the min transformed - c_delta
  tmp_min <- min(df_otu_clr)
  flat_otu <- as.vector(df_otu_clr)
  c_delta <- sd(flat_otu[flat_otu != 0])/10
  #c_delta <- 0.25
  df_otu_clr_t <- ifelse(df_otu_clr==0, tmp_min - c_delta, df_otu_clr)
  
  if (ps.prop@otu_table@taxa_are_rows) {
    otu_new <- data.frame(t(df_otu_clr_t))
  } else{
    otu_new <- data.frame(df_otu_clr_t)
  }
  
  ps.new <- phyloseq(otu_table(data.frame(otu_new), taxa_are_rows = ps.prop@otu_table@taxa_are_rows),
                     tax_table(ps.prop),
                     sample_data(ps.prop))
  return(ps.new)
}


fix_names <- function(x){
  colnames(x) <- str_replace(colnames(x), "-", ".")
  return(x)
}

get_friendly_names <- function(ps, unfriendly_names){
  friendly_names <- taxa_names(ps)
  names(friendly_names) <- unfriendly_names
  return(friendly_names)
}


get_same <- function(input){
  res <- input
  names(res) <- input
  
  return(res)
}

integrate_data_taxa <- function(src_dirs){
  ps_studies <- list()
  for(src_dir in src_dirs){
    print(src_dir)
    ps_curr <- readRDS(file.path(src_dir, "V4_ps.rds"))
    ps_curr <- phyloseq::filter_taxa(ps_curr, function(x) sum(x > 3) > 1, TRUE)
    ps_studies[src_dir]<- ps_curr
  }
  
  merged_otus <- do.call(merge_phyloseq, lapply(ps_studies, otu_table))
  merged_otus@.Data <- fix_names(merged_otus@.Data)
  merged_taxa <- do.call(merge_phyloseq, lapply(ps_studies, tax_table))
  merged_sam <- data.frame("name" = colnames(merged_otus@.Data))
  rownames(merged_sam) <- merged_sam$name
  
  ps <- phyloseq(otu_table(merged_otus), 
                 sample_data(merged_sam), 
                 tax_table(merged_taxa))
  
  return(ps)
}

load_GA_map_abundances <- function(src_dir, df_metadata){
  # Load microbiome data
  df <- read.table(file.path(src_dir, "pre_diet.csv"), header = TRUE, sep = ",")
  rownames(df) <- sprintf("%s_%s", df$study_name, df$ID)
  
  # Load probe names
  df_names <- read.table(file.path(src_dir, "probe_names.csv"), header = TRUE, sep = ",")
  friendly_names <- df_names$Bacteria
  names(friendly_names) <- df_names$ProbeID
  
  return(list(df_otu = df[, c(3:ncol(df))], friendly_names = friendly_names))
}

load_taxa_abundances <- function(src_dir, df_metadata, rank = "genus"){
  # Load
  src_filepath <- file.path(src_dir, "DADA2.rds")
  ps <- readRDS(src_filepath)
  
  # Prune
  ps <- prune_samples(rownames(df_metadata), ps)
  df_metadata <- df_metadata[sample_names(ps),]
  sample_data(ps) <- df_metadata
  
  # Transform
  ps <- transform_sample_counts(ps, function(otu) otu/sum(otu))
  ps <- tax_glom(ps, taxrank=rank)
  ps <- apply_clr(ps)
  
  # Data Frame
  df_otu <- data.frame(t(otu_table(ps))@.Data)
  colnames(df_otu) <- ps@tax_table@.Data[colnames(df_otu),c(rank)]
  
  # Friendly names
  friendly_names <- get_same(colnames(df_otu))
  
  return(list(df_otu = df_otu, friendly_names = friendly_names, ps = ps))
}

load_pathway_abundances <- function(src_dir, df_metadata){
  src_filepath <- file.path(src_dir, "PICRUSt_predicted_metagenomes_KEGG_L3.biom")
  predGenome <- read_biom(src_filepath)
  otumaty = fix_names(as(biom_data(predGenome), "matrix"))
  taxmaty = as.matrix(observation_metadata(predGenome), rownames.force=TRUE)
  df_metadata <- df_metadata[rownames(df_metadata) %in% colnames(otumaty),]
  
  ps = phyloseq(otu_table(otumaty, taxa_are_rows=TRUE), 
                tax_table(taxmaty), 
                sample_data(df_metadata))
  
  ps <- multi_study_filter(ps)
  ps0 = filter_taxa(ps, function(x) sum(x) > 2, TRUE)
  
  # * transform to relative abundances within sample *
  ps.prop <- transform_sample_counts(ps0, function(otu) otu/sum(otu))
  
  # CLR Transform
  ps.prop <- apply_clr(ps.prop)
  
  # Change to data frame
  ps.prop <- tax_glom(ps.prop, taxrank="KEGG_Pathways3")
  df_otu <- data.frame(t(otu_table(ps.prop))@.Data)
  
  # Friendly names
  friendly_names <- get_friendly_names(ps.prop, colnames(df_otu))

  return(list(df_otu = df_otu, friendly_names = friendly_names, ps = ps.prop))
}

get_selected_labels <- function(df_otu, df_meta, selected_labels = c("High", "No")){
  selected_rows <- rownames(df_meta[df_meta$Response %in% selected_labels,])
  
  return(list(df_otu = df_otu[selected_rows,],
              labels = df_meta[selected_rows, c("Response")]))
}

get_diff_abundances <- function(df_otu, labels)
{
  otu_names <- colnames(df_otu)
  df_comb <- cbind(df_otu, data.frame(Response=factor(labels)))
  label_levels <- levels(df_comb$Response)
  df_res_diff_abundances <- data.frame(name = character(), p.value=numeric())
  for (otu_name in otu_names){
    v1 <- df_comb[df_comb$Response == label_levels[1], otu_name]
    v2 <- df_comb[df_comb$Response == label_levels[2], otu_name]
    
    df_curr <- data.frame(name = otu_name, p.value=wilcox.test(v1, v2)$p.value)
    df_res_diff_abundances <- rbind(df_res_diff_abundances, df_curr)
  }
  
  df_res_diff_abundances$q.value <- p.adjust(df_res_diff_abundances$p.value, method="fdr")
  
  return(df_res_diff_abundances)
}

plot_diff_abundances <- function(df_otu, labels, diff_abundance_info, friendly_names, top_n = 5, feature_group_name = "KEGG Pathways"){
  data.table::setorderv(diff_abundance_info, cols = c("q.value", "p.value", "name"), order = 1)
  top_otus <- diff_abundance_info[1:top_n, c("name")]
  df_comb <- cbind(df_otu[,top_otus], data.frame(Response=labels))
  df_melt <- melt(df_comb, id.vars = c("Response"))
  colnames(df_melt) <- c("Response", "Feature", "Abundance")
  levels(df_melt$Feature) <- str_wrap(friendly_names[levels(df_melt$Feature)], width = 18)
  
  gPlot <- ggplot(df_melt, aes(x=Feature, y=Abundance, fill=Response))+
    coord_flip()+
    geom_boxplot(outlier.shape = 21, outlier.size = 1, outlier.stroke = .3, lwd=.3)+
    scale_x_discrete("Response", breaks=c("No", "High"))+
    facet_wrap(~Feature, scales = "free", nrow = 1)+
    scale_y_continuous("", breaks = scales::pretty_breaks(n=3))+
    scale_fill_manual(values=c("No"=unname(colors_assigned["D"]),
                               "High"=unname(colors_assigned["A"])))+
    ggtitle(sprintf("CLR-Transformed Relative Abundances In Gut Microbime (Top %d %s)", top_n, feature_group_name))+
    my_base_theme %+replace% theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
                                    legend.position = "none")
  
  df_qvalues <- data.frame(name = factor(top_otus, levels = top_otus),
                           q.value = sprintf("(%.3f)", diff_abundance_info[diff_abundance_info$name %in% top_otus, c("q.value")]))
  label_dic <- as.character(df_qvalues$q.value)
  names(label_dic) <- as.character(df_qvalues$name)
  
  
  gText <- ggplot(df_qvalues, aes(x=name, label=q.value, y="q-value:"))+
    geom_text()+
    scale_x_discrete("", breaks = NULL, expand = c(.1,.1) )+
    scale_y_discrete("")+
    my_base_theme %+replace% theme(axis.line = element_blank(), axis.ticks = element_blank(),
                                    plot.margin = unit(c(0, 0, 0, 0), "cm"),
                                    panel.border = element_blank())
  
  gPlot_comb <- cowplot::plot_grid(gPlot, gText, ncol = 1,
                                   rel_heights = c(7,1), rel_widths = c(1, 1),
                                   align = "v", axis = "lr")
  # Not using q-values per suggestion
  return(gPlot)
}
