#' Utility functions for metadata management

combine_dataframes <- function(df_main, df_new){
  if(is.null(df_main) || nrow(df_main) == 0) 
    return (df_new)
  
  # Append to df_new if in df_main only
  new_cols_df_new <- colnames(df_main)[!(colnames(df_main) %in% colnames(df_new))]
  if (length(new_cols_df_new) > 0) {
    df_new[,new_cols_df_new] <- NA
  }
  
  # Append to df_main if in df_new only
  new_cols_df_main <- colnames(df_new)[!(colnames(df_new) %in% colnames(df_main))]
  if (length(new_cols_df_main) > 0) {
    df_main[,new_cols_df_main] <- NA
  }
  
  return(rbind(df_main, df_new[,colnames(df_main)]))
}

add_diet_response <- function(df){
  df_res <- df
  df$delta_IBS_SSS <- (df$pre_diet_IBS_SSS - df$post_diet_IBS_SSS)
  df_res$Response <- cut(df$delta_IBS_SSS,
                     breaks = c(-Inf, 
                                No_Response_Max_Improvement,
                                Low_Response_Max_Improvement, 
                                Inf),
                     labels = c("No", "Low", "High"))
  return(df_res)
}

load_study_metadata <- function(study_id, config){
  df_meta <- read.table(config$metadata_filepath, sep = ",", header = TRUE, 
                        na.strings = "")
  df_meta$study_name <- study_id
  df_meta$technology <- config$Technology
  df_meta$reference <- config$Reference
  df_meta$pre_mic_avail <- !is.na(df_meta$pre_diet_microbiome_sample_id)
  df_meta$host_id <- sprintf("%s_%s", df_meta$study_name, df_meta$host_id)
  
  return(df_meta)
}

load_metadata <- function(studies = STUDY_IDS){
  # 1) Combine metadata from all studies
  df_comb_metadata <- data.frame()
  study_configs <- load_study_configs(studies)
  for(sid in names(study_configs)){
    df_study_metadata <- load_study_metadata(sid, study_configs[[sid]])
    df_comb_metadata <- combine_dataframes(df_comb_metadata, df_study_metadata)
  }

  # 2) Identify response labels (No/Low/High)
  df_comb_metadata <- add_diet_response(df_comb_metadata)
  
  # 3) Return
  return(df_comb_metadata)
}

load_metadata_diet_all_16S <- function(df_metadata){
  df_metadata <- load_metadata(STUDY_IDS_16S)
  colnames_res <- c("name", "study_name", "host_id", "low_fodmap")
  
  # Pre:
  df_pre <- df_metadata[,c("pre_diet_microbiome_sample_id",
                           "study_name", "host_id")]
  df_pre$low_fodmap <- "FALSE"
  colnames(df_pre) <- colnames_res
  
  # Post:
  df_post <- df_metadata[,c("post_diet_microbiome_sample_id",
                            "study_name", "host_id")]
  df_post$low_fodmap <- "TRUE"
  colnames(df_post) <- colnames_res
  
  # Comb
  df_res <- rbind(df_pre, df_post)
  df_res <- df_res[!is.na(df_res$name),]
  rownames(df_res) <- df_res$name
  
  return(df_res)
}