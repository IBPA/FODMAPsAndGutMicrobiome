#' Main configuration file

STUDIES_DIR <- "preprocessing"
STUDY_CONFIG_FILENAME <- "config.ini"
STUDY_METADATA_FILENAME <- "metadata.csv"
STUDY_IDS <- c("s1", "s2", "s3", "s4", "s5")

# Maximum improvement in IBS-SSS score for the "No" response group 
No_Response_Max_Improvement <- 22
# Maximum improvement in IBS-SSS score for the "Low" response group
Low_Response_Max_Improvement <- 150

######## configuration management methods #######

load_study_config <- function(id){
  config_filepath <- file.path(STUDIES_DIR, id, STUDY_CONFIG_FILENAME)
  config <- ini::read.ini(config_filepath)
  config$info$metadata_filepath <- file.path(STUDIES_DIR, id, STUDY_METADATA_FILENAME)
  
  return(config$info)
}

load_study_configs <- function() {
  res <- list()
  for(id in STUDY_IDS){
    res[[id]] <- load_study_config(id)
  }
  
  return(res)
}
