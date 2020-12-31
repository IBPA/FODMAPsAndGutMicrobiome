# Convert metadata into our unified format

update_group <- function(df_sample_metadata, df_seq_file_info){
  # This function updates group info (I/II) from df_sample_metadata as it appears incorrect.
  df_sample_metadata$GROUP <- NA
  # Group I:
  host_ids_g1 <- unique(df_seq_file_info[!is.na(df_seq_file_info$host_diet) 
                                         & df_seq_file_info$host_diet == "LowFODMAP", c("host_id")])
  df_sample_metadata[df_sample_metadata$host_id %in% host_ids_g1, c("GROUP") ] <- 1
  
  # Group II:
  host_ids_g2 <- unique(df_seq_file_info[!is.na(df_seq_file_info$host_diet) 
                                         & df_seq_file_info$host_diet == "Control", c("host_id")])
  df_sample_metadata[df_sample_metadata$host_id %in% host_ids_g2, c("GROUP") ] <- 2
  
  return(df_sample_metadata)
}

STUDY_DIR <- "./preprocessing/s1/"
RAW_DATA_DIR <- file.path(STUDY_DIR, "raw_data")

# Load seq_file_info.csv:
df_seq_file_info <- read.table(file.path(RAW_DATA_DIR, "seq_file_info.csv"), sep = ",", header = TRUE)
df_seq_file_info$ID <- unlist(lapply(df_seq_file_info$alias, function(x) str_split_fixed(x, "_", 2)[2]))

# Load sample_metadata.csv
df_sample_metadata <- read.table(file.path(RAW_DATA_DIR, "sample_metadata.csv"), sep = ",", header = TRUE)
df_sample_metadata$TimePoint <- substr(df_sample_metadata$ID, 0, 2)
df_sample_metadata$host_id <- as.numeric(substr(df_sample_metadata$ID, 3, 99999999))

# Fix GROUP attribute:
df_sample_metadata <- update_group(df_sample_metadata, df_seq_file_info)

# Assign DietStatus attribute
df_sample_metadata$DietStatus <- NA
df_sample_metadata[!is.na(df_sample_metadata$GROUP) & (df_sample_metadata$GROUP == 1 & df_sample_metadata$TimePoint == "BS"), c("DietStatus")] <- "RightBeforeLowFODMAP"
df_sample_metadata[!is.na(df_sample_metadata$GROUP) & (df_sample_metadata$GROUP == 1 & df_sample_metadata$TimePoint == "3M"), c("DietStatus")] <- "RightAfterLowFODMAP"
df_sample_metadata[!is.na(df_sample_metadata$GROUP) & (df_sample_metadata$GROUP == 2 & df_sample_metadata$TimePoint == "3M"), c("DietStatus")] <- "RightBeforeLowFODMAP"
df_sample_metadata[!is.na(df_sample_metadata$GROUP) & (df_sample_metadata$GROUP == 2 & df_sample_metadata$TimePoint == "6M"), c("DietStatus")] <- "RightAfterLowFODMAP"
df_sample_metadata[!is.na(df_sample_metadata$GROUP) & (df_sample_metadata$GROUP == 2 & df_sample_metadata$TimePoint == "BS"), c("DietStatus")] <- "LongBeforeLowFODMAP"

# New Dataframe
host_ids <- unique(c(df_seq_file_info$host_id, df_sample_metadata$host_id))
host_ids <- host_ids[!is.na(host_ids)]
df_new <- data.frame(matrix(ncol = 1, nrow = length(host_ids)))
colnames(df_new) <- c("host_id")
df_new$host_id <- host_ids

# Fill in pre_diet_microbiome_sample_id, pre_diet_IBS_SSS
pre_treatment_features <- c("IBS_SUMSCORE", "AGE", "FORM_NORMAL", "FORM_HARD", "FORM_STRING", "FORM_PELLETS", "FORM_MUSHY", "FORM_WATERY")
df_tmp1 <- df_seq_file_info[df_seq_file_info$misc_param == "NotThawed", c("ID", "name")]
df_tmp2 <- df_sample_metadata[!is.na(df_sample_metadata$DietStatus) & df_sample_metadata$DietStatus == "RightBeforeLowFODMAP", 
                              c("host_id", "ID", pre_treatment_features)]
df_tmp3 <- merge(df_tmp1, df_tmp2, by = "ID", all.y = TRUE)
df_new <- merge(df_new, df_tmp3[c("host_id", "name", pre_treatment_features)], by = "host_id", all.x = TRUE)
df_new$pre_diet_microbiome_sample_id <- df_new$name
df_new$pre_diet_IBS_SSS <- df_new$IBS_SUMSCORE
df_new$name <- NULL
df_new$IBS_SUMSCORE <- NULL

# If pre_diet_microbiome_sample_id is NA, then use relevent sample from record with LongBeforeLowFODMAP DietStatus
host_ids <- df_new[is.na(df_new$pre_diet_microbiome_sample_id), c("host_id")]
df_tmp1 <- df_seq_file_info[df_seq_file_info$host_id %in% host_ids  &
                              df_seq_file_info$misc_param == "NotThawed", c("ID", "name")]
df_tmp2 <- df_sample_metadata[!is.na(df_sample_metadata$DietStatus) &
                                df_sample_metadata$DietStatus == "LongBeforeLowFODMAP" &
                                df_sample_metadata$host_id %in% host_ids, 
                              c("host_id", "ID")]
df_tmp3 <- merge(df_tmp1, df_tmp2, by = "ID", all.y = TRUE)
df_new <- merge(df_new, df_tmp3[c("host_id", "name")], by = "host_id", all.x = TRUE)
df_new[is.na(df_new$pre_diet_microbiome_sample_id), c("pre_diet_microbiome_sample_id")] <- df_new[is.na(df_new$pre_diet_microbiome_sample_id), c("name")]
df_new$name <- NULL


# Fill in post_diet_microbiome_sample_id, post_diet_IBS_SSS
df_tmp1 <- df_seq_file_info[df_seq_file_info$misc_param == "NotThawed", c("ID", "name")]
df_tmp2 <- df_sample_metadata[!is.na(df_sample_metadata$DietStatus) & df_sample_metadata$DietStatus == "RightAfterLowFODMAP", 
                              c("host_id", "ID", "IBS_SUMSCORE")]
df_tmp3 <- merge(df_tmp1, df_tmp2, by = "ID", all.y = TRUE)
df_new <- merge(df_new, df_tmp3[c("host_id", "name", "IBS_SUMSCORE")], by = "host_id", all.x = TRUE)
df_new$post_diet_microbiome_sample_id <- df_new$name
df_new$post_diet_IBS_SSS <- df_new$IBS_SUMSCORE
df_new$name <- NULL
df_new$IBS_SUMSCORE <- NULL
df_new$misc_param <- NULL

df_new$IBS_subtype = NA
df_new[!is.na(df_new$FORM_HARD) & !is.na(df_new$FORM_WATERY) &
         df_new$FORM_HARD < 2 & df_new$FORM_WATERY > 1, c("IBS_subtype")] <- "IBS-D"
df_new[!is.na(df_new$FORM_HARD) & !is.na(df_new$FORM_WATERY) &
         df_new$FORM_HARD > 1 & df_new$FORM_WATERY < 2, c("IBS_subtype")] <- "IBS-C"
df_new[!is.na(df_new$FORM_HARD) & !is.na(df_new$FORM_WATERY) &
         df_new$FORM_HARD > 1 & df_new$FORM_WATERY > 1, c("IBS_subtype")] <- "IBS-M"

# Save into metadata.csv
df_new <- df_new[c("host_id", "pre_diet_microbiome_sample_id", "pre_diet_IBS_SSS", 
                   "post_diet_microbiome_sample_id", "post_diet_IBS_SSS", "IBS_subtype")]
# write.table(df_new[!is.na(df_new$pre_diet_IBS_SSS - df_new$post_diet_IBS_SSS),],
write.table(df_new, file = file.path(STUDY_DIR, "metadata.csv"), 
            sep = ",", quote = FALSE, na = "", row.names = FALSE)
print("Saved into metadata.csv")
