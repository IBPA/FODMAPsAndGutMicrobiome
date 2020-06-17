# Convert microbiome data and metadata into our unified format

fix_col_names <- function(df){
  colnames(df) <- c(as.character(df[1,c(1:10)]),
                    colnames(df)[11:ncol(df)])
  
  return(df[c(2:nrow(df)),])
}

# 1) Load raw file
STUDY_DIR <- "./preprocessing/s4/"
RAW_DATA_DIR <- file.path(STUDY_DIR, "raw_data")
MICROBIOME_DATA_DIR <- file.path(STUDY_DIR, "microbiome_data")

src_filepath <- file.path(RAW_DATA_DIR, 
                          "IBS_SSS and Microbiota Data 20191125.xlsx")
df <- readxl::read_excel(src_filepath)

# 2) Fix column names
df <- fix_col_names(df)

# 3) Filter low-FODMAP rows
df <- df[!is.na(df$DIET) & df$DIET == "low-FODMAP diet",]

# 4) Fix column types
microbiome_cols <- colnames(df)[11:ncol(df)]
df[microbiome_cols] <- sapply(df[microbiome_cols], as.numeric)

# 5) Save pre-diet microbiome
df_pre_microbiome <- df[df$Intervention == "Before intervention",
                        c("ID", microbiome_cols)]
write.table(df_pre_microbiome, file.path(MICROBIOME_DATA_DIR, "pre_diet.csv"), 
            sep = ",", row.names = FALSE)

# 6) Save post-diet microbiome
df_post_microbiome <- df[df$Intervention == "After intervention",
                        c("ID", microbiome_cols)]
write.table(df_post_microbiome, file.path(MICROBIOME_DATA_DIR, "post_diet.csv"),
            sep = ",", row.names = FALSE)

# 7) Save metadata
# host_id,pre_diet_microbiome_sample_id,pre_diet_IBS_SSS,post_diet_microbiome_sample_id,post_diet_IBS_SSS
df_meta <- df[df$Intervention == "Before intervention",
              c("ID", "ID", "IBS-SSS_Day_0", "ID", "IBS-SSS_Day_29")]
df_meta[c("IBS-SSS_Day_0", "IBS-SSS_Day_29")] <- 
  sapply(df_meta[c("IBS-SSS_Day_0", "IBS-SSS_Day_29")], as.numeric)
colnames(df_meta) <- c("host_id", 
                       "pre_diet_microbiome_sample_id", 
                       "pre_diet_IBS_SSS", 
                       "post_diet_microbiome_sample_id",
                       "post_diet_IBS_SSS")

write.table(df_meta, file.path(STUDY_DIR, "metadata.csv"),
            sep = ",", row.names = FALSE)
