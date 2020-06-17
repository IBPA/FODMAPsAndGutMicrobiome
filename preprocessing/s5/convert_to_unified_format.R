# Convert microbiome data and metadata into our unified format

# 1) Load raw file
STUDY_DIR <- "./preprocessing/s5/"
RAW_DATA_DIR <- file.path(STUDY_DIR, "raw_data")
MICROBIOME_DATA_DIR <- file.path(STUDY_DIR, "microbiome_data")
src_filepath <- file.path(RAW_DATA_DIR, "all.xlsx")
df <- readxl::read_excel(src_filepath)
microbiome_cols <- colnames(df)[9:ncol(df)]

# 2) Save pre-diet microbiome
df_pre_microbiome <- df[df$LFD == "Before",
                        c("IID", microbiome_cols)]
write.table(df_pre_microbiome, file.path(MICROBIOME_DATA_DIR, "pre_diet.csv"), 
            sep = ",", row.names = FALSE)

# 3) Save post-diet microbiome
df_post_microbiome <- df[df$LFD == "After",
                         c("IID", microbiome_cols)]
write.table(df_post_microbiome, file.path(MICROBIOME_DATA_DIR, "post_diet.csv"),
            sep = ",", row.names = FALSE)

# 4) Combine pre and post metadata
df_meta_before <- df[df$LFD == "Before", c("IID", "IBS-SSS")]
df_meta_after <- df[df$LFD == "After", c("IID", "IBS-SSS")]
df_meta_comb <- merge(df_meta_before, df_meta_after, by = "IID", all = TRUE, 
                      suffixes = c("_pre", "_post"))
df_meta <- data.frame(host_id = df_meta_comb$IID,
                      pre_diet_microbiome_sample_id = df_meta_comb$IID,
                      pre_diet_IBS_SSS = df_meta_comb$`IBS-SSS_pre`,
                      post_diet_microbiome_sample_id = df_meta_comb$IID,
                      post_diet_IBS_SSS = df_meta_comb$`IBS-SSS_post`)
# 5) Save metadata
write.table(df_meta, file.path(STUDY_DIR, "metadata.csv"),
            sep = ",", row.names = FALSE)
