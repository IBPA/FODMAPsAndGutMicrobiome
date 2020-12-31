#!/usr/bin/env Rscript

#' Integrate GA-map microbiome data from multiple studies
#' 1) Standardize before integration
#' 2) Do not Standardize before integration

load_GA_map <- function(study_name){
  filename <- file.path("./preprocessing", study_name, "microbiome_data", "pre_diet.csv")
  df <- read.table(filename, header = TRUE, sep = ",")
  df <- cbind(data.frame(study_name=rep(study_name, nrow(df))), 
              df)

  # Set probe to zero if NA
  for(i in c(3: ncol(df))){
    if (is.na(sd(df[,i]))){
      df[,i] <- 0
    }
  }
  
  return(df)
}

# 1.A) Load data from individual studies
df_s4 <- load_GA_map("s4")
df_s5 <- load_GA_map("s5")

# 1.B) Combine data
df_comb <- rbind(df_s4, df_s5[, colnames(df_s4)])

# 1.C) Save
write.table(df_comb, file = "./integrated_analysis/integrated_GA_map/pre_diet.csv",
            sep = ",", col.names = TRUE, na = "", row.names = FALSE)
file.copy(from = file.path("./preprocessing/s4/microbiome_data/probe_names.csv"),
          to = "./integrated_analysis/integrated_GA_map/probe_names.csv",
          overwrite = TRUE)

