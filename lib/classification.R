#' Utility functions for classification

# Perform sparce PCA and return the top features
get_top_otus_spc1 <- function(df_otu, num_dim, scale. = FALSE){
  num_EigWs <- num_dim
  
  # Run Sparse PCA
  pca_res <- nsprcomp(df_otu, k=num_EigWs, ncomp=7, center = TRUE, scale. = )
  dim_ranks <- pca_res$rotation[,"PC1"]
  
  # Identify the top features to return
  df_eigWs <- data.frame(dim_ranks)
  colnames(df_eigWs) <- c("value")
  df_eigWs$abs <- abs(df_eigWs$value)
  df_eigWs$sign <- as.factor(sign(df_eigWs$value))
  df_eigWs_ordered <- df_eigWs[order(df_eigWs$abs, decreasing = TRUE),]
  eigWs_top <- rownames(df_eigWs_ordered[c(1:num_EigWs),])
  
  return(eigWs_top)
}

# Create a binary 'Label' column based on the specified 'label_col' and 
# the associated label grouping 'label_group'.
s_get_binary_df <- function(df_comb, label_col, label_group){
  idx_true <- which(df_comb[,label_col] %in% label_group$true)
  df_sub_true <- df_comb[idx_true,]
  df_sub_true$Label <- "TRUE"
  
  idx_false <- which(df_comb[,label_col] %in% label_group$false)
  df_sub_false <- df_comb[idx_false,]
  df_sub_false$Label <- "FALSE"
  
  df_group <- rbind(df_sub_false, df_sub_true)
  df_group$Label <- factor(df_group$Label)
  
  return(df_group)
}

# Calculate the area under the ROC curve given LOO-CV results 
get_auc_roc <- function(df_loo_cv){
  roc_obj <- pROC::roc(as.numeric(df_loo_cv$Label), 
                       df_loo_cv$Score, direction="<", levels=c(0,1))
  return(roc_obj$auc[1])
}

# Calculate the area under the PR curve given LOO-CV results
get_auc_pr <- function(df_loo_cv){
  pr_obj <- pr.curve(df_loo_cv[df_loo_cv$Label == TRUE, c("Score")],
                     df_loo_cv[df_loo_cv$Label == FALSE, c("Score")],
                     curve = TRUE, rand.compute = TRUE)
  
  return(pr_obj$auc.integral)
}

# Calculate the F1-score given LOO-CV results
get_F1_score <- function(df_loo_cv, positive = TRUE){
  score_threshold <- dplyr::nth(df_loo_cv$Score, sum(df_loo_cv$Label == FALSE), order_by = df_loo_cv$Score)
  f1_score <- MLmetrics::F1_Score(df_loo_cv$Label, df_loo_cv$Score > score_threshold, positive = positive)
  
  return(f1_score)
}

# Classify using Random-Forest with LOO-CV
get_loo_CR_scores_RF <- function(x, labels){
  set.seed(RND_SEED) # to ensure reproducability
  df_res <- data.frame(Score=numeric(), Label=logical())
  df_weights <- data.frame()
  
  # LOO-CR
  for (i in c(1:nrow(x))) {
    # Create train and test
    train.x <- x[-i,,drop=FALSE]
    train.labels <- factor(labels[-i], levels = c(FALSE, TRUE)) # RF needs y to be factor 
    test.x <- x[i, ,drop=FALSE]
    test.labels <- labels[i,drop=FALSE]
    
    # Create RF model
    rf_res <- randomForest::randomForest(x=train.x, y= train.labels,
                                         ntree=RF_NTREE, maxnodes=RF_MaxNodes,
                                         strata=train.labels)
    
    # Predict for test.x
    rf_pred <- predict(rf_res, test.x, type="prob")
    
    # Append model's prediction(i.e. score) and actual Label
    df_row <- cbind(rf_pred[,c("TRUE")], test.labels)
    colnames(df_row) <- c("Score", "Label")
    
    # Append
    df_res <- rbind(df_res, df_row)
    df_weights <- rbind(df_weights, t(rf_res$importance))
  }
  df_res$Label <- as.logical(df_res$Label)
  
  # Weights
  df_weights_aggr <- data.frame(colMeans(df_weights))
  df_weights_aggr$name <- rownames(df_weights_aggr)
  colnames(df_weights_aggr) <- c("weight", "name")
  
  return( list(scores=df_res, 
               model_weights = df_weights_aggr, 
               auc = get_auc_roc(df_res),
               auc_pr = get_auc_pr(df_res),
               f1_score = get_F1_score(df_res)))
}

# Perform recurvise feature elimination using Random-Forests
select_otus_RF <- function(x, labels){
  all_res <- list()
  best_auc <- -1
  best_features <- NULL
  best_idx <- -1
  x_current <- x
  df_featues = data.frame(matrix(nrow=0, ncol = ncol(x), dimnames = list(NULL, colnames(x))))
  
  # Loop for each remaining feature
  for(i in c(1:ncol(x))){
    # Classify
    loo_CR_res <- get_loo_CR_scores_RF(x_current, labels)
    all_res[[i]] <- loo_CR_res
    
    # Update features used
    df_features_current <- matrix(rep(FALSE, ncol(x)), nrow=1, ncol = ncol(x), dimnames = list(NULL, colnames(x)))
    df_features_current[, colnames(x_current)] <- TRUE
    df_featues<- rbind(df_featues, df_features_current)
    
    # Update best
    if(loo_CR_res$auc > best_auc){
      best_idx <- i
      best_auc <- loo_CR_res$auc
      best_features <- loo_CR_res$model_weights$name
    }
    
    # Remove least important otu
    cols_new <- loo_CR_res$model_weights[-which.min(loo_CR_res$model_weights$weight), c("name")]
    all_res[[i]]$next_dropped_feature <- loo_CR_res$model_weights[which.min(loo_CR_res$model_weights$weight), c("name")]
    
    
    x_current <- x[,cols_new, drop=FALSE]
    print(sprintf("i:%d,top:%s, auc:%.3f, auc_pr:%.3f, f1:%.3f",
                  i,
                  loo_CR_res$model_weights[which.max(loo_CR_res$model_weights$weight), c("name")],
                  loo_CR_res$auc,
                  loo_CR_res$auc_pr,
                  loo_CR_res$f1_score))
  }
  
  # Summarize results
  df_summary <- data.frame(t(data.frame(lapply(all_res, function(x) {c(x$auc, x$auc_pr, x$f1_score, nrow(x$model_weights), x$next_dropped_feature)}))))
  df_summary <- data.frame(lapply(df_summary, as.character), stringsAsFactors=FALSE)
  colnames(df_summary) <- c("auc","auc_pr", "f1_score", "num_features", "next_dropped_feature")
  df_summary$auc <- as.numeric(df_summary$auc)
  df_summary$auc_pr <- as.numeric(df_summary$auc_pr)
  df_summary$f1_score <- as.numeric(df_summary$f1_score)
  df_summary$num_features <- as.integer(df_summary$num_features)
  
  rownames(df_featues) <- df_summary$num_features
  
  return(list(best_auc = best_auc, 
              best_features = best_features, 
              df_summary = df_summary,
              df_featues = df_featues))
}

# Perform recursive feeature elimination using binary Random-Forsests with LOO-CV
s_classify_binary_RF_loo_iterative_feature_removal <- function(df_comb, x_cols, label_col, label_group){
  # Get the groups
  df_group <- s_get_binary_df(df_comb, label_col, label_group)
  
  # Train
  train_x <- df_group[, x_cols, drop=FALSE]
  train_label <- df_group[, c("Label"), drop=TRUE]
  rf_res <- select_otus_RF(train_x, as.logical(train_label))
  
  return(rf_res)
}

# Binary classification with random forests in LOO-CV setting
s_classify_binary_RF_loo_cr <- function(df_comb, x_cols, label_col, label_group){
  set.seed(RND_SEED) # to ensure reproducability
  
  # Get the groups
  df_group <- s_get_binary_df(df_comb, label_col, label_group)
  
  df_res <- data.frame(Score=numeric(), Label=factor(levels= levels(df_group$Label)))
  for (i in c(1:nrow(df_group))) {
    # Create train and test
    train_x <- df_group[-i, x_cols, drop=FALSE]
    train_label <- df_group[-i, c("Label")]
    test_x <- df_group[i, x_cols, drop=FALSE]
    test_label <- df_group[i, c("Label"), drop=FALSE]
    
    # Train
    rf_res <- randomForest::randomForest(x=train_x, 
                                         y=train_label,
                                         strata=train_label,
                                         ntree=RF_NTREE, 
                                         maxnodes=RF_MaxNodes)
    
    # Test
    # Predict for test_x
    rf_pred <- predict(rf_res, test_x, type="prob")
    
    # Append model's prediction(i.e. score) and actual Label
    df_row <- cbind(rf_pred[,c("TRUE")], test_label)
    colnames(df_row) <- c("Score", "Label")
    
    # Append
    df_res <- rbind(df_res, df_row)
  }
  df_res$Label <- as.logical(df_res$Label)
  
  return(df_res)
}

# Plot Precision-Recall curve given binary classifiation results
get_pr_plot <- function(df_res, color_code, add_label = TRUE, positive = TRUE, add_star = TRUE){
  if(!positive){
    df_res$Score <- - df_res$Score
  }
  
  pr_obj <- pr.curve(df_res[df_res$Label == positive, c("Score")],
                     df_res[df_res$Label != positive, c("Score")],
                     curve = TRUE, rand.compute = TRUE)
  
  df_plot_data <- data.frame(pr_obj$curve)
  df_plot_data$Random <- FALSE
  df_plot_data <- rbind(df_plot_data, data.frame(pr_obj$rand$curve, Random=TRUE))
  
  gPlot <- ggplot(df_plot_data, aes(x=X1, y=X2, linetype=Random, colour=Random))+
    geom_line(colour=colors_assigned[color_code])+
    geom_ribbon(aes(fill=Random, ymin=0, ymax=X2), alpha=0.3)+
    scale_x_continuous("Recall", limits = c(0, 1), expand = c(0, 0), breaks = c(0, 0.5, 1))+
    scale_y_continuous("Precision", limits = c(0, 1), expand = c(0, 0), breaks = c(0, 0.5, 1))+
    scale_fill_manual("", values = c("FALSE"=unname(colors_assigned[color_code]),
                                     "TRUE"=NA))+
    scale_colour_manual("", values = c("FALSE"=unname(colors_assigned[color_code]),
                                       "TRUE"="black"))+
    coord_equal()+
    my_base_theme %+replace% theme(legend.position = "none")
  
  if(add_label){
    gPlot <- gPlot + geom_text(x=0.5, y=0.1, size = 3,
                               label=sprintf("AUC=%.2f", pr_obj$auc.integral))
  }
  
  if(add_star){
    # Identify Star threshold
    star_threshold <- dplyr::nth(df_res$Score, sum(df_res$Label != positive), order_by = df_res$Score)
    #    star_threshold <- if (!positive) -star_threshold else star_threshold
    star_point <- df_plot_data[df_plot_data$X3 == star_threshold & df_plot_data$Random == FALSE, 
                               c("X1", "X2")]
    print(sprintf("Star Recall: %.3f, Precision: %.3f", star_point$X1, star_point$X2))
    
    gPlot <- gPlot + geom_text(aes(label="\u2605", x=star_point$X1, y=star_point$X2),# nudge_y = 0.01,
                               color="red", size=3, , family = "HiraKakuProN-W3")
  }
  
  message(sprintf("** AUPR: %.3f (baseline: %.3f)", pr_obj$auc.integral, pr_obj$rand$auc.integral))
  return(gPlot)
}

# Plot Receiver Operating characteristic curve given binary classifiation results
get_roc_plot <- function(df_res, color_code, add_label = TRUE, positive = TRUE, add_star = TRUE){
  if(!positive){
    df_res$Score <- - df_res$Score
  }
  
  
  roc_obj <- roc.curve(df_res[df_res$Label == positive, c("Score")],
                       df_res[df_res$Label != positive, c("Score")],
                       curve = TRUE, rand.compute = TRUE)
  
  df_plot_data <- data.frame(roc_obj$curve)
  
  # Integrate random classifier (based on majority)
  df_plot_data$Random <- FALSE
  df_plot_data <- rbind(df_plot_data, data.frame(roc_obj$rand$curve, Random=TRUE))
  
  gPlot <- ggplot(df_plot_data, aes(x=X1, y=X2, linetype=Random, colour=Random))+
    geom_line(colour=colors_assigned[color_code])+
    geom_ribbon(aes(fill=Random, ymin=0, ymax=X2), alpha=0.3)+
    scale_x_continuous("1-Specificity", limits = c(0, 1), expand = c(0, 0), breaks = c(0, 0.5, 1))+
    scale_y_continuous("Sensitivity", limits = c(0, 1), expand = c(0, 0), breaks = c(0, 0.5, 1))+
    scale_fill_manual("", values = c("FALSE"=unname(colors_assigned[color_code]),
                                     "TRUE"=NA))+
    scale_colour_manual("", values = c("FALSE"=unname(colors_assigned[color_code]),
                                       "TRUE"="black"))+
    coord_equal()+
    my_base_theme %+replace% theme(legend.position = "none")
  
  if(add_label){
    gPlot <- gPlot + geom_text(x=0.5, y=0.1, size = 3,
                               label=sprintf("AUC=%.2f", roc_obj$auc))
  }
  
  if(add_star){
    
    # Identify Star threshold
    star_threshold <- dplyr::nth(df_res$Score, sum(df_res$Label != positive), order_by = df_res$Score)
    #    star_threshold <- if (!positive) -star_threshold else star_threshold
    
    star_point <- df_plot_data[df_plot_data$X3 == star_threshold & df_plot_data$Random == FALSE, c("X1", "X2")]
    
    gPlot <- gPlot + geom_text(aes(label="\u2605", x=star_point$X1, y=star_point$X2),# nudge_y = 0.01,
                               color="red", size=3, , family = "HiraKakuProN-W3")
  }
  
  message(sprintf("** AUROC: %.3f (baseline: %.3f)", roc_obj$auc, roc_obj$rand$auc))
  
  return(gPlot)
}

# Perform binary classification
#   'df_comb': data (dataframe)
#   'microbiome_cols': classification input
#   'binary_label_groups': binary grouping of labels
#   'friendly_names': friendly names used for microbiome_cols when plotting
#   'entity_label': label that identifies the data type (taxa, pathwas, or ga-map probes) 
#   'label_col': classification output column
plot_classification_perf <- function(df_comb, microbiome_cols, binary_label_groups, friendly_names, entity_label, label_col = "Response")
{
  # 1.B) Select Features (unsupervised)
  set.seed(RND_SEED)
  num_dim <- floor(nrow(df_comb) * DIM_REDUCTION_RATIO)
  otu_picks <- get_top_otus_spc1(df_comb[,microbiome_cols], num_dim = num_dim)
  names(otu_picks) <- friendly_names[otu_picks]
  otu_picks_labels <- friendly_names[otu_picks]
  
  # 1.C) Train Random Forests (LOO-CR), model info (all)
  rf_model <- s_classify_binary_RF_loo_iterative_feature_removal(df_comb, x_cols = otu_picks, label_col = label_col,
                                                                 label_group = binary_label_groups)
  
  # 1.D) Plot: for each classifier (which differ based on selected pathways), show: pathways, F1 measure
  num_top_features <- min(6, num_dim)
  df_summary <- rf_model$df_summary[rf_model$df_summary$num_features <= num_top_features,]
  df_summary <- df_summary[order(df_summary$num_features, decreasing = TRUE),]
  df_summary$next_dropped_feature <- factor(df_summary$next_dropped_feature, ordered = TRUE,
                                            levels = df_summary[order(df_summary$num_features, decreasing = TRUE), c("next_dropped_feature")])
  df_summary <- df_summary[order(df_summary$num_features),] # So that when there are multiple max, the simplest model is selected
  df_best <- df_summary[which.max(df_summary$f1_score), c("f1_score", "next_dropped_feature")]
  y_range <- max(df_summary$f1_score) - min(df_summary$f1_score)
  
  df_summary <- df_summary[order(df_summary$num_features, decreasing = TRUE),]
  gPlot_features <- ggplot(df_summary, aes(x=next_dropped_feature, y=f1_score, group=1))+
    geom_point()+
    geom_line()+
    geom_text(aes(label="\u2605", x=df_best$next_dropped_feature, y=df_best$f1_score), nudge_y = 0.15 * y_range,
              color="red", size=3, , family = "HiraKakuProN-W3")+
    ggtitle(sprintf("No. of cases = %d", nrow(df_comb)))+
    scale_y_continuous("F1 Score", breaks = scales::pretty_breaks(n=3))+ #expand = expansion(add = 0.2), 
    scale_x_discrete(NULL, breaks = df_summary$next_dropped_feature, 
                     labels = str_wrap(otu_picks_labels[levels(df_summary$next_dropped_feature)], width = 45))+
    geom_text(x = nrow(df_summary)/2 + 0.5,
              y = min(df_summary$f1_score),
              fontface = "bold", size = 3.3,
              label = sprintf("\nIncrementally removed %s", entity_label))+
    expand_limits(y= c(min(df_summary$f1_score) - 0.25*y_range, 
                       max(df_summary$f1_score) + 0.2*y_range))+
    my_base_theme %+replace% theme(axis.text.x=element_text(angle=90, hjust=1),
                                    panel.border = element_rect(fill=NA, size=.5),
                                    axis.line = element_line(arrow=arrow(type = "closed",
                                                                         angle = 10,
                                                                         length = unit(10, "points"))))
  print(gPlot_features)
  
  # 1.E) Plot ROC, PR curves relating to best F1 score
  df_summary <- df_summary[order(df_summary$num_features),] # So that when there are multiple max, the simplest model is selected
  best_num_features <- df_summary[which.max(df_summary$f1_score), c("num_features")]
  best_model_features <- t(rf_model$df_featues[rownames(rf_model$df_featues) == best_num_features,, drop=FALSE])
  otu_picks <- rownames(best_model_features[best_model_features == TRUE, 1, drop=FALSE])
  
  df_loo_cv <- s_classify_binary_RF_loo_cr(df_comb, x_cols = otu_picks, label_col = label_col,
                                           label_group = binary_label_groups)
  message(sprintf("F1-Score: %.3f", get_F1_score(df_loo_cv)))
  #message(sprintf("F1-Score: %.3f", get_F1_score_baseline(df_loo_cv$Label)))
  gPlot_roc <- get_roc_plot(df_loo_cv, color_code = "D", add_label = FALSE, add_star = TRUE)
  gPlot_pr <- get_pr_plot(df_loo_cv, color_code = "D", add_label = FALSE, add_star = TRUE)
  
  
  align_l <- cowplot::align_plots(gPlot_roc, gPlot_features, align = "v", axis = "l" )
  align_r <- cowplot::align_plots(gPlot_pr, align_l[[2]], align = "v", axis = "r" )
  
  gPlot_curves <- cowplot::plot_grid(align_l[[1]], align_r[[1]], nrow=1, ncol=2, align = "h")
  gPlot_ML_Pathways <- cowplot::plot_grid(gPlot_curves, align_l[[2]], nrow=2, ncol=1,
                                          rel_heights = c(1,2.2))
  
  return(list("main" = gPlot_ML_Pathways, 
              "top_left" = align_l[[1]], "top_right" = align_r[[1]],
              "bottom" = align_l[[2]]))
}