library(readxl)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

set.seed(125)

message("=== Loading data for binary analysis ===")
data_raw <- read_excel("5f_roc.xlsx")
feature_pool <- c("NCAN", "FIB4", "LiverRisk_score", "aMAP", "mPAGEB")
mean_fpr <- seq(0, 1, length.out = 100)
model_colors <- c("#D7191CFF", "#FDAE61FF", "#ABD9E9FF", "#2C7BB6FF", "#313695")

get_selected_features <- function(model_name) {
  switch(
    model_name,
    "NCAN" = c("NCAN"),
    "FIB4" = c("FIB4"),
    "NCAN_FIB4" = c("NCAN", "FIB4"),
    "LiverRisk_score" = c("LiverRisk_score"),
    "NCAN_LiverRisk_score" = c("NCAN", "LiverRisk_score"),
    "aMAP" = c("aMAP"),
    "NCAN_aMAP" = c("NCAN", "aMAP"),
    "mPAGEB" = c("mPAGEB"),
    "NCAN_mPAGEB" = c("NCAN", "mPAGEB")
  )
}

as_binary_01 <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.logical(x)) return(as.integer(x))
  ux <- sort(unique(x))
  ux <- ux[!is.na(ux)]
  if (length(ux) != 2) stop("Outcome is not binary (exactly 2 unique values required).")
  if (all(ux %in% c(0, 1))) return(as.integer(x))
  as.integer(x == ux[2])
}

calculate_roc <- function(X_train, Y_train, X_test, Y_test, mean_fpr) {
  fit <- glm(Y_train ~ ., data = X_train, family = binomial())
  y_pred <- predict(fit, newdata = X_test, type = "response")
  roc_curve <- roc(Y_test, y_pred, direction = "<", quiet = TRUE)
  fpr <- rev(1 - roc_curve$specificities)
  tpr <- rev(roc_curve$sensitivities)
  tpr_interp <- approx(fpr, tpr, xout = mean_fpr, ties = mean, rule = 2)$y
  tpr_interp[1] <- 0
  list(tpr_interp = tpr_interp, auc_value = as.numeric(auc(roc_curve)), y_pred = y_pred, Y_test = Y_test)
}

calculate_youden_metrics <- function(Y_true, Y_pred) {
  best_threshold <- 0
  best_youden_index <- -1
  for (threshold in seq(0, 1, by = 0.01)) {
    pred <- ifelse(Y_pred > threshold, 1, 0)
    tp <- sum((Y_true == 1) & (pred == 1))
    tn <- sum((Y_true == 0) & (pred == 0))
    fp <- sum((Y_true == 0) & (pred == 1))
    fn <- sum((Y_true == 1) & (pred == 0))
    sensitivity <- ifelse((tp + fn) == 0, NA_real_, tp / (tp + fn))
    specificity <- ifelse((tn + fp) == 0, NA_real_, tn / (tn + fp))
    youden_index <- sensitivity + specificity - 1
    if (is.finite(youden_index) && youden_index > best_youden_index) {
      best_youden_index <- youden_index
      best_threshold <- threshold
    }
  }

  pred <- ifelse(Y_pred > best_threshold, 1, 0)
  tp <- sum((Y_true == 1) & (pred == 1))
  tn <- sum((Y_true == 0) & (pred == 0))
  fp <- sum((Y_true == 0) & (pred == 1))
  fn <- sum((Y_true == 1) & (pred == 0))
  sensitivity <- ifelse((tp + fn) == 0, NA_real_, tp / (tp + fn))
  specificity <- ifelse((tn + fp) == 0, NA_real_, tn / (tn + fp))
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  precision <- ifelse((tp + fp) == 0, NA_real_, tp / (tp + fp))
  f1_score <- ifelse(is.na(precision) || is.na(sensitivity) || (precision + sensitivity) == 0, NA_real_, 2 * (precision * sensitivity) / (precision + sensitivity))
  ppv <- precision
  npv <- ifelse((tn + fn) == 0, NA_real_, tn / (tn + fn))

  list(
    sensitivity = sensitivity,
    specificity = specificity,
    accuracy = accuracy,
    precision = precision,
    f1_score = f1_score,
    youden_index = best_youden_index,
    ppv = ppv,
    npv = npv,
    best_cut_off = best_threshold
  )
}

make_cal_bin <- function(y_true, y_pred, n_groups = 5) {
  df <- data.frame(y = y_true, p = y_pred)
  df <- df[is.finite(df$p) & !is.na(df$y), , drop = FALSE]
  df$p_noise <- df$p + runif(nrow(df), 0, 1e-9)
  df$group <- dplyr::ntile(df$p_noise, n_groups)
  df %>%
    group_by(group) %>%
    summarize(pred = mean(p), obs = mean(y), n = n(), .groups = "drop")
}

calc_nri_idi_bin <- function(y_true, p_old, p_new) {
  p_new_1 <- p_new[y_true == 1]
  p_old_1 <- p_old[y_true == 1]
  p_new_0 <- p_new[y_true == 0]
  p_old_0 <- p_old[y_true == 0]
  idi <- (mean(p_new_1) - mean(p_old_1)) - (mean(p_new_0) - mean(p_old_0))
  up_1 <- sum(p_new_1 > p_old_1)
  down_1 <- sum(p_new_1 < p_old_1)
  up_0 <- sum(p_new_0 > p_old_0)
  down_0 <- sum(p_new_0 < p_old_0)
  nri_event <- (up_1 - down_1) / length(p_new_1)
  nri_nonevent <- (down_0 - up_0) / length(p_new_0)
  c(IDI = idi, NRI = nri_event + nri_nonevent)
}

calc_nri_idi_boot <- function(y_true, p_old, p_new, n_boot = 200) {
  est <- calc_nri_idi_bin(y_true, p_old, p_new)
  boot_res <- replicate(n_boot, {
    idx <- sample.int(length(y_true), replace = TRUE)
    calc_nri_idi_bin(y_true[idx], p_old[idx], p_new[idx])
  })
  idi_se <- sd(boot_res["IDI", ])
  nri_se <- sd(boot_res["NRI", ])
  idi_z <- ifelse(is.finite(idi_se) && idi_se > 0, est["IDI"] / idi_se, NA_real_)
  nri_z <- ifelse(is.finite(nri_se) && nri_se > 0, est["NRI"] / nri_se, NA_real_)
  idi_p <- ifelse(is.finite(idi_z), 2 * pnorm(-abs(idi_z)), NA_real_)
  nri_p <- ifelse(is.finite(nri_z), 2 * pnorm(-abs(nri_z)), NA_real_)
  data.frame(
    IDI = est["IDI"],
    IDI_lower = est["IDI"] - 1.96 * idi_se,
    IDI_upper = est["IDI"] + 1.96 * idi_se,
    IDI_pval = idi_p,
    NRI = est["NRI"],
    NRI_lower = est["NRI"] - 1.96 * nri_se,
    NRI_upper = est["NRI"] + 1.96 * nri_se,
    NRI_pval = nri_p
  )
}

run_binary_scenario <- function(data_raw, label_col, model_names, scenario_name, out_dir) {
  message(sprintf("\n==============================\nScenario: %s\nOutcome: %s\n==============================", scenario_name, label_col))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  data <- data_raw %>% filter(!is.na(.data[[label_col]]))
  data <- data[, c(label_col, feature_pool)]
  data <- data[complete.cases(data), , drop = FALSE]
  label <- as_binary_01(data[[label_col]])
  x_data <- data[, feature_pool, drop = FALSE]

  message(sprintf("Sample size after filtering: %d", nrow(data)))
  message("Creating 5-fold CV...")
  folds <- createFolds(label, k = 5, list = TRUE, returnTrain = TRUE)

  tprs <- list()
  aucs <- list()
  model_metrics <- list()
  model_predictions <- list()

  for (model_name in model_names) {
    message(sprintf("Running model: %s", model_name))
    selected_features <- get_selected_features(model_name)
    model_tprs <- list()
    model_aucs <- c()
    all_predicted_probs <- c()
    all_true_labels <- c()
    model_predictions[[model_name]] <- list(y_pred = c(), Y_test = c())

    for (i in seq_along(folds)) {
      train_index <- folds[[i]]
      test_index <- setdiff(seq_along(label), train_index)
      X_train <- x_data[train_index, selected_features, drop = FALSE]
      Y_train <- label[train_index]
      X_test <- x_data[test_index, selected_features, drop = FALSE]
      Y_test <- label[test_index]

      roc_results <- calculate_roc(X_train, Y_train, X_test, Y_test, mean_fpr)
      model_aucs <- c(model_aucs, roc_results$auc_value)
      model_tprs[[i]] <- roc_results$tpr_interp
      all_predicted_probs <- c(all_predicted_probs, roc_results$y_pred)
      all_true_labels <- c(all_true_labels, Y_test)
      model_predictions[[model_name]]$y_pred <- c(model_predictions[[model_name]]$y_pred, roc_results$y_pred)
      model_predictions[[model_name]]$Y_test <- c(model_predictions[[model_name]]$Y_test, Y_test)
    }

    mean_auc <- mean(model_aucs)
    mean_tpr <- rowMeans(do.call(cbind, model_tprs))
    std_tpr <- apply(do.call(cbind, model_tprs), 1, sd)
    aucs[[model_name]] <- mean_auc
    tprs[[model_name]] <- list(mean_tpr = mean_tpr, std_tpr = std_tpr)

    metrics <- calculate_youden_metrics(all_true_labels, all_predicted_probs)
    model_metrics[[model_name]] <- c(
      mean_auc,
      metrics$sensitivity,
      metrics$specificity,
      metrics$accuracy,
      metrics$precision,
      metrics$f1_score,
      metrics$youden_index,
      metrics$ppv,
      metrics$npv,
      metrics$best_cut_off
    )

    message(sprintf("  AUC = %.4f | Sens = %.4f | Spec = %.4f", mean_auc, metrics$sensitivity, metrics$specificity))
  }

  model_metrics_df <- do.call(rbind, lapply(model_metrics, function(x) as.data.frame(t(x))))
  colnames(model_metrics_df) <- c("AUC", "Sensitivity", "Specificity", "Accuracy", "Precision", "F1_score", "Youden_index", "PPV", "NPV", "Best_cutoff")
  model_metrics_df <- cbind(Model = rownames(model_metrics_df), model_metrics_df)
  row.names(model_metrics_df) <- NULL
  print(model_metrics_df)
  write.csv(model_metrics_df, file.path(out_dir, "model_metrics.csv"), row.names = FALSE)

  roc_objects <- list()
  for (model_name in names(model_predictions)) {
    roc_objects[[model_name]] <- roc(model_predictions[[model_name]]$Y_test, model_predictions[[model_name]]$y_pred, quiet = TRUE)
  }

  delong_df <- data.frame()
  model_names_loop <- names(model_predictions)
  for (i in 1:(length(model_names_loop) - 1)) {
    for (j in (i + 1):length(model_names_loop)) {
      model1 <- model_names_loop[i]
      model2 <- model_names_loop[j]
      test_result <- roc.test(roc_objects[[model1]], roc_objects[[model2]])
      row <- data.frame(
        model_1 = model1,
        model_2 = model2,
        z = as.numeric(test_result$statistic),
        p_value = as.numeric(test_result$p.value)
      )
      delong_df <- bind_rows(delong_df, row)
      message(sprintf("Delong %s vs %s: Z=%.4f, p=%.4g", model1, model2, row$z, row$p_value))
    }
  }
  write.csv(delong_df, file.path(out_dir, "delong_pairwise.csv"), row.names = FALSE)
  print(delong_df)

  p_roc <- ggroc(roc_objects, aes = "colour", linewidth = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(x = "False Positive Rate", y = "True Positive Rate", color = "Model", title = paste0(scenario_name, " ROC (OOF)")) +
    scale_color_manual(
      values = setNames(model_colors[seq_along(model_names_loop)], model_names_loop),
      labels = sapply(model_names_loop, function(name) paste0(name, " (AUC = ", round(aucs[[name]], 3), ")"))
    ) +
    theme_bw(base_size = 12)
  print(p_roc)
  ggsave(file.path(out_dir, "roc_ggroc.png"), p_roc, width = 9, height = 6, dpi = 300)

  plot_data <- data.frame(
    mean_fpr = rep(mean_fpr, length(model_names_loop)),
    mean_tpr = unlist(lapply(tprs, function(x) x$mean_tpr)),
    model_name = rep(model_names_loop, each = length(mean_fpr))
  )
  plot_data$model_name <- factor(plot_data$model_name, levels = model_names_loop)
  p_roc_mean <- ggplot(plot_data, aes(x = mean_fpr, y = mean_tpr, color = model_name, group = model_name)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = model_colors[seq_along(model_names_loop)]) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(x = "False Positive Rate", y = "True Positive Rate", title = paste0(scenario_name, " ROC (Mean Curve)")) +
    theme_bw(base_size = 12)
  print(p_roc_mean)
  ggsave(file.path(out_dir, "roc_mean_curve.png"), p_roc_mean, width = 9, height = 6, dpi = 300)

  message("Calculating OOF-based Brier/Calibration/DCA/NRI/IDI...")
  oof_df <- data.frame(Y_true = model_predictions[[model_names_loop[1]]]$Y_test)
  for (m in model_names_loop) oof_df[[m]] <- model_predictions[[m]]$y_pred
  write.csv(oof_df, file.path(out_dir, "oof_predictions.csv"), row.names = FALSE)
  print(head(oof_df))

  brier_all <- data.frame(Model = model_names_loop, Brier = NA_real_)
  for (i in seq_along(model_names_loop)) {
    m <- model_names_loop[i]
    brier_all$Brier[i] <- mean((oof_df$Y_true - oof_df[[m]])^2)
  }
  print(brier_all)
  write.csv(brier_all, file.path(out_dir, "brier_score.csv"), row.names = FALSE)

  cal_list <- list()
  for (m in model_names_loop) {
    cal_tmp <- make_cal_bin(oof_df$Y_true, oof_df[[m]], n_groups = 6)
    cal_tmp$Model <- m
    cal_list[[m]] <- cal_tmp
  }
  cal_all <- bind_rows(cal_list)
  print(cal_all)
  write.csv(cal_all, file.path(out_dir, "calibration_data.csv"), row.names = FALSE)
  p_cal <- ggplot(cal_all, aes(x = pred, y = obs, color = Model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2) +
    scale_color_manual(values = setNames(model_colors[seq_along(model_names_loop)], model_names_loop)) +
    labs(title = paste0("Calibration Curve (", scenario_name, ", OOF)"), x = "Predicted Probability", y = "Observed Proportion") +
    theme_bw(base_size = 12)
  print(p_cal)
  ggsave(file.path(out_dir, "calibration_curve.png"), p_cal, width = 9, height = 6, dpi = 300)

  risk_values <- unlist(oof_df[, model_names_loop, drop = FALSE], use.names = FALSE)
  q <- suppressWarnings(quantile(risk_values, probs = c(0.05, 0.95), na.rm = TRUE, names = FALSE))
  thr_low <- max(0.01, min(0.99, q[1]))
  thr_high <- max(0.01, min(0.99, q[2]))
  if (thr_high - thr_low < 0.02) {
    thresholds <- seq(0.01, 0.99, by = 0.01)
  } else {
    thresholds <- unique(round(seq(thr_low, thr_high, by = 0.01), 2))
  }

  N <- nrow(oof_df)
  positive_n <- sum(oof_df$Y_true == 1, na.rm = TRUE)
  prev <- positive_n / N
  message(sprintf("DCA prevalence (positive/total): %d/%d = %.4f", positive_n, N, prev))
  dca_prevalence <- data.frame(
    scenario = scenario_name,
    positive_n = positive_n,
    total_n = N,
    positive_ratio = prev
  )
  print(dca_prevalence)
  write.csv(dca_prevalence, file.path(out_dir, "dca_prevalence.csv"), row.names = FALSE)
  dca_list <- list()
  for (m in model_names_loop) {
    nb <- sapply(thresholds, function(pt) {
      pred_pos <- oof_df[[m]] >= pt
      tp <- sum(oof_df$Y_true == 1 & pred_pos)
      fp <- sum(oof_df$Y_true == 0 & pred_pos)
      (tp / N) - (fp / N) * (pt / (1 - pt))
    })
    dca_list[[m]] <- data.frame(threshold = thresholds, net_benefit = nb, Model = m)
  }
  dca_list[["Treat all"]] <- data.frame(
    threshold = thresholds,
    net_benefit = sapply(thresholds, function(pt) prev - (1 - prev) * (pt / (1 - pt))),
    Model = "Treat all"
  )
  dca_list[["Treat none"]] <- data.frame(threshold = thresholds, net_benefit = 0, Model = "Treat none")
  dca_all <- bind_rows(dca_list)
  print(head(dca_all, 20))
  write.csv(dca_all, file.path(out_dir, "dca_data.csv"), row.names = FALSE)
  dca_colors <- c(
    setNames(model_colors[seq_along(model_names_loop)], model_names_loop),
    "Treat all" = "grey30",
    "Treat none" = "black"
  )
  p_dca <- ggplot(dca_all, aes(x = threshold, y = net_benefit, color = Model)) +
    geom_line(linewidth = 1.1) +
    coord_cartesian(ylim = c(0, max(dca_all$net_benefit, na.rm = TRUE) * 1.05)) +
    scale_color_manual(values = dca_colors) +
    labs(title = paste0("Decision Curve (", scenario_name, ", OOF)"), x = "Threshold Probability", y = "Net Benefit") +
    theme_bw(base_size = 12)
  print(p_dca)
  ggsave(file.path(out_dir, "dca_curve.png"), p_dca, width = 9, height = 6, dpi = 300)

  comparisons <- list(
    c("FIB4", "NCAN_FIB4"),
    c("LiverRisk_score", "NCAN_LiverRisk_score"),
    c("aMAP", "NCAN_aMAP"),
    c("mPAGEB", "NCAN_mPAGEB")
  )
  nri_idi_res <- data.frame()
  for (pair in comparisons) {
    m_old <- pair[1]
    m_new <- pair[2]
    if (m_old %in% model_names_loop && m_new %in% model_names_loop) {
      message(sprintf("NRI/IDI comparing %s vs %s ...", m_new, m_old))
      res <- calc_nri_idi_boot(oof_df$Y_true, oof_df[[m_old]], oof_df[[m_new]], n_boot = 200)
      res$Base_Model <- m_old
      res$New_Model <- m_new
      nri_idi_res <- bind_rows(nri_idi_res, res)
    }
  }
  if (nrow(nri_idi_res) > 0) {
    print(nri_idi_res)
    write.csv(nri_idi_res, file.path(out_dir, "nri_idi_summary.csv"), row.names = FALSE)
  }

  message(sprintf("[SUCCESS] Scenario '%s' completed. Output directory: %s", scenario_name, out_dir))
}

scenario_configs <- list(
  list(
    scenario_name = "HCC",
    label_col = "HCC",
    model_names = c("NCAN", "aMAP", "NCAN_aMAP", "mPAGEB", "NCAN_mPAGEB"),
    out_dir = "model_compare_outputs_binary_HCC"
  ),
  list(
    scenario_name = "FIB",
    label_col = "FIB",
    model_names = c("NCAN", "FIB4", "NCAN_FIB4", "LiverRisk_score", "NCAN_LiverRisk_score"),
    out_dir = "model_compare_outputs_binary_FIB"
  )
)

message("=== Start running two binary scenarios ===")
for (cfg in scenario_configs) {
  run_binary_scenario(
    data_raw = data_raw,
    label_col = cfg$label_col,
    model_names = cfg$model_names,
    scenario_name = cfg$scenario_name,
    out_dir = cfg$out_dir
  )
}
message("=== All scenarios completed ===")
