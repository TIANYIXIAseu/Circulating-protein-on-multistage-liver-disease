library(caret)
library(dplyr)
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(survival)
library(timeROC)
library(riskRegression)
library(rmda)
library(survIDINRI)


#### timeROC 5fold ####

# ================= Configuration =================
# Choose which outcome to run: "liver_disease" or "liver_death"
TARGET_OUTCOME <- "liver_disease" 

if (TARGET_OUTCOME == "liver_disease") {
  time_col <- "Any_liver_disease_time"
  status_col <- "Any_liver_disease"
  out_dir_prefix <- "results_liver_disease"
} else if (TARGET_OUTCOME == "liver_death") {
  time_col <- "Time_of_death"
  status_col <- "Cause of liver death: ICD10 | Instance 0"
  out_dir_prefix <- "results_liver_death"
} else {
  stop("Invalid TARGET_OUTCOME")
}

message("=== Running analysis for: ", TARGET_OUTCOME, " ===")
message("Time column: ", time_col)
message("Status column: ", status_col)

df <- read_excel("disposal data0107.xlsx", sheet = "Sheet3")
df <- as.data.frame(df)
df <- df[,c("Participant ID", "NCAN;Neurocan core protein", "FIB-4", "liverRisk score", "aMAP", "mPAGE-B", "Time_of_death", "Cause of liver death: ICD10 | Instance 0", "Any_liver_disease", "Any_liver_disease_time")]
names(df)

# Drop missing or <=0 time values for the target outcome
df <- df[!is.na(df[[time_col]]) & df[[time_col]] > 0, ]
df <- df[!is.na(df[[status_col]]), ]

df$`Time_of_death` <- df$`Time_of_death` / 365.25
df$`Any_liver_disease_time` <- df$`Any_liver_disease_time` / 365.25
time_points <- seq(1, 16)
analysis_df <- df

message("--- Starting 5-fold Cross-Validation for TimeROC ---")
set.seed(123) 
folds_cv <- createFolds(analysis_df[[status_col]], k = 5, list = TRUE)

cross_validate_roc <- function(df, time_points, folds, t_col, s_col) {
  roc_results <- list()
  
  for (fold_idx in seq_along(folds)) {
    fold <- folds[[fold_idx]]
    message(sprintf("Processing fold %d/%d for TimeROC...", fold_idx, length(folds)))
    train_data <- df[-fold, ]
    test_data <- df[fold, ]
    
    formula_ncan <- as.formula(paste0("Surv(`", t_col, "`, `", s_col, "`) ~ `NCAN;Neurocan core protein`"))
    formula_ncan_fib4 <- as.formula(paste0("Surv(`", t_col, "`, `", s_col, "`) ~ `NCAN;Neurocan core protein` + `FIB-4`"))
    formula_ncan_liverrisk <- as.formula(paste0("Surv(`", t_col, "`, `", s_col, "`) ~ `NCAN;Neurocan core protein` + `liverRisk score`"))
    formula_ncan_amap <- as.formula(paste0("Surv(`", t_col, "`, `", s_col, "`) ~ `NCAN;Neurocan core protein` + `aMAP`"))
    formula_ncan_mpageb <- as.formula(paste0("Surv(`", t_col, "`, `", s_col, "`) ~ `NCAN;Neurocan core protein` + `mPAGE-B`"))

    cox_model_NCAN <- coxph(formula_ncan, data=train_data)
    cox_model_NCAN_FIB4 <- coxph(formula_ncan_fib4, data=train_data)
    cox_model_NCAN_Liverrisk <- coxph(formula_ncan_liverrisk, data=train_data)
    cox_model_NCAN_aMAP <- coxph(formula_ncan_amap, data=train_data)
    cox_model_NCAN_mPAGEB <- coxph(formula_ncan_mpageb, data=train_data)
    
    test_data$NCAN <- predict(cox_model_NCAN, newdata=test_data, type="risk")
    test_data$NCAN_FIB4 <- predict(cox_model_NCAN_FIB4, newdata=test_data, type="risk")
    test_data$NCAN_Liverrisk <- predict(cox_model_NCAN_Liverrisk, newdata=test_data, type="risk")
    test_data$NCAN_aMAP <- predict(cox_model_NCAN_aMAP, newdata=test_data, type="risk")
    test_data$NCAN_mPAGEB <- predict(cox_model_NCAN_mPAGEB, newdata=test_data, type="risk")
    
    roc_list <- list()
    
    for (marker in c("NCAN", "FIB-4", "NCAN_FIB4", "liverRisk score", "NCAN_Liverrisk", "aMAP", "NCAN_aMAP", "mPAGE-B", "NCAN_mPAGEB")) {
      roc_res <- timeROC(T=test_data[[t_col]],
                         delta=test_data[[s_col]],
                         marker=test_data[[marker]],
                         cause=1,
                         times=time_points,
                         iid=FALSE)
      roc_list[[marker]] <- roc_res$AUC
    }
    
    roc_results[[length(roc_results) + 1]] <- roc_list
  }
  return(roc_results)
}

auroc_values <- cross_validate_roc(analysis_df, time_points, folds_cv, time_col, status_col)
message("--- TimeROC Cross-Validation Complete ---")

final_df <- data.frame(Time = integer(), fold = integer(), auc = numeric(), Model = character())

for (fold_idx in 1:length(auroc_values)) {
  fold_data <- auroc_values[[fold_idx]]
  for (model_idx in 1:length(fold_data)) {
    model_data <- fold_data[[model_idx]]
    temp_df <- data.frame(
      Time = 1:length(model_data),  
      fold = fold_idx,              
      auc = model_data,             
      Model = names(auroc_values[[fold_idx]])[model_idx]
    )
    final_df <- rbind(final_df, temp_df)
  }
}

df <- final_df

df_stats <- df %>%
  group_by(Time, Model) %>%
  summarise(
    mean_auc = mean(auc),
    sd_auc = sd(auc),
    min_auc = min(auc),
    max_auc = max(auc),
    .groups = 'drop'
  )

model_order <- c("NCAN", "FIB-4", "NCAN_FIB4", "liverRisk score", 
                 "NCAN_Liverrisk", "aMAP", "NCAN_aMAP", "mPAGE-B", "NCAN_mPAGEB")
df$Model <- factor(df$Model, levels = model_order)
df_stats$Model <- factor(df_stats$Model, levels = model_order)

df <- df %>%
  mutate(Model = case_when(
    Model == "NCAN_FIB4" ~ "NCAN_FIB-4",
    Model == "liverRisk score" ~ "LiverRisk score",
    Model == "NCAN_Liverrisk" ~ "NCAN_LiverRisk score",
    Model == "NCAN_mPAGEB" ~ "NCAN_mPAGE-B",
    TRUE ~ Model
  ))
df_stats <- df_stats %>%
  mutate(Model = case_when(
    Model == "NCAN_FIB4" ~ "NCAN_FIB-4",
    Model == "liverRisk score" ~ "LiverRisk score",
    Model == "NCAN_Liverrisk" ~ "NCAN_LiverRisk score",
    Model == "NCAN_mPAGEB" ~ "NCAN_mPAGE-B",
    TRUE ~ Model
  ))

model_order <- c("NCAN", "FIB-4", "NCAN_FIB-4", "LiverRisk score", 
                 "NCAN_LiverRisk score", "aMAP", "NCAN_aMAP", "mPAGE-B", "NCAN_mPAGE-B")
df$Model <- factor(df$Model, levels = model_order)
df_stats$Model <- factor(df_stats$Model, levels = model_order)

colors <- brewer.pal(9, "RdYlBu")

ggplot() +
  geom_point(data = df, aes(x = factor(Time), y = auc, color = Model), 
             size = 2, alpha = 0.8, position = position_jitterdodge(jitter.width = 0, dodge.width = 1)) +
  geom_errorbar(data = df_stats, 
                aes(x = factor(Time), ymin = min_auc, ymax = max_auc, color = Model), 
                width = 0.4, position = position_dodge(1)) +
  geom_rect(data = df_stats,
            aes(xmin = as.numeric(factor(Time)) - 0.35, 
                xmax = as.numeric(factor(Time)) + 0.35,
                ymin = mean_auc - sd_auc,
                ymax = mean_auc + sd_auc,
                fill = Model),
            alpha = 0.4, position = position_dodge(1)) +
  geom_point(data = df_stats, 
             aes(x = factor(Time), y = mean_auc, fill = Model), 
             shape = 21, size = 3.5, 
             position = position_dodge(1)) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.2))) + 
  theme_minimal() +
  labs(x = "Time-point (Years)", y = "AUC", color = "Model", fill = "Model") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white")
  )

if (!dir.exists(out_dir_prefix)) dir.create(out_dir_prefix, recursive = TRUE)

write.csv(df_stats, file.path(out_dir_prefix, '5 fold_auc_summary.csv'), row.names = F)
message("Saved 5-fold AUC summary to '5 fold_auc_summary.csv' inside ", out_dir_prefix)

#### Calibration / Brier / DCA / NRI / IDI (1-16 years) ####

out_dir <- file.path(out_dir_prefix, "model_compare_outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Helper: subject-specific risk at time t from a coxph model
predict_risk_at_time <- function(fit, newdata, t0) {
  # Use predictRisk directly to avoid survfit summary alignment ambiguity.
  risk <- as.numeric(riskRegression::predictRisk(fit, newdata = newdata, times = t0))
  pmin(pmax(risk, 1e-6), 1 - 1e-6)
}

make_calibration_table <- function(risk, data, t0, n_groups = 10) {
  tmp <- data.frame(
    time = data$time,
    status = data$status,
    risk = risk
  )
  tmp <- tmp[is.finite(tmp$risk) & !is.na(tmp$time) & !is.na(tmp$status), , drop = FALSE]
  tmp$group <- dplyr::ntile(tmp$risk, n_groups)

  out <- lapply(sort(unique(tmp$group)), function(g) {
    sub <- tmp[tmp$group == g, , drop = FALSE]
    sf <- survfit(Surv(time, status) ~ 1, data = sub)
    ss <- summary(sf, times = t0, extend = FALSE)
    surv_t <- if (length(ss$surv) == 0 || is.na(ss$surv[1])) NA_real_ else ss$surv[1]
    data.frame(
      group = g,
      pred = mean(sub$risk),
      obs = 1 - surv_t,
      n = nrow(sub)
    )
  })
  bind_rows(out)
}

# IPCW survival DCA at time t
survival_dca_ipcw <- function(data, time_horizon, risk_cols, model_labels = risk_cols) {
  d <- data
  d <- d[complete.cases(d[, c("time", "status", risk_cols), drop = FALSE]), , drop = FALSE]
  d$event_t <- as.integer(d$time <= time_horizon & d$status == 1)
  d$nonevent_t <- as.integer(d$time > time_horizon)
  # Censoring survival G(u) = P(C > u), estimated from censoring process.
  sf_c <- survfit(Surv(time, status == 0) ~ 1, data = d)
  get_G <- function(x) {
    s <- summary(sf_c, times = x, extend = TRUE)$surv
    ifelse(is.na(s), 1, pmax(s, 1e-6))
  }
  d$G_t <- get_G(time_horizon)
  d$G_time <- get_G(d$time)
  d$w_event <- d$event_t / d$G_time
  d$w_nonevent <- d$nonevent_t / d$G_t

  risk_values <- unlist(d[, risk_cols, drop = FALSE], use.names = FALSE)
  risk_values <- risk_values[is.finite(risk_values)]
  q <- suppressWarnings(quantile(risk_values, probs = c(0.05, 0.95), na.rm = TRUE, names = FALSE))
  thr_low <- ifelse(length(q) >= 1 && is.finite(q[1]), max(0.01, min(0.99, q[1])), 0.01)
  thr_high <- ifelse(length(q) >= 2 && is.finite(q[2]), max(0.01, min(0.99, q[2])), 0.99)
  if (thr_high - thr_low < 0.02) {
    thresholds <- seq(0.01, 0.99, by = 0.01)
  } else {
    thresholds <- unique(round(seq(thr_low, thr_high, by = 0.01), 2))
    thresholds <- thresholds[thresholds > 0 & thresholds < 1]
    if (length(thresholds) < 3) thresholds <- seq(0.01, 0.99, by = 0.01)
  }

  calc_nb <- function(risk, pt) {
    pred_pos <- as.integer(risk >= pt)
    mean(d$w_event * pred_pos, na.rm = TRUE) -
      mean(d$w_nonevent * pred_pos, na.rm = TRUE) * (pt / (1 - pt))
  }

  dca_rows <- list()
  for (j in seq_along(risk_cols)) {
    rk <- d[[risk_cols[j]]]
    dca_rows[[j]] <- data.frame(
      threshold = thresholds,
      net_benefit = sapply(thresholds, function(pt) calc_nb(rk, pt)),
      Model = model_labels[j]
    )
  }

  nb_all <- data.frame(
    threshold = thresholds,
    net_benefit = sapply(thresholds, function(pt) {
      mean(d$w_event, na.rm = TRUE) - mean(d$w_nonevent, na.rm = TRUE) * (pt / (1 - pt))
    }),
    Model = "Treat all"
  )
  nb_none <- data.frame(
    threshold = thresholds,
    net_benefit = 0,
    Model = "Treat none"
  )

  bind_rows(dca_rows, list(nb_all, nb_none))
}

years_to_run <- if (TARGET_OUTCOME == "liver_death") 4:16 else 1:16

# ----------- Build 5-fold OOF predicted risks for each time-point -----------
message(sprintf("--- Building 5-fold OOF Predicted Risks (%d-%d Years) ---", min(years_to_run), max(years_to_run)))
oof_long <- data.frame()

for (fold_idx in seq_along(folds_cv)) {
  message(sprintf("Fitting models and predicting risks for fold %d/%d...", fold_idx, length(folds_cv)))
  test_idx <- folds_cv[[fold_idx]]
  train_data <- analysis_df[-test_idx, , drop = FALSE]
    test_data <- analysis_df[test_idx, , drop = FALSE]
  
    formula_ncan_liverrisk <- as.formula(paste0("Surv(`", time_col, "`, `", status_col, "`) ~ `NCAN;Neurocan core protein` + `liverRisk score`"))
    formula_liverrisk <- as.formula(paste0("Surv(`", time_col, "`, `", status_col, "`) ~ `liverRisk score`"))

    fit_ncan_liverrisk <- coxph(
      formula_ncan_liverrisk,
      data = train_data,
      x = TRUE, y = TRUE
    )
    fit_liverrisk <- coxph(
      formula_liverrisk,
      data = train_data,
      x = TRUE, y = TRUE
    )

  for (tt in years_to_run) {
    risk_new <- predict_risk_at_time(fit_ncan_liverrisk, test_data, tt)
    risk_old <- predict_risk_at_time(fit_liverrisk, test_data, tt)

    oof_long <- bind_rows(
        oof_long,
        data.frame(
          id = test_idx,
          fold = fold_idx,
          year = tt,
          time = test_data[[time_col]],
          status = test_data[[status_col]],
          risk_new = risk_new,
          risk_old = risk_old
        )
      )
  }
}

oof_long <- oof_long %>% arrange(year, id)
oof_file <- file.path(out_dir, paste0("oof_risk_year_", min(years_to_run), "_to_", max(years_to_run), "y.csv"))
write.csv(oof_long, oof_file, row.names = FALSE)
message("Saved OOF predictions to: ", oof_file)

# ----------- Year-wise metrics based on OOF risks -----------
message("--- Starting Year-wise Metrics Calculation (Calibration only) ---")
calibration_all <- data.frame()
cal_plots <- list()
brier_all <- data.frame()
nri_idi_all <- data.frame()
dca_summary <- data.frame()
dca_obj_all <- list()
analysis_notes <- data.frame()

get_metric <- function(obj, candidates) {
  for (nm in candidates) {
    if (!is.null(obj[[nm]])) return(obj[[nm]])
  }
  return(NA_real_)
}

for (tt in years_to_run) {
  message(sprintf("\n==============================\nProcessing Year: %d\n==============================", tt))

  d_year <- oof_long %>% filter(year == tt)

  # Calibration on OOF risks
  message("  -> Calculating Calibration...")
  cal_new <- make_calibration_table(d_year$risk_new, d_year, tt, n_groups = 10) %>%
    mutate(Model = "NCAN_LiverRisk", year = tt)
  cal_old <- make_calibration_table(d_year$risk_old, d_year, tt, n_groups = 10) %>%
    mutate(Model = "LiverRisk_only", year = tt)
  cal_df <- bind_rows(cal_new, cal_old)
  calibration_all <- bind_rows(calibration_all, cal_df)
  write.csv(cal_df, file.path(out_dir, paste0("calibration_year_", tt, ".csv")), row.names = FALSE)

  p_cal <- ggplot(cal_df, aes(x = pred, y = obs, color = Model)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2) +
    scale_color_manual(values = c("NCAN_LiverRisk" = "#D73027", "LiverRisk_only" = "#4575B4")) +
    labs(title = paste0("Calibration Year ", tt, " (OOF)"), x = "Predicted risk", y = "Observed risk (KM)") +
    theme_bw(base_size = 12)
  print(p_cal)
  ggsave(file.path(out_dir, paste0("calibration_year_", tt, ".png")), p_cal, width = 9, height = 6, dpi = 300)
  cal_plots[[as.character(tt)]] <- p_cal
  # NRI / IDI (and DCA/Brier) disabled: only keep calibration plots.
  if (FALSE) {
    # Survival DCA on OOF risks (IPCW)
    message("  -> Calculating Survival DCA (IPCW)...")
  dca_fit <- survival_dca_ipcw(
    data = d_year,
    time_horizon = tt,
    risk_cols = c("risk_new", "risk_old"),
    model_labels = c("NCAN_LiverRisk", "LiverRisk_only")
  )
  thr_info <- range(dca_fit$threshold, na.rm = TRUE)
  thr_n <- length(unique(dca_fit$threshold))
  dca_obj_all[[paste0("year", tt)]] <- dca_fit

  p_dca <- ggplot(dca_fit, aes(x = threshold, y = net_benefit, color = Model)) +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = c(
      "NCAN_LiverRisk" = "#D73027",
      "LiverRisk_only" = "#4575B4",
      "Treat all" = "grey30",
      "Treat none" = "black"
    )) +
    labs(
      title = paste0("Survival DCA Year ", tt, " (OOF, IPCW)"),
      x = "Threshold probability",
      y = "Net benefit"
    ) +
    theme_bw(base_size = 12)
  print(p_dca)
  ggsave(file.path(out_dir, paste0("dca_year_", tt, ".png")), p_dca, width = 9, height = 6, dpi = 300)

  keep_idx <- !(d_year$time <= tt & d_year$status == 0)
  d_eval <- d_year[keep_idx, , drop = FALSE]
  d_eval$event_t <- ifelse(d_eval$time <= tt & d_eval$status == 1, 1, 0)

  # Brier on OOF risks (subjects censored before t excluded)
  message("  -> Calculating Brier Score...")
  brier_all <- bind_rows(
    brier_all,
    data.frame(year = tt, model = "NCAN_LiverRisk", Brier = mean((d_eval$event_t - d_eval$risk_new)^2, na.rm = TRUE)),
    data.frame(year = tt, model = "LiverRisk_only", Brier = mean((d_eval$event_t - d_eval$risk_old)^2, na.rm = TRUE))
  )

  dca_summary <- bind_rows(
    dca_summary,
    data.frame(
      year = tt,
      n_eval = nrow(d_eval),
      prevalence = mean(d_eval$event_t, na.rm = TRUE),
      threshold_low = ifelse(is.null(thr_info), NA_real_, thr_info[1]),
      threshold_high = ifelse(is.null(thr_info), NA_real_, thr_info[2]),
      threshold_n = ifelse(is.null(thr_n), NA_integer_, thr_n)
    )
  )

  # NRI / IDI with survival-aware estimator (survIDINRI)
  message("  -> Calculating NRI and IDI (survIDINRI)...")
  indata <- as.matrix(data.frame(
    time = d_year$time,
    status = d_year$status
  ))
  covs0 <- as.matrix(data.frame(
    risk_old = d_year$risk_old
  ))
  covs1 <- as.matrix(data.frame(
    risk_new = d_year$risk_new
  ))
  idi_obj <- IDI.INF(
    indata = indata,
    covs0 = covs0,
    covs1 = covs1,
    t0 = tt,
    npert = 200
  )
  
  idi_est <- idi_obj$m1[1]
  idi_lower <- idi_obj$m1[2]
  idi_upper <- idi_obj$m1[3]
  idi_pval <- idi_obj$m1[4]
  
  nri_est <- idi_obj$m2[1]
  nri_lower <- idi_obj$m2[2]
  nri_upper <- idi_obj$m2[3]
  nri_pval <- idi_obj$m2[4]

  nri_idi_all <- bind_rows(
    nri_idi_all,
    data.frame(
      year = tt,
      IDI = idi_est,
      IDI_lower = idi_lower,
      IDI_upper = idi_upper,
      IDI_pval = idi_pval,
      NRI = nri_est,
      NRI_lower = nri_lower,
      NRI_upper = nri_upper,
      NRI_pval = nri_pval
    )
  )

  analysis_notes <- bind_rows(
    analysis_notes,
    data.frame(
      year = tt,
      note = "Time-specific NRI/IDI estimated independently at each horizon; no multiplicity correction applied."
    ),
    data.frame(
      year = tt,
      note = "Brier (binary-at-t), DCA (IPCW survival utility), and IDI/NRI (survival influence-function) are complementary estimands."
    ),
    data.frame(
      year = tt,
      note = "Calibration uses Kaplan-Meier observed risk, while DCA uses IPCW weighting. They are complementary evaluation frameworks."
    ),
    data.frame(
      year = tt,
      note = "DCA censoring model uses marginal KM (Surv(time, status == 0) ~ 1) instead of G(t|X), acceptable in clinical contexts."
    )
  )
  }
}

message("\n--- Calibration Complete. Creating 4x4 Panel ---")

write.csv(
  calibration_all,
  file.path(out_dir, paste0("calibration_oof_year_", min(years_to_run), "_to_", max(years_to_run), "y.csv")),
  row.names = FALSE
)

empty_plot <- function(title = "") {
  ggplot() + theme_void() + labs(title = title)
}

grid_plots <- list()
for (tt in years_to_run) {
  key <- as.character(tt)
  if (!is.null(cal_plots[[key]])) grid_plots[[length(grid_plots) + 1]] <- cal_plots[[key]]
}
while (length(grid_plots) < 16) grid_plots[[length(grid_plots) + 1]] <- empty_plot("")

draw_4x4 <- function(plots) {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(4, 4)))
  for (i in seq_len(16)) {
    r <- ceiling(i / 4)
    c <- i - (r - 1) * 4
    print(plots[[i]], vp = grid::viewport(layout.pos.row = r, layout.pos.col = c))
  }
}

panel_file <- file.path(out_dir, paste0("calibration_panel_4x4_year_", min(years_to_run), "_to_", max(years_to_run), ".png"))
png(panel_file, width = 3600, height = 2400, res = 300)
draw_4x4(grid_plots)
dev.off()
message("Saved 4x4 calibration panel to: ", panel_file)

draw_4x4(grid_plots)

message(sprintf("\n[SUCCESS] Calibration-only run finished! Output files are located in: %s", out_dir))
message("The script completed at: ", Sys.time())
