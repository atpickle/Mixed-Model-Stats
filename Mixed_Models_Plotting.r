# Packages (guard installs)
pkgs <- c("lmerTest","lme4","ggplot2","lattice","dplyr")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
lapply(pkgs, library, character.only = TRUE)

# Working directory
setwd("/Users/allisonpickle/Desktop")

# Data
input_file  <- "PD_Microglia Area_RMixedFormated.csv"
process_data <- read.csv(input_file)

# Factors (hard-coded order, reused for plots)
baseline_treatment <- "Ctrl"

process_data$Treatment <- factor(
  process_data$Treatment,
  levels = c(
    "Ctrl",
    "PD",
    "1ugLPS",
    "5ugLPS",
    "10ugLPS",
    "1ugLPSPD",
    "5ugLPSPD",
    "10ugLPSPD"
  )
)
process_data$Well <- factor(process_data$Well)

# Model 1 (raw)
RawFit <- lmerTest::lmer(Results ~ Treatment + (1|Treatment:Well), data = process_data)
anova(RawFit); summary(RawFit)

# macOS graphics device (optional) - COMMENTED OUT to prevent pop-up
# if (Sys.info()[["sysname"]] == "Darwin" && capabilities("aqua")) {
#   if (is.null(dev.list()) || !any(grepl("quartz", names(dev.list())))) quartz()
# }

# Helper plotting function: ONLY std residuals vs fitted + QQ
plot_diag <- function(fit, tag) {
  par(mfrow = c(1, 2))
  fitted_vals <- fitted(fit)
  std_resid   <- resid(fit, type = "pearson")
  plot(fitted_vals, std_resid,
       xlab = "Fitted Values",
       ylab = "Standardized (Pearson) Residuals",
       main = paste(tag, "Std Residuals vs Fitted"))
  abline(h = 0, col = "red", lwd = 2)
  print(lattice::qqmath(fit, id = 0.05,
                        main = paste(tag, "QQ Plot")))
  par(mfrow = c(1, 1))
}

# Diagnostics for raw model
# plot_diag(RawFit, "Raw")  # COMMENTED OUT to prevent pop-up

# Transform: log
process_data$logResults <- log(process_data$Results + 0.5)
LogFit <- lmerTest::lmer(logResults ~ Treatment + (1|Treatment:Well), data = process_data)
anova(LogFit); summary(LogFit)
# plot_diag(LogFit, "Log")  # COMMENTED OUT to prevent pop-up

# ---- Save ANOVA and summary outputs for RawFit and LogFit ----
save_model_summaries <- function() {
  csv_base <- tools::file_path_sans_ext(basename(input_file))
  out_dir  <- paste0("MM_Outputs_", csv_base)
  dir.create(out_dir, showWarnings = FALSE)

  out_file <- file.path(out_dir, "Model_Summaries.txt")

  sink(out_file)
  cat("=== RawFit (Results ~ Treatment + (1|Treatment:Well)) ===\n\n")
  cat("ANOVA:\n"); print(anova(RawFit))
  cat("\nSummary:\n"); print(summary(RawFit))

  cat("\n\n=== LogFit (logResults ~ Treatment + (1|Treatment:Well)) ===\n\n")
  cat("ANOVA:\n"); print(anova(LogFit))
  cat("\nSummary:\n"); print(summary(LogFit))
  sink()

  message("Saved ANOVA and summary outputs to: ", out_file)
}

# ---- Save diagnostics: ONLY std_resid vs fitted + QQ, plus residual data ----
save_diagnostics <- function(save_data = TRUE) {
  csv_base <- tools::file_path_sans_ext(basename(input_file))
  out_dir  <- paste0("MM_Outputs_", csv_base)
  dir.create(out_dir, showWarnings = FALSE)
  
  models <- list(
    RawFit = RawFit,
    LogFit = LogFit
  )
  
  all_resid <- list()
  
  for (nm in names(models)) {
    fit <- models[[nm]]
    fitted_vals <- fitted(fit)
    std_resid   <- resid(fit, type = "pearson")
    
    ## 1) STANDARDIZED (Pearson) residuals vs fitted
    png(file.path(out_dir, paste0(nm, "_residuals.png")),
        width = 1200, height = 1000, res = 150)
    plot(fitted_vals, std_resid,
         xlab = "Fitted Values",
         ylab = "Standardized (Pearson) Residuals",
         main = paste(nm, "Std Residuals vs Fitted"))
    abline(h = 0, col = "red", lwd = 2)
    dev.off()
    
    ## 2) QQ plot
    png(file.path(out_dir, paste0(nm, "_qq.png")),
        width = 1200, height = 1000, res = 150)
    print(lattice::qqmath(fit, id = 0.05,
                          main = paste(nm, "QQ Plot")))
    dev.off()
    
    if (save_data) {
      all_resid[[nm]] <- data.frame(
        Model       = nm,
        Fitted      = fitted_vals,
        StdResidual = std_resid,
        RawResidual = resid(fit)
      )
    }
  }
  
  if (save_data) {
    write.csv(
      do.call(rbind, all_resid),
      file.path(out_dir, "All_Calculated_Residuals.csv"),
      row.names = FALSE
    )
  }
  
  message("Saved standardized (Pearson) residual vs fitted and QQ plots for each model to: ", out_dir)
}

# ---------- Helpers: extract p-values and make barplots with stars ----------

extract_treatment_pvals <- function(fit, baseline_name = "Ctrl") {
  fe <- coef(summary(fit))

  rn <- rownames(fe)
  is_trt <- grepl("^Treatment", rn)
  if (!any(is_trt)) {
    stop("No Treatment* rows found in fixed effects.")
  }

  trt_rows  <- fe[is_trt, , drop = FALSE]
  trt_names <- sub("^Treatment", "", rownames(trt_rows))

  p_col <- grep("Pr\\(", colnames(trt_rows), value = TRUE)[1]
  if (is.na(p_col)) stop("Could not find p-value column in fixed effects.")

  pvals <- trt_rows[, p_col]
  names(pvals) <- trt_names

  data.frame(
    Treatment = c(baseline_name, trt_names),
    p.value   = c(NA, as.numeric(pvals)),
    stringsAsFactors = FALSE
  )
}

add_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  ""
}

build_bar_summary_with_sig <- function(data, fit, response_col, baseline_name = "Ctrl") {
  summary_df <- data |>
    dplyr::group_by(Treatment) |>
    dplyr::summarise(
      mean = mean(.data[[response_col]], na.rm = TRUE),
      se   = sd(.data[[response_col]], na.rm = TRUE) /
             sqrt(sum(!is.na(.data[[response_col]]))),
      .groups = "drop"
    )

  # Re-apply original factor order to keep plotting order
  summary_df$Treatment <- factor(
    summary_df$Treatment,
    levels = levels(data$Treatment)
  )

  p_df <- extract_treatment_pvals(fit, baseline_name)
  p_df$Treatment <- factor(p_df$Treatment, levels = levels(data$Treatment))

  out <- dplyr::left_join(summary_df, p_df, by = "Treatment")
  out$stars <- vapply(out$p.value, add_stars, character(1))
  out
}

plot_bar_p_values <- function(summ_df, y_label, title_text) {
  # Calculate the maximum height needed (mean + se + some padding for stars)
  max_y <- max(summ_df$mean + summ_df$se, na.rm = TRUE)
  
  ggplot(summ_df, aes(x = Treatment, y = mean)) +
    geom_col(fill = "gray70") +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    geom_text(
      aes(y = mean + se, label = stars), # Position at top of error bar
      vjust = -0.5,                      # Move up slightly from that point
      size = 5
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + # Increased to 20% extra space
    labs(
      x = "Treatment",
      y = y_label,
      title = title_text
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5) # Center the title
    )
}

# ---------- Build and save two barplots: RawFit and LogFit ----------

csv_base <- tools::file_path_sans_ext(basename(input_file))
out_dir  <- paste0("MM_Outputs_", csv_base)
dir.create(out_dir, showWarnings = FALSE)

# 1) Raw scale barplot, stars from RawFit fixed effects
raw_bar_summary <- build_bar_summary_with_sig(
  data         = process_data,
  fit          = RawFit,
  response_col = "Results",
  baseline_name = baseline_treatment
)

raw_bar_plot <- plot_bar_p_values(
  raw_bar_summary,
  y_label   = "Results (mean ± SE)",
  title_text = "Treatment (RawFit p-values)"
)

ggplot2::ggsave(
  filename = file.path(out_dir, "Barplot_RawFit.png"),
  plot     = raw_bar_plot,
  width    = 6,
  height   = 4,
  dpi      = 300
)

# 2) Log scale barplot, stars from LogFit fixed effects
log_bar_summary <- build_bar_summary_with_sig(
  data         = process_data,
  fit          = LogFit,
  response_col = "logResults",
  baseline_name = baseline_treatment
)

log_bar_plot <- plot_bar_p_values(
  log_bar_summary,
  y_label   = "log(Results + 0.5) (mean ± SE)",
  title_text = "Treatment (LogFit p-values)"
)

ggplot2::ggsave(
  filename = file.path(out_dir, "Barplot_LogFit.png"),
  plot     = log_bar_plot,
  width    = 6,
  height   = 4,
  dpi      = 300
)

# ---- Run saves at the end ----
save_diagnostics()
save_model_summaries()