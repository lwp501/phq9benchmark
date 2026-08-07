# Project: PHQ9 benchmarking
# Task 1: Estimate PHQ-9 benchmarking figures (using the Minami method) & do sense-check via metafor
# Task 2: Leave-one-out sensitivity analysis
# Task 3: Simulate critical values over sample sizes (N) and generate plot
# Author: Lewis Paton, University of York
# Updated: 7/8/26

library(readxl)
library(metafor)
library(ggplot2)
library(tidyr)


# --- Helper Functions ---
save_forest_tiff <- function(model, filename) {
  tiff(
    filename = filename,
    width = 10, height = 8, units = "in", res = 300, compression = "lzw"
  )
  par(mar = c(5, 4, 5, 2))
  forest(
    model,
    cex = 1.1,
    header = c("Study", "Effect Size [95% CI]"),
    top = 3
  )
  dev.off()
}


lambda_ncp <- function(samplesize, d.plus, dmin) {
  sqrt(samplesize) * (d.plus - dmin)
}

d_crit <- function(percentile, samplesize, dplus, dmin) {
  lambda <- lambda_ncp(samplesize, dplus, dmin)
  qt(percentile, df = samplesize - 1, ncp = lambda) / sqrt(samplesize)
}




## 1. Estimate benchmarks using Minami method

data <- read_excel("ITT_pure_dep_for_R_with_Mohr_update.xlsx")


# Estimate treatment SD in Mohr et al. (2012) using reported 95% CIs
data$treatment_pre_SD[4]  <- (sqrt(159) * (17.29 - 16.76)) / 1.96
data$treatment_post_SD[4] <- (sqrt(136) * (7.73 - 6.74))   / 1.96


# Estimate standardized effect sizes using Minami formula
data$treatment_d <- (1-(3/((4*data$treatment_pre_n-5))))*(
  (data$treatement_post_M-data$treatment_pre_M)/data$treatment_pre_SD)

data$control_d <- (1-(3/((4*data$control_pre_n-5))))*(
  (data$control_post_M-data$control_pre_M)/data$control_pre_SD)


# Calculate sampling variances
data$treatment_var <- (1 / data$treatment_pre_n) + ((data$treatment_d)^2 / (2 * data$treatment_pre_n))
data$control_var   <- (1 / data$control_pre_n)   + ((data$control_d)^2   / (2 * data$control_pre_n))

# Direct weighted benchmarks
benchmark_treatment <- sum(data$treatment_d / data$treatment_var) / sum(1 / data$treatment_var)
benchmark_control   <- sum(data$control_d   / data$control_var, na.rm = TRUE) / sum(1 / data$control_var, na.rm = TRUE)

message("Benchmark (Treatment): ", round(benchmark_treatment, 4))
message("Benchmark (Control):   ", round(benchmark_control, 4))


# Meta-Analysis Comparison using metafor

# Equal-effects models
model_treatment <- rma(yi = treatment_d, vi = treatment_var, data = data, slab = study, method = "EE")
model_control   <- rma(yi = control_d,   vi = control_var,   data = data, slab = study, method = "EE")

summary(model_treatment)
summary(model_control)

# Export forest plots
save_forest_tiff(model_treatment, "Fig_1_u.tiff")
save_forest_tiff(model_control,   "Fig_2_u.tiff")


## 2, Leave-one-out meta-analysis via metafor 
sensitivity_results <- leave1out(model_treatment)
print(sensitivity_results)


# 3. Simulation & Critical Value Plotting 


eff_bench_treatment_ITT <- -benchmark_treatment 
eff_bench_control_ITT   <- -benchmark_control 

ITT <- data.frame(
  n = c(100, 125, 150, 200, 250, 300, 400, 500, 600, 800, 900, 1000, 
        1250, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 8000, 10000)
)

ITT$crit_val_eff_control <- d_crit(0.95, ITT$n, eff_bench_control_ITT, -0.2)  
ITT$crit_val_eff_treat   <- d_crit(0.95, ITT$n, eff_bench_treatment_ITT, 0.2)

ITT_long <- pivot_longer(ITT, cols = starts_with("crit_val"), names_to = "variable", values_to = "value")

# Dynamic label formatting for benchmark reference lines
ys <- c(eff_bench_treatment_ITT, eff_bench_treatment_ITT - 0.2, eff_bench_control_ITT, eff_bench_control_ITT + 0.2)
ls <- c(
  sprintf("Efficacy Benchmark = %.2f", eff_bench_treatment_ITT),
  sprintf("Efficacy Benchmark - 0.2 = %.2f", eff_bench_treatment_ITT - 0.2),
  sprintf("TAU Benchmark = %.2f", eff_bench_control_ITT),
  sprintf("TAU Benchmark + 0.2 = %.2f", eff_bench_control_ITT + 0.2)
)

cuts <- data.frame(
  yintercept = ys,
  Linetype   = factor(ls, levels = ls),
  variable   = c("crit_val_eff_treat", "crit_val_eff_treat", "crit_val_eff_control", "crit_val_eff_control")
)

# Plot generation
p_sim <- ggplot(ITT_long, aes(x = n, y = value, color = variable)) +
  geom_line(linewidth = 1.1) +
  geom_hline(data = cuts, aes(yintercept = yintercept, linetype = Linetype, color = variable), linewidth = 0.8) + 
  scale_x_continuous(
    trans = "log10",
    breaks = c(100, 125, 150, 200, 250, 300, 400, 500, 600, 800, 1000, 
               1250, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 8000, 10000),
    limits = c(100, 10000)
  ) + 
  scale_y_continuous(
    breaks = seq(0.7, 1.5, by = 0.1),
    limits = c(0.7, 1.5)
  ) + 
  scale_color_manual(
    name = "Condition",
    values = c("crit_val_eff_control" = "#F8766D", "crit_val_eff_treat" = "#619CFF"),
    labels = c("TAU Critical Value", "Efficacy Critical Value")
  ) +
  scale_linetype_manual(
    name = "Benchmark Thresholds",
    values = c("dashed", "dotdash", "twodash", "dotted")
  ) +
  labs(x = "Sample Size N", y = "Effect Size d") +
  theme_classic(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# Save figure as high-res TIFF
ggsave(
  filename = "crit_vals_update_R1.tiff",
  plot = p_sim,
  width = 9.20,
  height = 7.85,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

