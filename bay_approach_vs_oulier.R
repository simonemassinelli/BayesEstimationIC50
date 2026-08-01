library(ggplot2)
library(dplyr)

# Define the target patient
target_patient <- "ACH-000008"

# Extract raw data for the target patient
df_pat <- data %>% filter(patient_id == target_patient)
p_idx <- df_pat$patient_idx[1]

# Extract deterministic parameters for this specific patient
params_pat <- fixed_params %>% filter(patient_id == target_patient)
min_p <- params_pat$min_val
max_p <- params_pat$max_val
h_p <- params_pat$h_slope

# Extract the Bayesian estimated IC50 using the posterior mean
ic50_est <- fit_patient$BUGSoutput$mean$IC50[p_idx]

# Generate a sequence of doses to draw the continuous theoretical curve
dose_seq <- exp(seq(log(min(df_pat$dose)), log(max(df_pat$dose)), length.out = 100))

# Calculate predicted viability for the dose sequence using the LL.4 equation
pred_viab <- min_p + (max_p - min_p) / (1 + (dose_seq / ic50_est)^h_p)
df_curve <- data.frame(dose = dose_seq, viability = pred_viab)

# Build the dose-response plot
plot_curve <- ggplot() +
  # Raw experimental data points
  geom_point(data = df_pat, aes(x = dose, y = viability), 
             color = "black", size = 3, alpha = 0.7) +
  
  # Smooth theoretical curve (Academic Blue)
  geom_line(data = df_curve, aes(x = dose, y = viability), 
            color = "blue", linewidth = 1.2) + 
  
  # Green dashed line for Bayesian IC50
  geom_vline(xintercept = ic50_est, linetype = "dashed", 
             color = "green", linewidth = 1.2) + 
  
  # Red dashed line for the outlier at dose = 1.0
  geom_vline(xintercept = 1.0, linetype = "dashed", 
             color = "red", linewidth = 1.2) +
  
  # Red circle highlighting the outlier data point
  geom_point(data = df_pat %>% filter(dose == 1), aes(x = dose, y = viability), 
             color = "red", size = 5, shape = 1, stroke = 1.5) +
  
  # Text labels for the dose values at the bottom (-Inf anchors it to the bottom axis)
  annotate("text", x = ic50_est, y = 0.31, label = sprintf("%.4f", ic50_est), 
           color = "green", vjust = -1, hjust = -0.1, fontface = "bold", size = 4.5) +
  annotate("text", x = 1.0, y = 0.31, label = "1.000", 
           color = "red", vjust = -1, hjust = 1.2, fontface = "bold", size = 4.5) +
  
  scale_x_log10() +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", color = "black"),
    plot.subtitle = element_text(color = "grey30"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = paste("Robust Bayesian Fit vs Outlier:", target_patient),
    subtitle = "Green: Bayesian IC50 | Red: Extreme isolated outlier ignored by the hierarchical model",
    x = "Dose [Log Scale]",
    y = "Cell Viability"
  )

# Display the plot
print(plot_curve)