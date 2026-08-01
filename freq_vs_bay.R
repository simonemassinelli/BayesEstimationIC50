library(ggplot2)
library(dplyr)
library(patchwork)
library(drc) 

# ==============================================================================
# 1. THEMES & COLORS
# ==============================================================================

# Human-readable, academic colors
color_global_fill <- "purple" # Muted Red for Global
color_pat_fill    <- "orchid" # Steel Blue for Specific Patients
color_line        <- "red" # Dark Slate for the reference line

# Clean academic theme
academic_theme <- theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 14, color = "black"),
    plot.subtitle = element_text(size = 11, color = "grey30", margin = margin(b = 10)),
    panel.grid.minor = element_blank(),
    legend.position = "none" # Legend removed to maximize plot width
  )

# ==============================================================================
# 2. GLOBAL MODEL COMPARISON
# ==============================================================================

# Calculate the global frequentist IC50 to use as our reference line
fit_gen_plot <- tryCatch(
  drm(viability ~ dose, data = data, fct = LL.4(names = c("h_slope", "lower_limit", "upper_limit", "IC50"))),
  error = function(e) NULL
)
freq_ic50_global <- as.numeric(coef(fit_gen_plot)["IC50:(Intercept)"])

# Extract Bayesian posterior samples for the global IC50
post_ic50_gen <- data.frame(IC50 = as.numeric(fit_generic$BUGSoutput$sims.list$IC50))

# Adaptive X-limit to prevent squishing
x_limit_global <- max(quantile(post_ic50_gen$IC50, 0.99), freq_ic50_global * 1.2)

# Plot Global Model: Bayesian Density vs Frequentist Point Estimate
plot_global <- ggplot(post_ic50_gen, aes(x = IC50)) +
  geom_density(fill = color_global_fill, alpha = 0.7, color = "black", linewidth = 0.5, adjust = 1.5) +
  geom_vline(xintercept = freq_ic50_global, linetype = "dashed", color = color_line, linewidth = 1.2) +
  coord_cartesian(xlim = c(0, x_limit_global)) +
  academic_theme +
  labs(
    title = "Global Population IC50",
    subtitle = sprintf("Density: Posterior | Dashed Line: Frequentist Estimate = %.4f", freq_ic50_global),
    x = "IC50 Estimate",
    y = "Density"
  )

# ==============================================================================
# 3. PATIENT-SPECIFIC MODEL COMPARISON
# ==============================================================================

# Select a couple of interesting patients to visualize
target_patients <- c("ACH-000007", "ACH-000008")

# Filter the theoretical parameters for these specific patients
candidates <- fixed_params %>%
  filter(patient_id %in% target_patients)

patient_plots <- list()

for (i in 1:nrow(candidates)) {
  p_id <- candidates$patient_id[i]
  p_idx <- candidates$patient_idx[i]
  
  # Retrieve the theoretical Frequentist IC50 from our CSV data
  freq_ic50 <- candidates$ic50_frequentist[i]
  
  # Extract Bayesian posterior samples for this specific patient
  post_samples <- as.numeric(fit_patient$BUGSoutput$sims.list$IC50[, p_idx])
  df_post <- data.frame(IC50 = post_samples)
  
  # Calculate a dynamic upper limit for the X-axis to cut off the long flat tail
  x_limit_pat <- max(quantile(df_post$IC50, 0.99), freq_ic50 * 1.2)
  
  # Plot Specific Model: Bayesian Density vs Frequentist Point Estimate
  p_plot <- ggplot(df_post, aes(x = IC50)) +
    geom_density(fill = color_pat_fill, alpha = 0.7, color = "black", linewidth = 0.5, adjust = 1.5) +
    geom_vline(xintercept = freq_ic50, linetype = "dashed", color = color_line, linewidth = 1.2) +
    coord_cartesian(xlim = c(0, x_limit_pat)) +
    academic_theme +
    labs(
      title = paste("Patient-Specific IC50:", p_id),
      subtitle = sprintf("Density: Posterior | Dashed Line: Frequentist Estimate = %.4f", freq_ic50),
      x = "IC50 Estimate",
      y = "Density"
    )
  
  patient_plots[[i]] <- p_plot
}

# ==============================================================================
# 4. COMBINE AND DISPLAY PLOTS
# ==============================================================================

# Stack the plots vertically (1 column, 3 rows) to maximize x-axis width
final_plot <- plot_global / patient_plots[[1]] / patient_plots[[2]] +
  plot_annotation(
    title = "Bayesian Posterior Densities vs. Frequentist Estimates",
    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 15)))
  )

print(final_plot)