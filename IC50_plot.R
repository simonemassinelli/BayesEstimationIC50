library(ggplot2)
library(dplyr)
library(patchwork)

# ==============================================================================
# 1. SETUP & EMPIRICAL CALCULATIONS
# ==============================================================================

# Find the global closest point assuming fixed generic limits of 0 and 1
global_closest <- data %>%
  mutate(dist_to_mid = abs(viability - 0.5)) %>%
  arrange(dist_to_mid) %>%
  slice(1)

emp_dose_gen <- global_closest$dose
emp_viab_gen <- global_closest$viability

post_ic50_gen <- data.frame(IC50 = as.numeric(fit_generic$BUGSoutput$sims.list$IC50))

# ==============================================================================
# 2. THEMES & COLORS
# ==============================================================================

# Human-readable, academic colors
color_global_fill <- "purple" # Muted Red for Global
color_pat_fill    <- "orchid" # Steel Blue for Specific Patients
color_line        <- "red" # Dark Slate for the reference line

# Clean academic theme without the bulky legend
academic_theme <- theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 14, color = "black"),
    plot.subtitle = element_text(size = 11, color = "grey30", margin = margin(b = 10)),
    panel.grid.minor = element_blank(),
    legend.position = "none" # Legend removed to maximize plot width
  )

# ==============================================================================
# 3. GLOBAL MODEL PLOT
# ==============================================================================

# Adaptive X-limit to prevent squishing
x_limit_global <- max(quantile(post_ic50_gen$IC50, 0.99), emp_dose_gen * 1.2)

plot_global <- ggplot(post_ic50_gen, aes(x = IC50)) +
  geom_density(fill = color_global_fill, alpha = 0.7, color = "black", linewidth = 0.5, adjust = 1.5) +
  # Note: xintercept is OUTSIDE aes() to fix the missing line bug
  geom_vline(xintercept = emp_dose_gen, linetype = "dashed", color = color_line, linewidth = 1.2) +
  coord_cartesian(xlim = c(0, x_limit_global)) +
  academic_theme +
  labs(
    title = "Global Population IC50",
    subtitle = sprintf("Density: Posterior | Dashed Line: Empirical Dose = %.4f (Obs. Viability: %.2f)", 
                       emp_dose_gen, emp_viab_gen),
    x = "IC50 Estimate",
    y = "Density"
  )

# ==============================================================================
# 4. PATIENT-SPECIFIC PLOTS
# ==============================================================================

# Calculate theoretical specific targets based on limits
patients_midpoint <- data %>%
  left_join(fixed_params, by = c("patient_id", "patient_idx")) %>%
  mutate(
    target_viab = (min_val + max_val) / 2,
    dist_to_target = abs(viability - target_viab)
  ) %>%
  group_by(patient_idx, patient_id) %>%
  arrange(dist_to_target) %>%
  slice(1) %>%
  ungroup()

# Select the top 2 best candidates
best_candidates <- patients_midpoint %>%
  filter(dist_to_target < 0.05) %>%
  head(2)

patient_plots <- list()

for (i in 1:nrow(best_candidates)) {
  p_id <- best_candidates$patient_id[i]
  p_idx <- best_candidates$patient_idx[i]
  emp_dose <- best_candidates$dose[i]
  emp_viab <- best_candidates$viability[i]
  
  post_samples <- as.numeric(fit_patient$BUGSoutput$sims.list$IC50[, p_idx])
  df_post <- data.frame(IC50 = post_samples)
  
  # Adaptive X-axis tailored to each patient
  x_limit_pat <- max(quantile(df_post$IC50, 0.99), emp_dose * 1.2)
  
  p_plot <- ggplot(df_post, aes(x = IC50)) +
    geom_density(fill = color_pat_fill, alpha = 0.7, color = "black", linewidth = 0.5, adjust = 1.5) +
    # Note: xintercept is OUTSIDE aes() to fix the missing line bug in loops
    geom_vline(xintercept = emp_dose, linetype = "dashed", color = color_line, linewidth = 1.2) +
    coord_cartesian(xlim = c(0, x_limit_pat)) +
    academic_theme +
    labs(
      title = paste("Patient-Specific IC50:", p_id),
      subtitle = sprintf("Density: Posterior | Dashed Line: Empirical Dose = %.4f (Obs. Viability: %.2f)", 
                         emp_dose, emp_viab),
      x = "IC50 Estimate",
      y = "Density"
    )
  
  patient_plots[[i]] <- p_plot
}

# ==============================================================================
# 5. COMBINE PLOTS WITH PATCHWORK
# ==============================================================================

# Stack vertically
final_plot <- plot_global / patient_plots[[1]] / patient_plots[[2]] +
  plot_annotation(
    title = "Bayesian Posterior Densities vs. Empirical Observations",
    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 15)))
  )

print(final_plot)