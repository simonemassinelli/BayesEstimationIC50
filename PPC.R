library(ggplot2)
library(patchwork)

# Extract the posterior mean of the IC50 from both JAGS models
ic50_generic_post <- fit_generic$BUGSoutput$mean$IC50
ic50_patient_post <- fit_patient$BUGSoutput$mean$IC50

# Create a dataframe mapping the posterior IC50 means to patient indices
posterior_ic50_df <- data.frame(
  patient_idx = 1:length(ic50_patient_post),
  ic50_est = ic50_patient_post
)

# Merge the posterior estimates and theoretical parameters back into the main dataset
data_ppc <- data %>%
  left_join(fixed_params, by = c("patient_id", "patient_idx")) %>%
  left_join(posterior_ic50_df, by = "patient_idx")

# Calculate expected viability for the Generic Pooled Model
data_ppc$pred_generic <- 0 + (1 - 0) / 
  (1 + (data_ppc$dose / ic50_generic_post)^data_ppc$generic_h_reference)

# Calculate expected viability for the Patient-Specific Model using theoretical limits
data_ppc$pred_specific <- data_ppc$min_val + (data_ppc$max_val - data_ppc$min_val) / 
  (1 + (data_ppc$dose / data_ppc$ic50_est)^data_ppc$h_slope)

# Build and display the diagnostic plot for the Generic Model
plot_generic <- ggplot(data_ppc, aes(x = pred_generic, y = viability)) +
  geom_point(alpha = 0.4, color = "purple") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1) +
  theme_minimal() +
  labs(
    title = "PPC: Generic Pooled Model",
    subtitle = paste("DIC:", round(fit_generic$BUGSoutput$DIC, 1)),
    x = "Predicted Viability",
    y = "Observed Viability"
  ) +
  coord_fixed(ratio = 1, xlim = c(0, 1.2), ylim = c(0, 1.2))

# Build and display the diagnostic plot for the Patient-Specific Model
plot_specific <- ggplot(data_ppc, aes(x = pred_specific, y = viability)) +
  geom_point(alpha = 0.4, color = "orchid") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1) +
  theme_minimal() +
  labs(
    title = "PPC: Patient-Specific Model",
    subtitle = paste("DIC:", round(fit_patient$BUGSoutput$DIC, 1)),
    x = "Predicted Viability",
    y = "Observed Viability"
  ) +
  coord_fixed(ratio = 1, xlim = c(0, 1.2), ylim = c(0, 1.2))

# Display both plots side-by-side
print(plot_generic + plot_specific)