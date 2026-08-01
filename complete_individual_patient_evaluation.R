library(coda)
library(ggplot2)
library(dplyr)
library(R2jags)

# ==============================================================================
# SHOWCASE PATIENT SETUP
# ==============================================================================
target_patient <- "ACH-000327"
p_idx <- fixed_params$patient_idx[fixed_params$patient_id == target_patient]
param_name <- paste0("IC50[", p_idx, "]")

cat("\n========================================================\n")
cat("IN-DEPTH ANALYSIS FOR PATIENT:", target_patient, "\n")
cat("========================================================\n")

# ==============================================================================
# 1. MCMC DIAGNOSTICS (Traceplot & Autocorrelation)
# ==============================================================================
# Extract MCMC chains
mcmc_samples <- as.mcmc(fit_patient)
combined_chains <- as.mcmc(do.call(rbind, mcmc_samples)) # Combine chains for autocorrelation

par(mfrow = c(1, 2))

# Traceplot (separated chains to visually check mixing)
traceplot(mcmc_samples[, param_name], 
          main = paste("Traceplot:", target_patient), 
          col = c("darkred", "darkblue", "yellow"))

# Autocorrelation (combined chains to avoid plot overlapping issues)
# 1. Estraiamo i valori di autocorrelazione
acf_vals <- autocorr(combined_chains[, param_name])
lags <- 0:(length(acf_vals) - 1)

# 2. Ricreiamo ESATTAMENTE lo stesso grafico di autocorr.plot (con type = "h")
plot(x = lags, y = acf_vals, 
     type = "h",          # Disegna le stesse identiche linee verticali
     main = paste("Autocorrelation:", target_patient),
     xlab = "Lag", 
     ylab = "Autocorrelation",
     ylim = c(-0.1, 1)) # Taglio dell'asse Y

# 3. Aggiungiamo i pallini e la linea tratteggiata dello zero
points(x = lags, y = acf_vals, pch = 16, cex = 0.7)
abline(h = 0, lty = 2)        

# Reset plotting layout
par(mfrow = c(1, 1))

# Extract and print formal diagnostic statistics
summary_stats <- fit_patient$BUGSoutput$summary[param_name, ]
cat("\n--- FORMAL DIAGNOSTICS ---\n")
cat("Gelman-Rubin (Rhat):", round(summary_stats["Rhat"], 4), "(Target: ~1.0)\n")
cat("Effective Sample Size (n.eff):", summary_stats["n.eff"], "\n")

# ==============================================================================
# 2. PRIOR VS POSTERIOR PLOT (FIXED LEGEND)
# ==============================================================================
# Extract posterior samples for the specific patient
post_ic50 <- as.numeric(fit_patient$BUGSoutput$sims.list$IC50[, p_idx])
df_post <- data.frame(IC50 = post_ic50)

# Calculate the maximum density height and optimal X limit
max_y <- max(density(post_ic50)$y)
x_max <- quantile(post_ic50, 0.99) * 1.5

# Create the plot with a forced unified legend
plot_prior_post <- ggplot(df_post, aes(x = IC50)) +
  # Posterior density
  geom_density(aes(color = "Posterior", fill = "Posterior"), alpha = 0.6, linewidth = 1.2) +
  
  # Prior density (added dummy fill mapping to force legend unification)
  stat_function(aes(color = "Prior (Vague)", fill = "Prior (Vague)"), 
                fun = dgamma, args = list(shape = 0.001, rate = 0.001), 
                linetype = "dashed", linewidth = 1.2) +
  
  # Unified legend mapping
  scale_color_manual(name = "Distribution", values = c("Posterior" = "purple", "Prior (Vague)" = "grey30")) +
  scale_fill_manual(name = "Distribution", values = c("Posterior" = "orchid", "Prior (Vague)" = "transparent")) +
  
  # Force X-axis to start strictly at 0 and extend to x_max, removing left padding
  scale_x_continuous(limits = c(0, x_max), expand = expansion(mult = c(0, 0.05))) +
  
  # Keep Y-axis zoom to frame the peak perfectly
  coord_cartesian(ylim = c(0, max_y * 1.05)) +
  
  # Elegant classic theme
  theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  ) +
  labs(
    title = paste("Prior vs Posterior Learning:", target_patient),
    x = "IC50 Estimate",
    y = "Density"
  )

print(plot_prior_post)

# ==============================================================================
# 3. PRIOR SENSITIVITY ANALYSIS
# ==============================================================================
# Subset raw data and parameters for the target patient
df_sens <- data %>% filter(patient_id == target_patient)
params_sens <- fixed_params %>% filter(patient_id == target_patient)

# Function to run JAGS rapidly with a specific Gamma prior on IC50
fit_sensitivity <- function(shape_val, rate_val) {
  jags_data <- list(
    viability = df_sens$viability, 
    dose = df_sens$dose, 
    N = nrow(df_sens),
    h_fixed = params_sens$h_slope, 
    min_fixed = params_sens$min_val, 
    max_fixed = params_sens$max_val
  )
  
  # Dynamic JAGS model string embedding the chosen prior parameters
  jags_model <- sprintf("
  model {
    for (i in 1:N) {
      viability[i] ~ dnorm(mu[i], tau)
      mu[i] <- min_fixed + (max_fixed - min_fixed) / (1 + pow(dose[i] / IC50, h_fixed))
    }
    IC50 ~ dgamma(%f, %f)
    tau ~ dgamma(0.001, 0.001)
  }", shape_val, rate_val)
  
  # Fit the temporary model
  fit <- jags(data = jags_data, parameters.to.save = c("IC50"),
              model.file = textConnection(jags_model),
              n.chains = 3, n.iter = 15000, n.burnin = 3000, quiet = TRUE)
  
  return(fit$BUGSoutput$summary["IC50", "mean"])
}

# Run the model under 3 structurally different priors
mean_vague  <- fit_sensitivity(0.001, 0.001) # Baseline vague prior
mean_weak   <- fit_sensitivity(1, 1)         # Weakly informative prior
mean_strong <- fit_sensitivity(10, 10)       # Strong prior pulling towards 1

# Compile and print the results
sensitivity_results <- data.frame(
  Prior_Type = c("Vague Gamma(0.001, 0.001)", "Weak Gamma(1, 1)", "Strong Gamma(10, 10)"),
  Posterior_Mean_IC50 = c(mean_vague, mean_weak, mean_strong)
)

cat("\n--- PRIOR SENSITIVITY RESULTS ---\n")
print(sensitivity_results) 