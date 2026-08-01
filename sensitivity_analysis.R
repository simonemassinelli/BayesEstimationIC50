# Prepare single patient data
target_patient <- "ACH-000008"
df_sens <- data %>% filter(patient_id == target_patient)
params_sens <- fixed_params %>% filter(patient_id == target_patient)

# Function to run JAGS with a specific prior on IC50
fit_sensitivity <- function(shape_val, rate_val) {
  jags_data <- list(
    viability = df_sens$viability, dose = df_sens$dose, N = nrow(df_sens),
    h_fixed = params_sens$h_slope, min_fixed = params_sens$min_val, max_fixed = params_sens$max_val
  )
  
  jags_model <- sprintf("
  model {
    for (i in 1:N) {
      viability[i] ~ dnorm(mu[i], tau)
      mu[i] <- min_fixed + (max_fixed - min_fixed) / (1 + pow(dose[i] / IC50, h_fixed))
    }
    IC50 ~ dgamma(%f, %f)
    tau ~ dgamma(0.001, 0.001)
  }", shape_val, rate_val)
  
  fit <- jags(data = jags_data, parameters.to.save = c("IC50"),
              model.file = textConnection(jags_model),
              n.chains = 3, n.iter = 15000, n.burnin = 3000, quiet = TRUE)
  
  return(fit$BUGSoutput$summary["IC50", "mean"])
}

# Run the model with 3 different priors
mean_vague <- fit_sensitivity(0.001, 0.001) # Vague Prior
mean_weak  <- fit_sensitivity(1, 1)         # Weakly Informative
mean_strong<- fit_sensitivity(10, 10)       # Stronger Prior (centered at 1, but biological truth is ~0.06)

# Compile results
sensitivity_results <- data.frame(
  Prior_Type = c("Vague Gamma(0.001, 0.001)", "Weak Gamma(1, 1)", "Strong Gamma(10, 10)"),
  Posterior_Mean_IC50 = c(mean_vague, mean_weak, mean_strong)
)

cat("\n=== PRIOR SENSITIVITY ANALYSIS FOR", target_patient, "===\n")
print(sensitivity_results)