# Clear the global environment
rm(list=ls())

# Load required libraries
library(dplyr)
library(R2jags)

# Load raw data and the newly generated theoretical parameters
raw_data <- TRAMETINIB_jags_data
params_data <- trametinib_patient_parameters_theoretical

# Filter raw data to strictly include patients with valid estimated parameters
data <- raw_data %>%
  filter(patient_id %in% params_data$patient_id)

# Create a sequential numeric index for JAGS
data$patient_idx <- as.numeric(as.factor(data$patient_id))

# Generate a mapping dataframe to link indices to patient IDs
patient_mapping <- data %>%
  select(patient_id, patient_idx) %>%
  distinct() %>%
  arrange(patient_idx)

colnames(params_data)

# Merge fixed theoretical parameters with the indices
fixed_params <- patient_mapping %>%
  inner_join(params_data, by = "patient_id")

# --- Generic Pooled Model ---
jags_data_generic <- list(
  viability = data$viability,
  dose = data$dose,
  N = nrow(data),
  h_fixed = fixed_params$generic_h_reference[1], 
  min_fixed = 0,                                 
  max_fixed = 1                                  
)

jags_model_generic <- "
model {
  for (i in 1:N) {
    viability[i] ~ dnorm(mu[i], tau)
    mu[i] <- min_fixed + (max_fixed - min_fixed) / (1 + pow(dose[i] / IC50, h_fixed))
  }
  
  # Non-informative Gamma prior for the global IC50
  IC50 ~ dgamma(0.001, 0.001)
  tau ~ dgamma(0.001, 0.001)
}
"

fit_generic <- jags(
  data = jags_data_generic,
  parameters.to.save = c("IC50", "tau"),
  model.file = textConnection(jags_model_generic),
  n.chains = 3,
  n.iter = 20000, 
  n.burnin = 5000,
  DIC = TRUE
)

# --- Patient-Specific Hierarchical Model ---
jags_data_patient <- list(
  viability = data$viability,
  dose = data$dose,
  patient = data$patient_idx,
  N = nrow(data),
  n_patients = nrow(fixed_params),
  h_fixed = fixed_params$h_slope,       
  min_fixed = fixed_params$min_val,     
  max_fixed = fixed_params$max_val      
)

jags_model_patient <- "
model {
  for (i in 1:N) {
    viability[i] ~ dnorm(mu[i], tau)
    mu[i] <- min_fixed[patient[i]] + 
             (max_fixed[patient[i]] - min_fixed[patient[i]]) / 
             (1 + pow(dose[i] / IC50[patient[i]], h_fixed[patient[i]]))
  }
  
  # Estimate IC50 for each individual patient using a hierarchical Gamma prior
  for (p in 1:n_patients) {
    IC50[p] ~ dgamma(shape_ic50, rate_ic50)
  }
  
  shape_ic50 ~ dgamma(0.01, 0.01)
  rate_ic50 ~ dgamma(0.01, 0.01)
  tau ~ dgamma(0.001, 0.001)
}
"

fit_patient <- jags(
  data = jags_data_patient,
  parameters.to.save = c("IC50", "shape_ic50", "rate_ic50", "tau"),
  model.file = textConnection(jags_model_patient),
  n.chains = 3,
  n.iter = 20000,
  n.burnin = 5000,
  DIC = TRUE
)

# Print Model Comparison
cat("\n=== DIC COMPARISON ===\n")
cat("Generic Pooled Model DIC:", fit_generic$BUGSoutput$DIC, "\n")
cat("Patient-Specific Model DIC:", fit_patient$BUGSoutput$DIC, "\n")