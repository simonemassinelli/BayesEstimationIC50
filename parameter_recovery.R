library(dplyr)
library(R2jags)

# 1. Define true parameters for our simulated patient
true_ic50 <- 0.08
true_h <- 1.5
true_min <- 0.2
true_max <- 1.0
true_tau <- 100 # Precision (implies a standard deviation of 0.1 for the noise)

# 2. Simulate experimental doses (log-spaced)
sim_doses <- exp(seq(log(0.001), log(1), length.out = 15))

# 3. Simulate viability using the LL.4 equation + Normal noise
set.seed(123)
mu_sim <- true_min + (true_max - true_min) / (1 + (sim_doses / true_ic50)^true_h)
sim_viability <- rnorm(length(sim_doses), mean = mu_sim, sd = sqrt(1/true_tau))

# 4. Prepare data for JAGS
jags_data_sim <- list(
  viability = sim_viability,
  dose = sim_doses,
  N = length(sim_doses),
  h_fixed = true_h,
  min_fixed = true_min,
  max_fixed = true_max
)

# 5. JAGS Model for a single patient
jags_model_sim <- "
model {
  for (i in 1:N) {
    viability[i] ~ dnorm(mu[i], tau)
    mu[i] <- min_fixed + (max_fixed - min_fixed) / (1 + pow(dose[i] / IC50, h_fixed))
  }
  IC50 ~ dgamma(0.001, 0.001)
  tau ~ dgamma(0.001, 0.001)
}
"

# 6. Fit the model
fit_sim <- jags(
  data = jags_data_sim,
  parameters.to.save = c("IC50", "tau"),
  model.file = textConnection(jags_model_sim),
  n.chains = 3, n.iter = 15000, n.burnin = 3000, quiet = TRUE
)

# 7. Print Recovery Results
ic50_posterior <- fit_sim$BUGSoutput$summary["IC50", ]
tau_posterior <- fit_sim$BUGSoutput$summary["tau", ]

cat("\n=== PARAMETER RECOVERY ===\n")
cat("--- IC50 ---\n")
cat("True IC50:", true_ic50, "\n")
cat("Estimated Posterior Mean:", round(ic50_posterior["mean"], 4), "\n")
cat("95% Credible Interval: [", round(ic50_posterior["2.5%"], 4), "-", round(ic50_posterior["97.5%"], 4), "]\n\n")

cat("--- Precision (tau) ---\n")
cat("True Precision:", true_tau, "\n")
cat("Estimated Posterior Mean:", round(tau_posterior["mean"], 4), "\n")
cat("95% Credible Interval: [", round(tau_posterior["2.5%"], 4), "-", round(tau_posterior["97.5%"], 4), "]\n")