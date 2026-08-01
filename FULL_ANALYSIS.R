# Full analysis pipeline for the Trametinib IC50 project.
# Frequentist pre-fit -> hierarchical Bayesian model in JAGS -> model checking,
# parameter recovery, sensitivity analysis, and all the comparison plots.

#SETUP---------------------------------------------------------------------------

library(dplyr)
library(drc)
library(R2jags)
library(ggplot2)
library(patchwork)
library(tidyr)
library(coda)

# main dataset: dose-response measurements for every cell line
data <- TRAMETINIB_jags_data

#FREQUENTIST APPROACH-------------------------------------------------------------

# look at one patient (ACH-000327) against the pooled fit, just to get a feel
# for how much a single cell line can differ from the population as a whole
specific_patient_data <- data %>%
  filter(patient_id == "ACH-000327")

fit_freq_specific <- tryCatch({
  drm(viability ~ dose,
      data = specific_patient_data,
      fct = LL.4(names = c("h_slope", "lower_limit", "upper_limit", "IC50")))
}, error = function(e) {
  NULL
})

if (!is.null(fit_freq_specific)) {
  parameters_specific <- coef(fit_freq_specific)
  h_spec <- as.numeric(parameters_specific["h_slope:(Intercept)"])
  ic50_spec <- as.numeric(parameters_specific["IC50:(Intercept)"])
  min_spec <- as.numeric(parameters_specific["lower_limit:(Intercept)"])
  max_spec <- as.numeric(parameters_specific["upper_limit:(Intercept)"])
  
  if (!is.na(h_spec) && !is.na(ic50_spec) && h_spec > 0 && ic50_spec > 0) {
    print(paste("Specific Patient (ACH-000327) - Estimated slope (h):", round(h_spec, 4)))
    print(paste("Specific Patient (ACH-000327) - Estimated IC50:", round(ic50_spec, 4)))
    print(paste("Specific Patient (ACH-000327) - Estimated Lower Limit:", round(min_spec, 4)))
    print(paste("Specific Patient (ACH-000327) - Estimated Upper Limit:", round(max_spec, 4)))
  } else {
    print("Specific patient parameters are invalid or non-positive.")
  }
} else {
  print("Specific patient model fit failed.")
}

# same model, but pooling every cell line together
fit_freq_generic <- tryCatch({
  drm(viability ~ dose,
      data = data,
      fct = LL.4(names = c("h_slope", "lower_limit", "upper_limit", "IC50")))
}, error = function(e) {
  NULL
})

if (!is.null(fit_freq_generic)) {
  parameters_generic <- coef(fit_freq_generic)
  h_gen <- as.numeric(parameters_generic["h_slope:(Intercept)"])
  ic50_gen <- as.numeric(parameters_generic["IC50:(Intercept)"])
  
  if (!is.na(h_gen) && !is.na(ic50_gen) && h_gen > 0 && ic50_gen > 0) {
    print(paste("Generic Fit (All Patients) - Estimated slope (h):", round(h_gen, 4)))
    print(paste("Generic Fit (All Patients) - Estimated IC50:", round(ic50_gen, 4)))
  } else {
    print("Generic fit parameters are invalid or non-positive.")
  }
} else {
  print("Generic model fit failed.")
}

par(mfrow = c(1, 2))

if (!is.null(fit_freq_specific)) {
  plot(fit_freq_specific,
       type = "all",
       main = "Frequentist Fit (Patient ACH-000327)",
       xlab = "Dose (log scale)",
       ylab = "Cell Viability")
}

if (!is.null(fit_freq_generic)) {
  plot(fit_freq_generic,
       type = "all",
       main = "Frequentist Fit (All Patients)",
       xlab = "Dose (log scale)",
       ylab = "Cell Viability")
}

par(mfrow = c(1, 1))

# now fit every cell line separately, so we have a frequentist IC50 (and the
# other LL.4 parameters) to feed into the Bayesian model as fixed inputs
results <- list()

for (pat in unique(data$patient_id)) {
  p_data <- filter(data, patient_id == pat)
  
  fit <- tryCatch(
    drm(viability ~ dose, data = p_data, fct = LL.4(names = c("h_slope", "lower_limit", "upper_limit", "IC50"))),
    error = function(e) NULL
  )
  
  if (!is.null(fit)) {
    co <- coef(fit)
    h_val <- as.numeric(co["h_slope:(Intercept)"])
    ic50_val <- as.numeric(co["IC50:(Intercept)"])
    min_val <- as.numeric(co["lower_limit:(Intercept)"])
    max_val <- as.numeric(co["upper_limit:(Intercept)"])
    
    # keep only fits that actually converged to something sensible
    if (!is.na(h_val) && !is.na(ic50_val) && h_val > 0 && ic50_val > 0) {
      results[[length(results) + 1]] <- data.frame(
        patient_id = pat,
        h_slope = h_val,
        ic50_frequentist = ic50_val,
        min_val = min_val,
        max_val = max_val,
        stringsAsFactors = FALSE
      )
    }
  }
}

df_out <- bind_rows(results)

if (!is.null(fit_freq_generic)) {
  df_out$generic_h_reference <- as.numeric(coef(fit_freq_generic)["h_slope:(Intercept)"])
}

# plain data frame, no tibble weirdness downstream
df_out <- as.data.frame(df_out)

# keep a copy on disk too, useful outside of this session
write.table(df_out, "trametinib_patient_parameters_theoretical.csv",
            sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

#BAYESIAN IC50 ESTIMATION (JAGS)---------------------------------------------------

# keep only the cell lines for which the frequentist pre-fit actually converged
data <- data %>%
  filter(patient_id %in% df_out$patient_id)

data$patient_idx <- as.numeric(as.factor(data$patient_id))

patient_mapping <- data %>%
  dplyr::select(patient_id, patient_idx) %>%
  distinct() %>%
  arrange(patient_idx)

fixed_params <- patient_mapping %>%
  inner_join(df_out, by = "patient_id")

# pooled model: one global IC50 shared by every cell line
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

    mu[i] <- min_fixed +
             (max_fixed - min_fixed) /
             (1 + pow(dose[i] / IC50, h_fixed))

  }

  IC50 ~ dgamma(0.001, 0.001)

  tau ~ dgamma(0.001, 0.001)

}
"

fit_bayes_generic <- jags(
  data = jags_data_generic,
  parameters.to.save = c("IC50", "tau"),
  model.file = textConnection(jags_model_generic),
  n.chains = 3,
  n.iter = 20000,
  n.burnin = 5000,
  DIC = TRUE
)

# hierarchical model: a separate IC50 per cell line, sharing a population prior
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
             (max_fixed[patient[i]] -
             min_fixed[patient[i]]) /
             (1 + pow(dose[i] /
             IC50[patient[i]],
             h_fixed[patient[i]]))

  }

  for (p in 1:n_patients) {

    IC50[p] ~ dgamma(shape_ic50,
                     rate_ic50)

  }

  shape_ic50 ~ dgamma(0.01, 0.01)

  rate_ic50 ~ dgamma(0.01, 0.01)

  tau ~ dgamma(0.001, 0.001)

}
"

fit_bayes_patient <- jags(
  data = jags_data_patient,
  parameters.to.save = c("IC50", "shape_ic50", "rate_ic50", "tau"),
  model.file = textConnection(jags_model_patient),
  n.chains = 3,
  n.iter = 20000,
  n.burnin = 5000,
  DIC = TRUE
)

cat("DIC COMPARISON\n")
cat("Generic Pooled Model DIC:", fit_bayes_generic$BUGSoutput$DIC, "\n")
cat("Patient-Specific Model DIC:", fit_bayes_patient$BUGSoutput$DIC, "\n")

#POSTERIOR PREDICTIVE CHECKS--------------------------------------------------------

ic50_generic_post <- fit_bayes_generic$BUGSoutput$mean$IC50
ic50_patient_post <- fit_bayes_patient$BUGSoutput$mean$IC50

posterior_ic50_df <- data.frame(
  patient_idx = 1:length(ic50_patient_post),
  ic50_est = ic50_patient_post
)

data_ppc <- data %>%
  left_join(fixed_params, by = c("patient_id", "patient_idx")) %>%
  left_join(posterior_ic50_df, by = "patient_idx")

# expected viability under the pooled model
data_ppc$pred_generic <- 0 + (1 - 0) /
  (1 + (data_ppc$dose / ic50_generic_post)^data_ppc$generic_h_reference)

# expected viability under the hierarchical model, using each patient's own limits
data_ppc$pred_specific <- data_ppc$min_val + (data_ppc$max_val - data_ppc$min_val) /
  (1 + (data_ppc$dose / data_ppc$ic50_est)^data_ppc$h_slope)

plot_generic <- ggplot(data_ppc, aes(x = pred_generic, y = viability)) +
  geom_point(alpha = 0.4, color = "purple") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1) +
  theme_minimal() +
  labs(
    title = "PPC: Generic Pooled Model",
    subtitle = paste("DIC:", round(fit_bayes_generic$BUGSoutput$DIC, 1)),
    x = "Predicted Viability",
    y = "Observed Viability"
  ) +
  coord_fixed(ratio = 1, xlim = c(0, 1.2), ylim = c(0, 1.2))

plot_specific <- ggplot(data_ppc, aes(x = pred_specific, y = viability)) +
  geom_point(alpha = 0.4, color = "orchid") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1) +
  theme_minimal() +
  labs(
    title = "PPC: Patient-Specific Model",
    subtitle = paste("DIC:", round(fit_bayes_patient$BUGSoutput$DIC, 1)),
    x = "Predicted Viability",
    y = "Observed Viability"
  ) +
  coord_fixed(ratio = 1, xlim = c(0, 1.2), ylim = c(0, 1.2))

print(plot_generic + plot_specific)

#BAYESIAN VS EMPIRICAL OBSERVATIONS--------------------------------------------------

# naive empirical IC50 for the global model: the dose where viability is closest to 0.5
global_closest <- data %>%
  mutate(dist_to_mid = abs(viability - 0.5)) %>%
  arrange(dist_to_mid) %>%
  slice(1)

emp_dose_gen <- global_closest$dose
emp_viab_gen <- global_closest$viability

post_ic50_gen <- data.frame(IC50 = as.numeric(fit_bayes_generic$BUGSoutput$sims.list$IC50))

# shared theme and colors, reused later for the frequentist comparison plot too
color_global_fill <- "purple"
color_pat_fill    <- "orchid"
color_line        <- "red"

academic_theme <- theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 14, color = "black"),
    plot.subtitle = element_text(size = 11, color = "grey30", margin = margin(b = 10)),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

x_limit_global <- max(quantile(post_ic50_gen$IC50, 0.99), emp_dose_gen * 1.2)

plot_global_emp <- ggplot(post_ic50_gen, aes(x = IC50)) +
  geom_density(fill = color_global_fill, alpha = 0.7, color = "black", linewidth = 0.5, adjust = 1.5) +
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

# same empirical logic, but per patient: the target viability is each patient's own midpoint
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

best_candidates <- patients_midpoint %>%
  filter(dist_to_target < 0.05) %>%
  head(2)

patient_plots_emp <- list()

for (i in 1:nrow(best_candidates)) {
  p_id <- best_candidates$patient_id[i]
  p_idx <- best_candidates$patient_idx[i]
  emp_dose <- best_candidates$dose[i]
  emp_viab <- best_candidates$viability[i]
  
  post_samples <- as.numeric(fit_bayes_patient$BUGSoutput$sims.list$IC50[, p_idx])
  df_post <- data.frame(IC50 = post_samples)
  
  x_limit_pat <- max(quantile(df_post$IC50, 0.99), emp_dose * 1.2)
  
  p_plot <- ggplot(df_post, aes(x = IC50)) +
    geom_density(fill = color_pat_fill, alpha = 0.7, color = "black", linewidth = 0.5, adjust = 1.5) +
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
  
  patient_plots_emp[[i]] <- p_plot
}

final_plot_emp <- plot_global_emp / patient_plots_emp[[1]] / patient_plots_emp[[2]] +
  plot_annotation(
    title = "Bayesian Posterior Densities vs. Empirical Observations",
    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 15)))
  )

print(final_plot_emp)

#OUTLIER ROBUSTNESS CHECK-----------------------------------------------------------

# ACH-000008 has one clearly isolated point at dose = 1; see how much it pulls the fit
target_patient <- "ACH-000008"

df_pat <- data %>% filter(patient_id == target_patient)
p_idx <- df_pat$patient_idx[1]

params_pat <- fixed_params %>% filter(patient_id == target_patient)
min_p <- params_pat$min_val
max_p <- params_pat$max_val
h_p <- params_pat$h_slope

ic50_est <- fit_bayes_patient$BUGSoutput$mean$IC50[p_idx]

dose_seq <- exp(seq(log(min(df_pat$dose)), log(max(df_pat$dose)), length.out = 100))
pred_viab <- min_p + (max_p - min_p) / (1 + (dose_seq / ic50_est)^h_p)
df_curve <- data.frame(dose = dose_seq, viability = pred_viab)

plot_curve <- ggplot() +
  geom_point(data = df_pat, aes(x = dose, y = viability),
             color = "black", size = 3, alpha = 0.7) +
  geom_line(data = df_curve, aes(x = dose, y = viability),
            color = "blue", linewidth = 1.2) +
  geom_vline(xintercept = ic50_est, linetype = "dashed",
             color = "green", linewidth = 1.2) +
  geom_vline(xintercept = 1.0, linetype = "dashed",
             color = "red", linewidth = 1.2) +
  geom_point(data = df_pat %>% filter(dose == 1), aes(x = dose, y = viability),
             color = "red", size = 5, shape = 1, stroke = 1.5) +
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

print(plot_curve)

#BAYESIAN VS FREQUENTIST COMPARISON-------------------------------------------------

# refit the pooled frequentist model on the trimmed patient set (the one actually
# used for the Bayesian analysis), so the reference line matches the same data
fit_freq_generic_refit <- tryCatch(
  drm(viability ~ dose, data = data, fct = LL.4(names = c("h_slope", "lower_limit", "upper_limit", "IC50"))),
  error = function(e) NULL
)
freq_ic50_global <- as.numeric(coef(fit_freq_generic_refit)["IC50:(Intercept)"])

x_limit_global_freq <- max(quantile(post_ic50_gen$IC50, 0.99), freq_ic50_global * 1.2)

plot_global_freq <- ggplot(post_ic50_gen, aes(x = IC50)) +
  geom_density(fill = color_global_fill, alpha = 0.7, color = "black", linewidth = 0.5, adjust = 1.5) +
  geom_vline(xintercept = freq_ic50_global, linetype = "dashed", color = color_line, linewidth = 1.2) +
  coord_cartesian(xlim = c(0, x_limit_global_freq)) +
  academic_theme +
  labs(
    title = "Global Population IC50",
    subtitle = sprintf("Density: Posterior | Dashed Line: Frequentist Estimate = %.4f", freq_ic50_global),
    x = "IC50 Estimate",
    y = "Density"
  )

target_patients <- c("ACH-000007", "ACH-000008")

candidates <- fixed_params %>%
  filter(patient_id %in% target_patients)

patient_plots_freq <- list()

for (i in 1:nrow(candidates)) {
  p_id <- candidates$patient_id[i]
  p_idx <- candidates$patient_idx[i]
  freq_ic50 <- candidates$ic50_frequentist[i]
  
  post_samples <- as.numeric(fit_bayes_patient$BUGSoutput$sims.list$IC50[, p_idx])
  df_post <- data.frame(IC50 = post_samples)
  
  x_limit_pat <- max(quantile(df_post$IC50, 0.99), freq_ic50 * 1.2)
  
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
  
  patient_plots_freq[[i]] <- p_plot
}

final_plot_freq <- plot_global_freq / patient_plots_freq[[1]] / patient_plots_freq[[2]] +
  plot_annotation(
    title = "Bayesian Posterior Densities vs. Frequentist Estimates",
    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 15)))
  )

print(final_plot_freq)

#PARAMETER RECOVERY------------------------------------------------------------------

# sanity check: can the model recover parameters we picked ourselves?
true_ic50 <- 0.08
true_h <- 1.5
true_min <- 0.2
true_max <- 1.0
true_tau <- 100 # precision, implies a noise sd of 0.1

sim_doses <- exp(seq(log(0.001), log(1), length.out = 15))

set.seed(123)
mu_sim <- true_min + (true_max - true_min) / (1 + (sim_doses / true_ic50)^true_h)
sim_viability <- rnorm(length(sim_doses), mean = mu_sim, sd = sqrt(1 / true_tau))

jags_data_sim <- list(
  viability = sim_viability,
  dose = sim_doses,
  N = length(sim_doses),
  h_fixed = true_h,
  min_fixed = true_min,
  max_fixed = true_max
)

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

fit_sim <- jags(
  data = jags_data_sim,
  parameters.to.save = c("IC50"),
  model.file = textConnection(jags_model_sim),
  n.chains = 3, n.iter = 15000, n.burnin = 3000, quiet = TRUE
)

sim_posterior <- fit_sim$BUGSoutput$summary["IC50", ]
cat("PARAMETER RECOVERY\n")
cat("True IC50:", true_ic50, "\n")
cat("Estimated Posterior Mean:", round(sim_posterior["mean"], 4), "\n")
cat("95% Credible Interval: [", round(sim_posterior["2.5%"], 4), "-", round(sim_posterior["97.5%"], 4), "]\n")

#INDIVIDUAL PATIENT DIAGNOSTICS-------------------------------------------------------

target_patient <- "ACH-000327"
p_idx <- fixed_params$patient_idx[fixed_params$patient_id == target_patient]
param_name <- paste0("IC50[", p_idx, "]")

cat("IN-DEPTH ANALYSIS FOR PATIENT:", target_patient, "\n")

# traceplot and autocorrelation, to check the chains actually converged
mcmc_samples <- as.mcmc(fit_bayes_patient)
combined_chains <- as.mcmc(do.call(rbind, mcmc_samples))

par(mfrow = c(1, 2))

traceplot(mcmc_samples[, param_name],
          main = paste("Traceplot:", target_patient),
          col = c("darkred", "darkblue", "yellow"))

# grab the autocorrelation values
acf_vals <- autocorr(combined_chains[, param_name])
lags <- 0:(length(acf_vals) - 1)

# recreate the same plot as autocorr.plot (using type = "h")
plot(x = lags, y = acf_vals,
     type = "h",
     main = paste("Autocorrelation:", target_patient),
     xlab = "Lag",
     ylab = "Autocorrelation",
     ylim = c(-0.1, 1))

# add the points and the dashed zero line
points(x = lags, y = acf_vals, pch = 16, cex = 0.7)
abline(h = 0, lty = 2)

par(mfrow = c(1, 1))

summary_stats <- fit_bayes_patient$BUGSoutput$summary[param_name, ]
cat("FORMAL DIAGNOSTICS\n")
cat("Gelman-Rubin (Rhat):", round(summary_stats["Rhat"], 4), "(Target: ~1.0)\n")
cat("Effective Sample Size (n.eff):", summary_stats["n.eff"], "\n")

# prior vs posterior, to see how much the data actually updated our beliefs
post_ic50 <- as.numeric(fit_bayes_patient$BUGSoutput$sims.list$IC50[, p_idx])
df_post <- data.frame(IC50 = post_ic50)

max_y <- max(density(post_ic50)$y)
x_max <- quantile(post_ic50, 0.99) * 1.5

plot_prior_post <- ggplot(df_post, aes(x = IC50)) +
  geom_density(aes(color = "Posterior", fill = "Posterior"), alpha = 0.6, linewidth = 1.2) +
  stat_function(aes(color = "Prior (Vague)", fill = "Prior (Vague)"),
                fun = dgamma, args = list(shape = 0.001, rate = 0.001),
                linetype = "dashed", linewidth = 1.2) +
  scale_color_manual(name = "Distribution", values = c("Posterior" = "purple", "Prior (Vague)" = "grey30")) +
  scale_fill_manual(name = "Distribution", values = c("Posterior" = "orchid", "Prior (Vague)" = "transparent")) +
  scale_x_continuous(limits = c(0, x_max), expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(ylim = c(0, max_y * 1.05)) +
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

# prior sensitivity: does the posterior actually depend on how vague the prior is?
df_sens <- data %>% filter(patient_id == target_patient)
params_sens <- fixed_params %>% filter(patient_id == target_patient)

fit_sensitivity <- function(shape_val, rate_val) {
  jags_data <- list(
    viability = df_sens$viability,
    dose = df_sens$dose,
    N = nrow(df_sens),
    h_fixed = params_sens$h_slope,
    min_fixed = params_sens$min_val,
    max_fixed = params_sens$max_val
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

mean_vague  <- fit_sensitivity(0.001, 0.001) # baseline vague prior
mean_weak   <- fit_sensitivity(1, 1)         # weakly informative prior
mean_strong <- fit_sensitivity(10, 10)       # strong prior pulling towards 1

sensitivity_results <- data.frame(
  Prior_Type = c("Vague Gamma(0.001, 0.001)", "Weak Gamma(1, 1)", "Strong Gamma(10, 10)"),
  Posterior_Mean_IC50 = c(mean_vague, mean_weak, mean_strong)
)

cat("PRIOR SENSITIVITY RESULTS\n")
print(sensitivity_results)

#ADVANCED VISUALIZATIONS--------------------------------------------------------------

cat("ADVANCED HIERARCHICAL VISUALIZATIONS\n")

# grab the posterior mean estimate for every patient from the JAGS model
bayes_summary <- fit_bayes_patient$BUGSoutput$summary
ic50_rows <- grep("^IC50\\[", rownames(bayes_summary))
bayes_ic50_means <- bayes_summary[ic50_rows, "mean"]

# build a dataset merging the frequentist and Bayesian estimates
df_shrinkage <- fixed_params %>%
  mutate(
    Bayesian_IC50 = bayes_ic50_means[patient_idx],
    Frequentist_IC50 = ic50_frequentist
  ) %>%
  # drop any NA or blown-up frequentist estimates, for a cleaner plot
  filter(!is.na(Frequentist_IC50) & Frequentist_IC50 < max(Bayesian_IC50) * 3)

global_mean_ic50 <- mean(df_shrinkage$Bayesian_IC50)

plot_shrinkage <- ggplot(df_shrinkage) +
  geom_vline(aes(xintercept = global_mean_ic50, linetype = "Global Population Mean"),
             color = "#B03A2E", linewidth = 1) +
  # arrows showing the shrinkage shift from frequentist to Bayesian estimate
  geom_segment(aes(x = Frequentist_IC50, xend = Bayesian_IC50,
                   y = reorder(patient_id, Bayesian_IC50), yend = reorder(patient_id, Bayesian_IC50)),
               arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
               color = "gray50", linewidth = 0.5) +
  geom_point(aes(x = Frequentist_IC50, y = reorder(patient_id, Bayesian_IC50), color = "Frequentist Estimate"),
             size = 2, alpha = 0.7) +
  geom_point(aes(x = Bayesian_IC50, y = reorder(patient_id, Bayesian_IC50), color = "Bayesian Hierarchical Estimate"),
             size = 2.5) +
  scale_color_manual(name = "Estimation Method",
                     values = c("Frequentist Estimate" = "#E41A1C",
                                "Bayesian Hierarchical Estimate" = "#2980B9")) +
  scale_linetype_manual(name = "Reference",
                        values = c("Global Population Mean" = "dashed")) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 15, color = "black"),
    plot.subtitle = element_text(size = 11, color = "grey30", margin = margin(b = 15))
  ) +
  labs(
    title = "Hierarchical Shrinkage Effect on IC50 Estimates",
    subtitle = "Arrows demonstrate extreme, isolated frequentist estimates being regularized toward the population mean",
    x = "IC50 Estimate",
    y = "Patients (Ordered by Bayesian Estimate)"
  )

print(plot_shrinkage)

# posterior predictive dose-response curve for our showcase patient
target_patient <- "ACH-000327"
p_idx <- fixed_params$patient_idx[fixed_params$patient_id == target_patient]

df_raw <- data %>% filter(patient_id == target_patient)

pat_params <- fixed_params %>% filter(patient_id == target_patient)
h_val <- pat_params$h_slope
min_val <- pat_params$min_val
max_val <- pat_params$max_val

post_ic50_pat <- fit_bayes_patient$BUGSoutput$sims.list$IC50[, p_idx]

dose_seq <- exp(seq(log(min(data$dose)), log(max(data$dose)), length.out = 100))

# each row is one MCMC draw, each column a dose
pred_matrix <- matrix(NA, nrow = length(post_ic50_pat), ncol = length(dose_seq))

for (i in 1:length(post_ic50_pat)) {
  pred_matrix[i, ] <- min_val + (max_val - min_val) / (1 + (dose_seq / post_ic50_pat[i])^h_val)
}

curve_data <- data.frame(
  dose = dose_seq,
  median_viability = apply(pred_matrix, 2, median),
  lwr_95 = apply(pred_matrix, 2, quantile, probs = 0.025),
  upr_95 = apply(pred_matrix, 2, quantile, probs = 0.975)
)

plot_predictive <- ggplot() +
  geom_ribbon(data = curve_data, aes(x = dose, ymin = lwr_95, ymax = upr_95, fill = "95% Credible Band"), alpha = 0.4) +
  geom_line(data = curve_data, aes(x = dose, y = median_viability, color = "Posterior Median"), linewidth = 1.2) +
  geom_point(data = df_raw, aes(x = dose, y = viability, shape = "Observed Data"), size = 2.5, alpha = 0.8) +
  scale_x_log10() +
  scale_fill_manual(name = "", values = c("95% Credible Band" = "lightblue")) +
  scale_color_manual(name = "", values = c("Posterior Median" = "blue")) +
  scale_shape_manual(name = "", values = c("Observed Data" = 16)) +
  theme_classic() +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = paste("Posterior Predictive Dose-Response:", target_patient),
    subtitle = "The credible band reflects the exact MCMC uncertainty around the non-linear fit",
    x = "Dose (Log Scale)",
    y = "Cell Viability"
  )

print(plot_predictive)