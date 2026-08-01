library(coda)

# 1. Estrai le catene MCMC
mcmc_samples <- as.mcmc(fit_patient)

# NUOVA RIGA: Unisci le 3 catene in un'unica catena per l'autocorrelazione
combined_chains <- as.mcmc(do.call(rbind, mcmc_samples))

# Seleziona il paziente
target_patient <- "ACH-000008"
p_idx <- fixed_params$patient_idx[fixed_params$patient_id == target_patient]
param_name <- paste0("IC50[", p_idx, "]")

# Imposta la griglia 1x2
par(mfrow = c(1, 2))

# 1. Traceplot (Usiamo le catene separate per vederle di colori diversi)
traceplot(mcmc_samples[, param_name], 
          main = paste("Traceplot:", target_patient), 
          col = c("#1b9e77", "#d95f02", "#7570b3"))

# 2. Autocorrelation plot (Usiamo le catene UNITE così disegna 1 solo grafico)
autocorr.plot(combined_chains[, param_name], 
              main = paste("Autocorrelation:", target_patient),
              auto.layout = FALSE)

# Resetta la griglia grafica
par(mfrow = c(1, 1))

# 3. Calcola e stampa le statistiche
summary_stats <- fit_patient$BUGSoutput$summary[param_name, ]
cat("\n=== FORMAL DIAGNOSTICS FOR", target_patient, "===\n")
cat("Gelman-Rubin (Rhat):", round(summary_stats["Rhat"], 4), "(Target: ~1.0)\n")
cat("Effective Sample Size (n.eff):", summary_stats["n.eff"], "\n")