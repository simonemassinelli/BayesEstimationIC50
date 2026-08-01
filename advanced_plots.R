library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork) 

cat("\n========================================================\n")
cat("ADVANCED HIERARCHICAL VISUALIZATIONS\n")
cat("========================================================\n")

# ==============================================================================
# 1. DATA PREPARATION (Extracting and Merging)
# ==============================================================================

# Estraiamo le stime medie a posteriori per tutti i pazienti dal modello JAGS
bayes_summary <- fit_patient$BUGSoutput$summary
ic50_rows <- grep("^IC50\\[", rownames(bayes_summary))
bayes_ic50_means <- bayes_summary[ic50_rows, "mean"]

# Creiamo un dataset unendo le stime Frequentiste e quelle Bayesiane
df_shrinkage <- fixed_params %>%
  mutate(
    Bayesian_IC50 = bayes_ic50_means[patient_idx],
    Frequentist_IC50 = ic50_frequentist
  ) %>%
  # Rimuoviamo eventuali NA o stime esplose frequentiste per pulizia visiva
  filter(!is.na(Frequentist_IC50) & Frequentist_IC50 < max(Bayesian_IC50) * 3)

# Calcoliamo la media globale per tracciare la linea di riferimento
global_mean_ic50 <- mean(df_shrinkage$Bayesian_IC50)

# ==============================================================================
# 2. SHRINKAGE PLOT (Academic Style)
# ==============================================================================

plot_shrinkage <- ggplot(df_shrinkage) +
  
  # Linea verticale (inserita in aes per forzarla nella legenda unificata)
  geom_vline(aes(xintercept = global_mean_ic50, linetype = "Global Population Mean"), 
             color = "#B03A2E", linewidth = 1) +
  
  # Frecce che mostrano lo spostamento (shrinkage)
  geom_segment(aes(x = Frequentist_IC50, xend = Bayesian_IC50, 
                   y = reorder(patient_id, Bayesian_IC50), yend = reorder(patient_id, Bayesian_IC50)), 
               arrow = arrow(length = unit(0.12, "cm"), type = "closed"), 
               color = "gray50", linewidth = 0.5) +
  
  # Pallini delle stime Frequentiste (Rossi)
  geom_point(aes(x = Frequentist_IC50, y = reorder(patient_id, Bayesian_IC50), color = "Frequentist Estimate"), 
             size = 2, alpha = 0.7) +
  
  # Pallini delle stime Bayesiane (Blu Acciaio)
  geom_point(aes(x = Bayesian_IC50, y = reorder(patient_id, Bayesian_IC50), color = "Bayesian Hierarchical Estimate"), 
             size = 2.5) +
  
  # Associazione dei colori
  scale_color_manual(name = "Estimation Method", 
                     values = c("Frequentist Estimate" = "#E41A1C", 
                                "Bayesian Hierarchical Estimate" = "#2980B9")) +
  
  # Associazione della linea nella legenda
  scale_linetype_manual(name = "Reference", 
                        values = c("Global Population Mean" = "dashed")) +
  
  # Tema accademico
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical", # Incolonna le due legende (metodo e linea)
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    axis.text.y = element_blank(), # Rimuove l'ammasso illeggibile dei nomi dei pazienti
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(), # Pulisce le righe orizzontali
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

# ==============================================================================
# 2. POSTERIOR PREDICTIVE CURVE (Dose-Response con 95% Credible Bands)
# ==============================================================================
# Usiamo il nostro paziente "da vetrina"
target_patient <- "ACH-000327"
p_idx <- fixed_params$patient_idx[fixed_params$patient_id == target_patient]

# Dati grezzi del paziente (i "pallini" sperimentali)
df_raw <- data %>% filter(patient_id == target_patient)

# Parametri fissi per la curva LL.4 (dal CSV teorico)
pat_params <- fixed_params %>% filter(patient_id == target_patient)
h_val <- pat_params$h_slope
min_val <- pat_params$min_val
max_val <- pat_params$max_val

# Campioni a posteriori dell'IC50 per questo paziente
post_ic50_pat <- fit_patient$BUGSoutput$sims.list$IC50[, p_idx]

# Creiamo una griglia di 100 dosi (in scala logaritmica) per tracciare una curva fluida
dose_seq <- exp(seq(log(min(data$dose)), log(max(data$dose)), length.out = 100))

# Matrice per salvare le predizioni: ogni riga è un'iterazione MCMC, ogni colonna è una dose
pred_matrix <- matrix(NA, nrow = length(post_ic50_pat), ncol = length(dose_seq))

# Calcoliamo la curva vitale attesa per ogni singolo campionamento a posteriori
for (i in 1:length(post_ic50_pat)) {
  pred_matrix[i, ] <- min_val + (max_val - min_val) / (1 + (dose_seq / post_ic50_pat[i])^h_val)
}

# Calcoliamo Mediana e intervalli di credibilità al 95% (2.5% e 97.5%)
curve_data <- data.frame(
  dose = dose_seq,
  median_viability = apply(pred_matrix, 2, median),
  lwr_95 = apply(pred_matrix, 2, quantile, probs = 0.025),
  upr_95 = apply(pred_matrix, 2, quantile, probs = 0.975)
)

# Creazione del grafico Dose-Risposta
plot_predictive <- ggplot() +
  # Nastro dell'incertezza Bayesiana (95% Credible Band)
  geom_ribbon(data = curve_data, aes(x = dose, ymin = lwr_95, ymax = upr_95, fill = "95% Credible Band"), alpha = 0.4) +
  
  # Curva mediana a posteriori
  geom_line(data = curve_data, aes(x = dose, y = median_viability, color = "Posterior Median"), linewidth = 1.2) +
  
  # Dati grezzi osservati
  geom_point(data = df_raw, aes(x = dose, y = viability, shape = "Observed Data"), size = 2.5, alpha = 0.8) +
  
  scale_x_log10() + # Scala X logaritmica obbligatoria per le dosi farmacologiche
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