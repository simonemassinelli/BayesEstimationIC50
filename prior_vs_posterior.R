library(ggplot2)

# 1. Seleziona il paziente
target_patient <- "ACH-000008"
p_idx <- fixed_params$patient_idx[fixed_params$patient_id == target_patient]

# 2. Estrai ESCLUSIVAMENTE i campioni a posteriori
post_ic50 <- as.numeric(fit_patient$BUGSoutput$sims.list$IC50[, p_idx])
df_post <- data.frame(IC50 = post_ic50)

# Calcola il limite massimo dell'asse Y per inquadrare perfettamente la campana
max_y <- max(density(post_ic50)$y)

# 3. Creiamo il grafico usando la formula matematica per la prior
plot_prior_post <- ggplot(df_post, aes(x = IC50)) +
  
  # Disegna la campana della Posterior
  geom_density(aes(fill = "Posterior", color = "Posterior"), alpha = 0.5, linewidth = 1) +
  
  # Disegna la linea teorica della Prior Gamma (senza simulare dati!)
  stat_function(aes(color = "Prior (Vague Gamma)"), 
                fun = dgamma, args = list(shape = 0.001, rate = 0.001), 
                linetype = "dashed", linewidth = 1) +
  
  scale_fill_manual(name = "Distribution", values = c("Posterior" = "steelblue")) +
  scale_color_manual(name = "Distribution", values = c("Posterior" = "darkblue", "Prior (Vague Gamma)" = "grey50")) +
  
  # Zoom sull'area di interesse, nascondendo il picco infinito a x=0 della Gamma
  coord_cartesian(xlim = c(0, quantile(post_ic50, 0.99) * 1.5), 
                  ylim = c(0, max_y * 1.1)) +
  theme_minimal() +
  labs(
    title = paste("Prior vs Posterior Learning:", target_patient),
    subtitle = "The flat dashed line is the vague prior; the sharp peak is the data-driven posterior",
    x = "IC50 Estimate",
    y = "Density"
  )

print(plot_prior_post)