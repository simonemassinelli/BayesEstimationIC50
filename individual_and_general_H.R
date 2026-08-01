library(drc)
library(dplyr)

# Load the target dataset containing cell viability and dose-response measurements
data <- TRAMETINIB_jags_data
results <- list()

# Loop through each unique patient to fit individual dose-response curves
for (pat in unique(data$patient_id)) {
  p_data <- filter(data, patient_id == pat)
  
  # Fit the 4-parameter log-logistic model safely, handling any convergence errors
  fit <- tryCatch(
    drm(viability ~ dose, data = p_data, fct = LL.4(names = c("h_slope", "lower_limit", "upper_limit", "IC50"))), 
    error = function(e) NULL
  )
  
  # If the model converges, extract parameters and validate biological constraints
  if (!is.null(fit)) {
    co <- coef(fit)
    h_val <- as.numeric(co["h_slope:(Intercept)"])
    ic50_val <- as.numeric(co["IC50:(Intercept)"])
    
    # Keep only successful fits where both slope and IC50 are strictly positive
    if (!is.na(h_val) && !is.na(ic50_val) && h_val > 0 && ic50_val > 0) {
      results[[length(results) + 1]] <- data.frame(
        patient_id = pat,
        h_slope = h_val,
        ic50 = ic50_val
      )
    }
  }
}

# Combine all individual patient results into a single clean data frame
df_out <- bind_rows(results)

# Perform a generic fit on the entire dataset to establish a reference slope
fit_gen <- drm(viability ~ dose, data = data, fct = LL.4(names = c("h_slope", "lower_limit", "upper_limit", "IC50")))
df_out$generic_h_reference <- as.numeric(coef(fit_gen)["h_slope:(Intercept)"])

# Export the final parameters table to a CSV file and preview the first few rows
write.csv(df_out, "trametinib_patient_parameters.csv", row.names = FALSE)
print(head(df_out))