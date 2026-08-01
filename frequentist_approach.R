library(drc)
library(dplyr)

# Load the target dataset containing cell viability and dose-response measurements
data <- TRAMETINIB_jags_data

# --- Part 1: Visualizing a Specific Patient vs Generic Fit ---

# Filter data for the specific patient
specific_patient_data <- data %>% 
  filter(patient_id == "ACH-000327")

# Fit the 4-parameter model safely for the specific patient
fit_specific <- tryCatch({
  drm(viability ~ dose, 
      data = specific_patient_data, 
      fct = LL.4(names = c("h_slope", "lower_limit", "upper_limit", "IC50")))
}, error = function(e) {
  NULL
})

# Validate and extract all parameters for the specific patient if the fit succeeded
if (!is.null(fit_specific)) {
  parameters_specific <- coef(fit_specific)
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

# Fit the 4-parameter model for all patients aggregated
fit_generic <- tryCatch({
  drm(viability ~ dose, 
      data = data, 
      fct = LL.4(names = c("h_slope", "lower_limit", "upper_limit", "IC50")))
}, error = function(e) {
  NULL
})

# Validate and extract parameters for the generic fit if it succeeded
if (!is.null(fit_generic)) {
  parameters_generic <- coef(fit_generic)
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

# Set plot layout to 1 row and 2 columns to compare both fits side by side
par(mfrow = c(1, 2))

if (!is.null(fit_specific)) {
  plot(fit_specific, 
       type = "all", 
       main = "Frequentist Fit (Patient ACH-000327)",
       xlab = "Dose (log scale)", 
       ylab = "Cell Viability")
}

if (!is.null(fit_generic)) {
  plot(fit_generic, 
       type = "all", 
       main = "Frequentist Fit (All Patients)",
       xlab = "Dose (log scale)", 
       ylab = "Cell Viability")
}

# Reset plot layout back to default single window
par(mfrow = c(1, 1))


# --- Part 2: Extracting All Parameters for Bayesian Analysis ---

results <- list()

# Loop through each unique patient to fit individual dose-response curves
for (pat in unique(data$patient_id)) {
  p_data <- filter(data, patient_id == pat)
  
  # Fit the 4-parameter log-logistic model safely
  fit <- tryCatch(
    drm(viability ~ dose, data = p_data, fct = LL.4(names = c("h_slope", "lower_limit", "upper_limit", "IC50"))), 
    error = function(e) NULL
  )
  
  # If the model converges, extract all four parameters
  if (!is.null(fit)) {
    co <- coef(fit)
    h_val <- as.numeric(co["h_slope:(Intercept)"])
    ic50_val <- as.numeric(co["IC50:(Intercept)"])
    min_val <- as.numeric(co["lower_limit:(Intercept)"])
    max_val <- as.numeric(co["upper_limit:(Intercept)"])
    
    # Keep only successful fits where both slope and IC50 are strictly positive
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

# Combine all individual patient results into a single clean data frame
df_out <- bind_rows(results)

# Append the generic slope reference calculated earlier
if (!is.null(fit_generic)) {
  df_out$generic_h_reference <- as.numeric(coef(fit_generic)["h_slope:(Intercept)"])
}

# Forziamo il dataframe in un formato base standard per distruggere eventuali strutture tibble
df_out <- as.data.frame(df_out)

# Export the final parameters table forcing explicit column names and removing quotes
write.table(df_out, "trametinib_patient_parameters_theoretical.csv", 
            sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

print("Parameter extraction completed and saved to CSV with explicit headers.")