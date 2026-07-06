# ============================================================
# PILOT RF MODEL - OXISOLS
# SOC total inference from electrochemical signatures
# and auxiliary soil covariates
# ============================================================

# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

packages <- c(
  "readxl", "readr", "dplyr", "stringr", "janitor",
  "ranger", "ggplot2", "purrr", "tidyr"
)

installed <- rownames(installed.packages())
to_install <- setdiff(packages, installed)

if (length(to_install) > 0) {
  install.packages(to_install)
}

library(readxl)
library(readr)
library(dplyr)
library(stringr)
library(janitor)
library(ranger)
library(ggplot2)
library(purrr)
library(tidyr)

# ------------------------------------------------------------
# 1. Input paths
# ------------------------------------------------------------

dir_soil <- "D:/1_PhD_CIENCIAS_AGRARIAS/7_CAPÍTULOS_TESIS/CAPÍTULO_I/1_MEDICIONES_SUELOS/LECTURAS_PTO_LOPEZ_PTO_GAITAN_CORREGIDAS/SOIL_COVARIATES"

dir_features <- "D:/1_PhD_CIENCIAS_AGRARIAS/7_CAPÍTULOS_TESIS/CAPÍTULO_I/1_MEDICIONES_SUELOS/LECTURAS_PTO_LOPEZ_PTO_GAITAN_CORREGIDAS/PROCESADOS_CORRIENTE/ANALISIS_VOLTAGRAMAS_FEATURES"

soil_file <- file.path(dir_soil, "Soil_Covariates.xlsx")
features_file <- file.path(dir_features, "FEATURES_PARA_CALIBRACION_COS_WIDE.csv")

dir_out <- file.path(dir_features, "RF_OX_PILOT_SOC_TOTAL")
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 2. Modeling options
# ------------------------------------------------------------
# Recommended first option:
# "a_only"  = use the standard saturated paste condition
# "b_only"  = use the +10% water condition
# "mean_ab" = average a and b features by sample ID

feature_strategy <- "a_only"

set.seed(123)

# ------------------------------------------------------------
# 3. Read soil covariates
# ------------------------------------------------------------

soil <- read_excel(
  soil_file,
  na = c("", "NA", "N.A.", "-", "NaN")
) %>%
  clean_names()

# ------------------------------------------------------------
# Standardize possible pH column names after clean_names()
# ------------------------------------------------------------

if ("p_h_colorimetric_field" %in% names(soil) &&
    !"ph_colorimetric_field" %in% names(soil)) {
  
  soil <- soil %>%
    rename(ph_colorimetric_field = p_h_colorimetric_field)
}

# Check column names
print(names(soil))

# Recalculate depth variables to avoid formula inconsistencies
soil <- soil %>%
  mutate(
    muestra_id = as.integer(muestra_id),
    sitio_norm = str_to_lower(sitio),
    sitio_norm = case_when(
      str_detect(sitio_norm, "she") ~ "sheefer",
      str_detect(sitio_norm, "nubes") ~ "las_nubes",
      str_detect(sitio_norm, "carimagua") ~ "carimagua",
      TRUE ~ sitio_norm
    ),
    order = str_to_upper(order),
    depthmid = (sup_limit + inferior_limit) / 2,
    thickness = inferior_limit - sup_limit,
    horizon = as.factor(horizon),
    horizon_group = case_when(
      str_detect(str_to_lower(as.character(horizon)), "^a|^ap$") ~ "A",
      str_detect(str_to_lower(as.character(horizon)), "^ab$") ~ "AB",
      str_detect(str_to_lower(as.character(horizon)), "^b") ~ "B",
      TRUE ~ "Other"
    ),
    horizon_group = as.factor(horizon_group)
  )

# ------------------------------------------------------------
# 4. Read electrochemical features
# ------------------------------------------------------------

features_raw <- read_csv(
  features_file,
  na = c("", "NA", "N.A.", "-", "NaN"),
  show_col_types = FALSE
) %>%
  clean_names()

features <- features_raw %>%
  mutate(
    muestra_id = as.integer(str_extract(codigo_muestra, "\\d+")),
    treatment_paste = str_extract(codigo_muestra, "[A-Za-z]+$"),
    treatment_paste = str_to_lower(treatment_paste),
    sitio_features_norm = str_to_lower(sitio),
    sitio_features_norm = case_when(
      str_detect(sitio_features_norm, "she") ~ "sheefer",
      str_detect(sitio_features_norm, "nubes") ~ "las_nubes",
      str_detect(sitio_features_norm, "carimagua") ~ "carimagua",
      TRUE ~ sitio_features_norm
    )
  )

# ------------------------------------------------------------
# 5. Select feature strategy: a_only, b_only, or mean_ab
# ------------------------------------------------------------

feature_cols <- names(features)[str_detect(names(features), "^g[0-3]_")]

if (feature_strategy == "a_only") {
  
  features_model <- features %>%
    filter(treatment_paste == "a") %>%
    select(muestra_id, sitio_features_norm, all_of(feature_cols))
  
} else if (feature_strategy == "b_only") {
  
  features_model <- features %>%
    filter(treatment_paste == "b") %>%
    select(muestra_id, sitio_features_norm, all_of(feature_cols))
  
} else if (feature_strategy == "mean_ab") {
  
  features_model <- features %>%
    group_by(muestra_id, sitio_features_norm) %>%
    summarise(
      across(all_of(feature_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
} else {
  stop("Invalid feature_strategy. Use: 'a_only', 'b_only', or 'mean_ab'.")
}

# ------------------------------------------------------------
# 6. Join soil covariates and electrochemical features
# ------------------------------------------------------------

data_all <- soil %>%
  left_join(features_model, by = "muestra_id") %>%
  mutate(
    site_match_check = sitio_norm == sitio_features_norm
  )

# Check site matching
print(data_all %>% select(muestra_id, sitio, sitio_norm, sitio_features_norm, site_match_check))

# ------------------------------------------------------------
# 7. Filter Oxisols
# ------------------------------------------------------------

data_ox <- data_all %>%
  filter(order == "OX")

# Important:
# SOC_total_lab is the response variable.
# It must NOT be included as predictor.

# ------------------------------------------------------------
# 8. Define predictors
# ------------------------------------------------------------
# For this pilot, avoid variables with no variability:
# order and landscape are constant in this OX subset.
# Longitude/Latitude are not used initially because n is very small.

soil_predictors <- c(
  "clay",
  "bulk_density",
  "ph_colorimetric_field",
  "tsample_paste",
  "depthmid",
  "thickness"
)

# Recommended electrochemical features for small n:
# use G2 and G3 core features. You can expand later.
#Change: g[0-3] if we need de four G.

electro_predictors <- names(data_ox)[
  str_detect(
    names(data_ox),
    "^g[23]_(di_mean_na|di_sd_na|peak_abs_magnitude_na|e_peak_abs_v|auc_abs_na_v|pendiente_)"
  )
]

predictors <- c(soil_predictors, electro_predictors)

# Keep only predictors that exist
predictors <- predictors[predictors %in% names(data_ox)]

# Remove predictors with all NA
predictors <- predictors[
  map_lgl(data_ox[predictors], ~ !all(is.na(.x)))
]

# Remove zero-variance numeric predictors
zero_var <- names(data_ox[predictors])[
  map_lgl(data_ox[predictors], function(x) {
    if (is.numeric(x)) {
      length(unique(na.omit(x))) <= 1
    } else {
      FALSE
    }
  })
]

predictors <- setdiff(predictors, zero_var)

cat("Predictors used:\n")
print(predictors)

# ------------------------------------------------------------
# 9. Prepare modeling table
# ------------------------------------------------------------

model_data <- data_ox %>%
  select(
    muestra_id, sitio, horizon, horizon_group,
    soc_total_lab,
    all_of(predictors)
  )

# Impute missing predictor values using median.
# For the current OX pilot, this is mainly a safety step.

impute_median <- function(x) {
  if (!is.numeric(x)) return(x)
  if (all(is.na(x))) return(x)
  x[is.na(x)] <- median(x, na.rm = TRUE)
  return(x)
}

model_data <- model_data %>%
  mutate(across(all_of(predictors), impute_median))

# Training data: only rows with laboratory SOC
train_data <- model_data %>%
  filter(!is.na(soc_total_lab))

# Prediction data: all 11 OX rows
pred_data <- model_data

cat("Number of OX rows:", nrow(model_data), "\n")
cat("Rows with SOC_total_lab for training:", nrow(train_data), "\n")
cat("Rows without SOC_total_lab, only prediction:", sum(is.na(model_data$soc_total_lab)), "\n")

# ------------------------------------------------------------
# 10. Fit exploratory Random Forest
# ------------------------------------------------------------
# Warning: n is very small. This is a pilot workflow, not final calibration.

rf_formula <- as.formula(
  paste("soc_total_lab ~", paste(predictors, collapse = " + "))
)

p <- length(predictors)
mtry_value <- max(1, floor(sqrt(p)))

rf_model <- ranger(
  formula = rf_formula,
  data = train_data,
  num.trees = 1000,
  mtry = mtry_value,
  min.node.size = 1,
  importance = "permutation",
  seed = 123
)

print(rf_model)

# ------------------------------------------------------------
# 11. Predict all OX samples
# ------------------------------------------------------------

pred_all <- predict(rf_model, data = pred_data)$predictions

results_all <- data_ox %>%
  select(
    any_of(c(
      "muestra_id", "sitio", "order", "coverage", "landscape",
      "horizon", "sup_limit", "inferior_limit", "depthmid", "thickness",
      "clay", "bulk_density", "soc_total_lab",
      "ph_colorimetric_field", "tsample_paste"
    ))
  ) %>%
  mutate(
    feature_strategy = feature_strategy,
    soc_pred_rf = pred_all,
    residual = soc_total_lab - soc_pred_rf,
    abs_error = abs(residual),
    row_use = if_else(
      is.na(soc_total_lab),
      "prediction_only_no_lab_reference",
      "training_reference_available"
    )
  )

write_csv(
  results_all,
  file.path(dir_out, paste0("RF_OX_predictions_all_11_", feature_strategy, ".csv"))
)

# ------------------------------------------------------------
# 12. Leave-one-out validation for measured SOC rows
# ------------------------------------------------------------
# This is exploratory because n is small.

loocv_results <- map_dfr(seq_len(nrow(train_data)), function(i) {
  
  test_i <- train_data[i, ]
  train_i <- train_data[-i, ]
  
  # Drop predictors with zero variance inside this fold
  predictors_i <- predictors[
    map_lgl(train_i[predictors], function(x) {
      if (is.numeric(x)) length(unique(na.omit(x))) > 1 else TRUE
    })
  ]
  
  rf_formula_i <- as.formula(
    paste("soc_total_lab ~", paste(predictors_i, collapse = " + "))
  )
  
  mtry_i <- max(1, floor(sqrt(length(predictors_i))))
  
  mod_i <- ranger(
    formula = rf_formula_i,
    data = train_i,
    num.trees = 1000,
    mtry = mtry_i,
    min.node.size = 1,
    seed = 123 + i
  )
  
  pred_i <- predict(mod_i, data = test_i)$predictions
  
  tibble(
    muestra_id = test_i$muestra_id,
    sitio = test_i$sitio,
    horizon = test_i$horizon,
    soc_obs = test_i$soc_total_lab,
    soc_pred_loocv = pred_i,
    residual_loocv = soc_obs - soc_pred_loocv,
    abs_error_loocv = abs(residual_loocv)
  )
})

rmse_loocv <- sqrt(mean(loocv_results$residual_loocv^2, na.rm = TRUE))
mae_loocv <- mean(abs(loocv_results$residual_loocv), na.rm = TRUE)

metrics_loocv <- tibble(
  feature_strategy = feature_strategy,
  n_train = nrow(train_data),
  n_predictors = length(predictors),
  RMSE_LOOCV = rmse_loocv,
  MAE_LOOCV = mae_loocv
)

write_csv(
  loocv_results,
  file.path(dir_out, paste0("RF_OX_LOOCV_results_", feature_strategy, ".csv"))
)

write_csv(
  metrics_loocv,
  file.path(dir_out, paste0("RF_OX_LOOCV_metrics_", feature_strategy, ".csv"))
)

print(metrics_loocv)

# ------------------------------------------------------------
# 13. Variable importance
# ------------------------------------------------------------

importance_df <- tibble(
  predictor = names(rf_model$variable.importance),
  importance = as.numeric(rf_model$variable.importance)
) %>%
  arrange(desc(importance))

write_csv(
  importance_df,
  file.path(dir_out, paste0("RF_OX_variable_importance_", feature_strategy, ".csv"))
)

# ------------------------------------------------------------
# 14. Plots
# ------------------------------------------------------------

# Observed vs predicted using LOOCV
p1 <- ggplot(loocv_results, aes(x = soc_obs, y = soc_pred_loocv)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = paste("Pilot RF OX - LOOCV observed vs predicted SOC", feature_strategy),
    x = "Observed SOC_total_lab",
    y = "Predicted SOC_total_lab"
  ) +
  theme_bw()

ggsave(
  filename = file.path(dir_out, paste0("RF_OX_LOOCV_observed_vs_predicted_", feature_strategy, ".png")),
  plot = p1,
  width = 7,
  height = 5,
  dpi = 300
)

# Residuals by horizon
p2 <- ggplot(loocv_results, aes(x = horizon, y = residual_loocv)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  labs(
    title = paste("Pilot RF OX - LOOCV residuals by horizon", feature_strategy),
    x = "Soil horizon",
    y = "Residual: observed - predicted"
  ) +
  theme_bw()

ggsave(
  filename = file.path(dir_out, paste0("RF_OX_LOOCV_residuals_by_horizon_", feature_strategy, ".png")),
  plot = p2,
  width = 7,
  height = 5,
  dpi = 300
)

# Variable importance
p3 <- importance_df %>%
  slice_max(order_by = importance, n = 20) %>%
  ggplot(aes(x = reorder(predictor, importance), y = importance)) +
  geom_col() +
  coord_flip() +
  labs(
    title = paste("Pilot RF OX - Variable importance", feature_strategy),
    x = "Predictor",
    y = "Permutation importance"
  ) +
  theme_bw()

ggsave(
  filename = file.path(dir_out, paste0("RF_OX_variable_importance_", feature_strategy, ".png")),
  plot = p3,
  width = 8,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------
# 15. Final message
# ------------------------------------------------------------

cat("\nProcessing finished.\n")
cat("Outputs saved in:\n")
cat(dir_out, "\n")