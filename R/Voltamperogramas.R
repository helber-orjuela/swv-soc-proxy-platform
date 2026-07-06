# ============================================================
# ANALISIS DE VOLTAMPEROGRAMAS Y EXTRACCION DE FEATURES
# Sensor proxy COS - SWV
# Archivo de entrada: BASE_CONSOLIDADA_CORRIENTES_SWV.csv
# ============================================================

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)

# ============================================================
# 1. Cargar archivo consolidado
# ============================================================

# Opción 1: seleccionar manualmente el archivo
archivo_base <- file.choose()

# Opción 2: escribir la ruta directa
# archivo_base <- "D:/.../PROCESADOS_CORRIENTE/BASE_CONSOLIDADA_CORRIENTES_SWV.csv"

base <- read_csv(archivo_base, show_col_types = FALSE)

# Carpeta de salida
carpeta_salida <- file.path(dirname(archivo_base), "ANALISIS_VOLTAGRAMAS_FEATURES")
dir.create(carpeta_salida, showWarnings = FALSE)

carpeta_figuras <- file.path(carpeta_salida, "FIGURAS")
dir.create(carpeta_figuras, showWarnings = FALSE)

carpeta_individual <- file.path(carpeta_figuras, "VOLTAMPEROGRAMAS_INDIVIDUALES")
dir.create(carpeta_individual, showWarnings = FALSE)

# ============================================================
# 2. Verificación de columnas necesarias
# ============================================================

columnas_necesarias <- c(
  "Archivo", "Sitio", "Codigo_muestra",
  "Indice", "Tiempo_s", "E_plot_V",
  "G0_dI_nA", "G1_dI_nA", "G2_dI_nA", "G3_dI_nA"
)

faltantes <- setdiff(columnas_necesarias, names(base))

if (length(faltantes) > 0) {
  stop("Faltan columnas necesarias en la base: ", paste(faltantes, collapse = ", "))
}

# ============================================================
# 3. Convertir base de formato ancho a formato largo
# ============================================================

volt_long <- base %>%
  select(
    Archivo, Sitio, Codigo_muestra,
    Indice, Tiempo_s, E_plot_V,
    G0_dI_nA, G1_dI_nA, G2_dI_nA, G3_dI_nA
  ) %>%
  pivot_longer(
    cols = matches("^G[0-3]_dI_nA$"),
    names_to = "Ganancia",
    values_to = "dI_nA"
  ) %>%
  mutate(
    Ganancia = str_extract(Ganancia, "^G[0-3]")
  )

# ============================================================
# 4. Suavizado ligero para extracción de variables
# ============================================================
# Se conserva dI_nA como señal original.
# dI_smooth_nA se usa para reducir ruido en picos, áreas y pendientes.

suavizar <- function(x, k = 7) {
  if (length(x) < k) return(x)
  y <- as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))
  y[is.na(y)] <- x[is.na(y)]
  return(y)
}

volt_long <- volt_long %>%
  arrange(Archivo, Ganancia, E_plot_V) %>%
  group_by(Archivo, Sitio, Codigo_muestra, Ganancia) %>%
  mutate(
    dI_smooth_nA = suavizar(dI_nA, k = 7)
  ) %>%
  ungroup()

write_csv(
  volt_long,
  file.path(carpeta_salida, "BASE_LARGA_VOLTAGRAMAS.csv")
)

# ============================================================
# 5. Función para área bajo la curva por método trapezoidal
# ============================================================

trapz <- function(x, y) {
  orden <- order(x)
  x <- x[orden]
  y <- y[orden]
  
  ok <- complete.cases(x, y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 2) return(NA_real_)
  
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

# ============================================================
# 6. Extracción de variables electroquímicas por muestra y ganancia
# ============================================================

extraer_features <- function(df) {
  
  df <- df %>%
    arrange(E_plot_V)
  
  idx_abs <- which.max(abs(df$dI_smooth_nA))
  idx_min <- which.min(df$dI_smooth_nA)
  idx_max <- which.max(df$dI_smooth_nA)
  
  tibble(
    n_puntos = nrow(df),
    
    E_min_V = min(df$E_plot_V, na.rm = TRUE),
    E_max_V = max(df$E_plot_V, na.rm = TRUE),
    
    dI_mean_nA = mean(df$dI_smooth_nA, na.rm = TRUE),
    dI_sd_nA = sd(df$dI_smooth_nA, na.rm = TRUE),
    dI_min_nA = min(df$dI_smooth_nA, na.rm = TRUE),
    dI_max_nA = max(df$dI_smooth_nA, na.rm = TRUE),
    
    # Pico de mayor magnitud absoluta
    peak_abs_nA = df$dI_smooth_nA[idx_abs],
    peak_abs_magnitude_nA = abs(df$dI_smooth_nA[idx_abs]),
    E_peak_abs_V = df$E_plot_V[idx_abs],
    
    # Pico negativo, útil si la señal dominante es catódica/reductiva
    peak_min_nA = df$dI_smooth_nA[idx_min],
    E_peak_min_V = df$E_plot_V[idx_min],
    
    # Pico positivo, útil si aparece señal anódica/oxidativa
    peak_max_nA = df$dI_smooth_nA[idx_max],
    E_peak_max_V = df$E_plot_V[idx_max],
    
    # Área bajo la curva
    AUC_signed_nA_V = trapz(df$E_plot_V, df$dI_smooth_nA),
    AUC_abs_nA_V = trapz(df$E_plot_V, abs(df$dI_smooth_nA))
  )
}

features <- volt_long %>%
  group_by(Archivo, Sitio, Codigo_muestra, Ganancia) %>%
  group_modify(~ extraer_features(.x)) %>%
  ungroup()

write_csv(
  features,
  file.path(carpeta_salida, "FEATURES_ELECTROQUIMICAS_LONG.csv")
)

# ============================================================
# 7. Pendientes por intervalos de potencial
# ============================================================
# Puedes modificar los intervalos según el rango electroquímico de interés.

breaks_E <- seq(0, 1.5, by = 0.3)

labels_E <- paste0(
  sprintf("%.1f", head(breaks_E, -1)),
  "_",
  sprintf("%.1f", tail(breaks_E, -1)),
  "V"
)

calcular_pendiente <- function(df) {
  
  df <- df %>%
    filter(!is.na(E_plot_V), !is.na(dI_smooth_nA))
  
  if (nrow(df) < 3 || length(unique(df$E_plot_V)) < 2) {
    return(
      tibble(
        n_intervalo = nrow(df),
        E_ini_V = min(df$E_plot_V, na.rm = TRUE),
        E_fin_V = max(df$E_plot_V, na.rm = TRUE),
        pendiente_nA_V = NA_real_,
        intercepto_nA = NA_real_,
        r2_intervalo = NA_real_
      )
    )
  }
  
  modelo <- lm(dI_smooth_nA ~ E_plot_V, data = df)
  
  tibble(
    n_intervalo = nrow(df),
    E_ini_V = min(df$E_plot_V, na.rm = TRUE),
    E_fin_V = max(df$E_plot_V, na.rm = TRUE),
    pendiente_nA_V = coef(modelo)[2],
    intercepto_nA = coef(modelo)[1],
    r2_intervalo = summary(modelo)$r.squared
  )
}

pendientes_intervalos <- volt_long %>%
  mutate(
    intervalo_E = cut(
      E_plot_V,
      breaks = breaks_E,
      labels = labels_E,
      include.lowest = TRUE,
      right = FALSE
    )
  ) %>%
  filter(!is.na(intervalo_E)) %>%
  group_by(Archivo, Sitio, Codigo_muestra, Ganancia, intervalo_E) %>%
  group_modify(~ calcular_pendiente(.x)) %>%
  ungroup()

write_csv(
  pendientes_intervalos,
  file.path(carpeta_salida, "PENDIENTES_POR_INTERVALOS.csv")
)

# Versión ancha de pendientes
pendientes_wide <- pendientes_intervalos %>%
  select(
    Archivo, Sitio, Codigo_muestra, Ganancia,
    intervalo_E, pendiente_nA_V
  ) %>%
  pivot_wider(
    names_from = intervalo_E,
    values_from = pendiente_nA_V,
    names_prefix = "pendiente_"
  )

# ============================================================
# 8. Base final de variables para calibración de COS
# ============================================================

features_modelo_long <- features %>%
  left_join(
    pendientes_wide,
    by = c("Archivo", "Sitio", "Codigo_muestra", "Ganancia")
  )

write_csv(
  features_modelo_long,
  file.path(carpeta_salida, "FEATURES_PARA_CALIBRACION_COS_LONG.csv")
)

# Versión ancha: una fila por muestra, variables separadas por ganancia
cols_valores <- setdiff(
  names(features_modelo_long),
  c("Archivo", "Sitio", "Codigo_muestra", "Ganancia")
)

features_modelo_wide <- features_modelo_long %>%
  pivot_wider(
    id_cols = c(Archivo, Sitio, Codigo_muestra),
    names_from = Ganancia,
    values_from = all_of(cols_valores),
    names_glue = "{Ganancia}_{.value}"
  )

write_csv(
  features_modelo_wide,
  file.path(carpeta_salida, "FEATURES_PARA_CALIBRACION_COS_WIDE.csv")
)

# ============================================================
# 9. Gráfico general: voltamperogramas por sitio y ganancia
# ============================================================

p_general <- ggplot(
  volt_long,
  aes(x = E_plot_V, y = dI_nA, group = Archivo)
) +
  geom_line(alpha = 0.35, linewidth = 0.35) +
  facet_grid(Sitio ~ Ganancia, scales = "free_y") +
  labs(
    title = "Voltamperogramas diferenciales SWV por sitio y ganancia",
    subtitle = "Señal original: dI = Forward - Reverse",
    x = "Potential, E_plot (V)",
    y = expression(Delta*"I (nA)")
  ) +
  theme_bw()

ggsave(
  filename = file.path(carpeta_figuras, "01_VOLTAGRAMAS_TODAS_MUESTRAS_POR_SITIO_GANANCIA.pdf"),
  plot = p_general,
  width = 14,
  height = 9
)

# ============================================================
# 10. Gráfico promedio por sitio y ganancia
# ============================================================

promedio_sitio <- volt_long %>%
  group_by(Sitio, Ganancia, Indice, E_plot_V) %>%
  summarise(
    dI_mean_nA = mean(dI_smooth_nA, na.rm = TRUE),
    dI_sd_nA = sd(dI_smooth_nA, na.rm = TRUE),
    .groups = "drop"
  )

p_promedio <- ggplot(
  promedio_sitio,
  aes(x = E_plot_V, y = dI_mean_nA)
) +
  geom_line(linewidth = 0.65) +
  facet_grid(Sitio ~ Ganancia, scales = "free_y") +
  labs(
    title = "Voltamperogramas promedio por sitio y ganancia",
    subtitle = "Señal suavizada usada para extracción de variables",
    x = "Potential, E_plot (V)",
    y = expression(Delta*"I promedio (nA)")
  ) +
  theme_bw()

ggsave(
  filename = file.path(carpeta_figuras, "02_VOLTAGRAMAS_PROMEDIO_POR_SITIO_GANANCIA.pdf"),
  plot = p_promedio,
  width = 14,
  height = 9
)

# ============================================================
# 11. Voltamperogramas individuales por muestra
# ============================================================

muestras <- unique(volt_long$Archivo)

walk(muestras, function(muestra_i) {
  
  df_i <- volt_long %>%
    filter(Archivo == muestra_i)
  
  peaks_i <- features %>%
    filter(Archivo == muestra_i)
  
  p_i <- ggplot(df_i, aes(x = E_plot_V)) +
    geom_line(
      aes(y = dI_nA),
      alpha = 0.45,
      linewidth = 0.35
    ) +
    geom_line(
      aes(y = dI_smooth_nA),
      linewidth = 0.75
    ) +
    geom_point(
      data = peaks_i,
      aes(x = E_peak_abs_V, y = peak_abs_nA),
      inherit.aes = FALSE,
      size = 2
    ) +
    facet_wrap(~ Ganancia, scales = "free_y", ncol = 2) +
    labs(
      title = paste("Voltamperograma diferencial SWV -", muestra_i),
      subtitle = "Línea delgada: señal original | Línea gruesa: señal suavizada | Punto: pico de mayor magnitud",
      x = "Potential, E_plot (V)",
      y = expression(Delta*"I (nA)")
    ) +
    theme_bw()
  
  nombre_limpio <- str_replace_all(muestra_i, "[^A-Za-z0-9_\\-]", "_")
  
  ggsave(
    filename = file.path(carpeta_individual, paste0(nombre_limpio, "_voltamperograma.png")),
    plot = p_i,
    width = 10,
    height = 7,
    dpi = 300
  )
})

# ============================================================
# 12. Gráfico de picos por muestra y ganancia
# ============================================================

p_picos <- ggplot(
  features,
  aes(x = Archivo, y = peak_abs_magnitude_nA)
) +
  geom_col() +
  facet_wrap(~ Ganancia, scales = "free_y", ncol = 2) +
  coord_flip() +
  labs(
    title = "Magnitud del pico electroquímico por muestra y ganancia",
    x = "Muestra",
    y = "|Corriente pico| (nA)"
  ) +
  theme_bw()

ggsave(
  filename = file.path(carpeta_figuras, "03_MAGNITUD_PICO_POR_MUESTRA_GANANCIA.pdf"),
  plot = p_picos,
  width = 11,
  height = 9
)

# ============================================================
# 13. Mensaje final
# ============================================================

cat("Análisis terminado correctamente.\n")
cat("Resultados guardados en:\n")
cat(carpeta_salida, "\n")