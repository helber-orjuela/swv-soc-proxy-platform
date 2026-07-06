library(readr)
library(dplyr)
library(purrr)
library(stringr)

# =====================================================
# 1. Ruta de entrada
# =====================================================

carpeta <- "D:/1_PhD_CIENCIAS_AGRARIAS/7_CAPÍTULOS_TESIS/CAPÍTULO_I/1_MEDICIONES_SUELOS/LECTURAS_PTO_LOPEZ_PTO_GAITAN_CORREGIDAS"



# Verificar que la carpeta existe
if (!dir.exists(carpeta)) {
  stop("La carpeta indicada no existe. Revisa la ruta: ", carpeta)
}

carpeta_salida <- file.path(carpeta, "PROCESADOS_CORRIENTE")
dir.create(carpeta_salida, showWarnings = FALSE)

# =====================================================
# 2. Parámetros del sistema
# =====================================================

Vref <- 4.096
ADCmax <- 65535
Voffset <- 1.65

frecuencia <- 50
dt <- 1 / frecuencia

resistencias <- c(
  G0 = 1e3,     # 1 kΩ
  G1 = 10e3,    # 10 kΩ
  G2 = 100e3,   # 100 kΩ
  G3 = 1e6      # 1 MΩ
)

# =====================================================
# 3. Buscar archivos CSV originales
# =====================================================

archivos <- list.files(
  path = carpeta,
  pattern = "\\.csv$",
  full.names = TRUE
)

# Evitar procesar archivos ya generados
archivos <- archivos[
  !grepl("BASE_CONSOLIDADA|RESUMEN|corrientes|PROCESADOS", archivos, ignore.case = TRUE)
]

print(archivos)
cat("Número de archivos encontrados:", length(archivos), "\n")

if (length(archivos) == 0) {
  stop("No se encontraron archivos CSV en la carpeta indicada. Revisa la ruta o la extensión de los archivos.")
}

# =====================================================
# 4. Función para procesar cada archivo
# =====================================================

procesar_archivo <- function(archivo) {
  
  df <- read_csv(archivo, show_col_types = FALSE)
  
  nombre_archivo <- tools::file_path_sans_ext(basename(archivo))
  
  muestra <- str_remove(nombre_archivo, "_COS$")
  sitio <- str_remove(muestra, " \\(.*\\)")
  codigo <- str_extract(muestra, "\\((.*?)\\)")
  codigo <- str_remove_all(codigo, "\\(|\\)")
  
  # Agregar metadata manualmente para evitar errores de group_by
  df$Archivo <- nombre_archivo
  df$Sitio <- sitio
  df$Codigo_muestra <- codigo
  
  # Tiempo y potencial de graficación
  df$Tiempo_s <- (df$Indice - 1) * dt
  df$E_plot_V <- (df$`V_Forward(V)` + df$`V_Reverse(V)`) / 2
  
  # Reconstrucción de corrientes por ganancia
  for (g in names(resistencias)) {
    
    Rf <- resistencias[[g]]
    
    col_fwd <- paste0(g, "_FWD")
    col_rev <- paste0(g, "_REV")
    
    # Verificación de columnas
    if (!(col_fwd %in% names(df))) {
      stop("No se encontró la columna: ", col_fwd, " en el archivo: ", nombre_archivo)
    }
    
    if (!(col_rev %in% names(df))) {
      stop("No se encontró la columna: ", col_rev, " en el archivo: ", nombre_archivo)
    }
    
    # Voltaje ADC
    df[[paste0(g, "_FWD_Vadc")]] <- df[[col_fwd]] / ADCmax * Vref
    df[[paste0(g, "_REV_Vadc")]] <- df[[col_rev]] / ADCmax * Vref
    
    # Voltaje corregido por offset
    df[[paste0(g, "_FWD_Vsignal")]] <- df[[paste0(g, "_FWD_Vadc")]] - Voffset
    df[[paste0(g, "_REV_Vsignal")]] <- df[[paste0(g, "_REV_Vadc")]] - Voffset
    
    # Corriente en amperios
    df[[paste0(g, "_FWD_A")]] <- df[[paste0(g, "_FWD_Vsignal")]] / Rf
    df[[paste0(g, "_REV_A")]] <- df[[paste0(g, "_REV_Vsignal")]] / Rf
    
    # Corriente en nanoamperios
    df[[paste0(g, "_FWD_nA")]] <- df[[paste0(g, "_FWD_A")]] * 1e9
    df[[paste0(g, "_REV_nA")]] <- df[[paste0(g, "_REV_A")]] * 1e9
    
    # Diferencia RAW
    df[[paste0(g, "_dRAW")]] <- df[[col_fwd]] - df[[col_rev]]
    
    # Corriente diferencial SWV en nanoamperios
    df[[paste0(g, "_dI_nA")]] <- df[[paste0(g, "_FWD_nA")]] - df[[paste0(g, "_REV_nA")]]
  }
  
  # Guardar archivo individual procesado
  salida_individual <- file.path(
    carpeta_salida,
    paste0(nombre_archivo, "_corrientes.csv")
  )
  
  write_csv(df, salida_individual)
  
  return(df)
}

# =====================================================
# 5. Procesar todos los archivos
# =====================================================

base_corrientes <- map_dfr(archivos, procesar_archivo)

# Verificación de columnas creadas
print(names(base_corrientes))

if (!all(c("Archivo", "Sitio", "Codigo_muestra") %in% names(base_corrientes))) {
  stop("No se crearon correctamente las columnas Archivo, Sitio y Codigo_muestra.")
}

# Guardar base consolidada
write_csv(
  base_corrientes,
  file.path(carpeta_salida, "BASE_CONSOLIDADA_CORRIENTES_SWV.csv")
)

# =====================================================
# 6. Resumen por muestra
# =====================================================

resumen <- base_corrientes %>%
  group_by(Archivo, Sitio, Codigo_muestra) %>%
  summarise(
    n_puntos = n(),
    tiempo_total_s = max(Tiempo_s, na.rm = TRUE),
    
    G0_dI_min_nA = min(G0_dI_nA, na.rm = TRUE),
    G0_dI_max_nA = max(G0_dI_nA, na.rm = TRUE),
    G0_dI_mean_nA = mean(G0_dI_nA, na.rm = TRUE),
    
    G1_dI_min_nA = min(G1_dI_nA, na.rm = TRUE),
    G1_dI_max_nA = max(G1_dI_nA, na.rm = TRUE),
    G1_dI_mean_nA = mean(G1_dI_nA, na.rm = TRUE),
    
    G2_dI_min_nA = min(G2_dI_nA, na.rm = TRUE),
    G2_dI_max_nA = max(G2_dI_nA, na.rm = TRUE),
    G2_dI_mean_nA = mean(G2_dI_nA, na.rm = TRUE),
    
    G3_dI_min_nA = min(G3_dI_nA, na.rm = TRUE),
    G3_dI_max_nA = max(G3_dI_nA, na.rm = TRUE),
    G3_dI_mean_nA = mean(G3_dI_nA, na.rm = TRUE),
    
    .groups = "drop"
  )

write_csv(
  resumen,
  file.path(carpeta_salida, "RESUMEN_CORRIENTES_POR_MUESTRA.csv")
)

cat("Procesamiento finalizado correctamente.\n")