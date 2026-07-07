
# 1. LIBRERÍAS
# ------------------------------------------------------------------------------
# Solo ejecutar install.packages si no están instaladas
# install.packages(c("sf", "dplyr", "terra", "blockCV", "spdep", "spatialreg", "mapview"))

library(spData)     
library(sf)         
library(dplyr)      
library(spdep)      
library(spatialreg) 
library(blockCV)    
library(terra)      
library(mapview)




# 1. Descargar el archivo .gpkg desde el GitHub oficial
url_gpkg <- "https://github.com/Nowosad/spData/raw/master/inst/shapes/boston_tracts.gpkg"

download.file(url = url_gpkg, 
              destfile = "boston_tracts.gpkg", 
              mode = "wb")

# 2. Cargar el mapa en R
library(sf)
map_boston <- st_read("boston_tracts.gpkg")


map_boston <- st_transform(map_boston, 4326) # 4326 es el código para WGS84


names(map_boston)

library(spdep)



# 3. Ver un resumen de las conexiones
summary(nb)


boston_sf_prep <- map_boston %>%
  # Proyectamos a metros (Massachusetts Mainland EPSG:26986)
  st_transform(26986) %>% 
  # Creamos las variables siguiendo el esquema de Harrison y Rubinfeld (1978)
  mutate(
    # Variable dependiente
    price = CMEDV,                
    logprice = log(price),        
    
    # Variables independientes (Aseguramos que coincidan con tu fórmula)
    CRIM = CRIM,                  # Tasa de criminalidad
    RM = RM,                      # Habitaciones (al cuadrado suele usarse, pero RM es estándar)
    LSTAT = LSTAT,                # % población estatus bajo
    NOX = NOX,                    # Concentración de óxido de nitrógeno
    DIS = DIS,                    # Distancia a centros de empleo
    PTRATIO = PTRATIO,            # Ratio alumno/profesor
    TAX = TAX,                    # Impuesto a la propiedad
    CHAS = as.numeric(CHAS)       # Variable del río (Dummy: 1 si toca el río, 0 si no)
  ) %>%
  # Limpieza de seguridad: eliminamos filas con NAs en estas variables
  filter(!is.na(logprice))
# Comprobamos que las nuevas columnas existan
head(boston_sf_prep)

# 1. Crear la lista de vecinos (quién toca a quién)
nb <- poly2nb(boston_sf_prep, queen = TRUE)

# 2. Convertir la lista en una matriz de pesos (estandarizada por filas)
lw <- nb2listw(nb, style = "W")



# Al hacer el Test de Moran o el modelo, SIEMPRE añade zero.policy = TRUE
moran.test(boston_sf_prep$logprice, lw, zero.policy = TRUE)




library(spatialreg)

# Definimos la fórmula base (puedes añadir más variables si quieres)
formula_boston <- logprice ~ CRIM + RM + LSTAT + NOX + DIS + PTRATIO + TAX + CHAS
#formula_boston <- logprice ~ CRIM + RM + TAX + LSTAT + river + dist_centro

# Asegúrate de tener la matriz de pesos lista (usaremos la del set de entrenamiento)
# Si tuviste problemas de islas, recuerda usar zero.policy = TRUE

# Añadimos la columna que falta a ambos sets
#train_df$dist_centro <- train_df$DIS
#test_df$dist_centro  <- test_df$DIS


# 5. PREDICCIONES Y SOLUCIÓN DE ERRORES DE IDENTIFICACIÓN
# ------------------------------------------------------------------------------

# Predicción SAR
# Estimación del modelo SAR
m_sar <- lagsarlm(formula_boston, data = boston_sf_prep, listw = lw, zero.policy = TRUE)

summary(m_sar)

#Predicción SEM
m_sem <- errorsarlm(formula_boston, data = boston_sf_prep, listw = lw, zero.policy = TRUE)
summary(m_sem)

# Estimación del Spatial Durbin Error Model
m_sdem <- errorsarlm(formula_boston, data =  = boston_sf_prep, 
                    listw = lw, Durbin = TRUE, zero.policy = TRUE)

# Resumen rápido
summary(m_sdem)

# Comprobar NAs en las variables de la fórmula
colSums(is.na(train_df[, c("logprice", "CRIM", "RM", "TAX", "LSTAT", "river", "dist_centro")]))



# Nota: lmSLX usa zero.policy dentro del argumento, no fuera
m_slx <- lmSLX(formula_boston, data = train_df, listw = lw_train, zero.policy = TRUE)
summary(m_slx)


# Comparación de AIC
AIC(m_sem, m_sdem)

# LR Test: ¿Aportan algo las variables de los vecinos (WX)?
anova(m_sdem, m_sem)


# Función para calcular el RMSE
rmse_calc <- function(actual, predicha) {
  sqrt(mean((actual - predicha)^2))
}

# 1. Aseguramos que los IDs coincidan perfectamente
row.names(test_df) <- as.character(1:nrow(test_df))
nb_test <- poly2nb(test_df, queen = TRUE)
attr(nb_test, "region.id") <- row.names(test_df)
lw_test <- nb2listw(nb_test, style = "W", zero.policy = TRUE)

# 2. Predicción para el SEM (Tu modelo óptimo)
# Usamos solo el SEM para evitar el error de índices del SDM
pred_sem <- predict(m_sem, newdata = test_df, listw = lw_test, zero.policy = TRUE)

# 3. Cálculo del RMSE en dólares para el SEM
reales_usd <- exp(test_df$logprice) * 1000
predichos_sem_usd <- exp(as.numeric(pred_sem)) * 1000
rmse_sem <- sqrt(mean((reales_usd - predichos_sem_usd)^2))

cat("RMSE del modelo ganador (SEM): $", round(rmse_sem, 2), "\n")



# 1. Extraer las variables del set de prueba (sin la geometría)
# Identificamos las variables originales (X) y sus retardos (WX)
n_vars <- (length(m_sdem$coefficients) - 1) / 2
nombres_x <- names(m_sdem$coefficients)[2:(n_vars + 1)]

# 2. Crear la matriz X de prueba y su retardo WX
X_test_mat <- as.matrix(st_drop_geometry(test_df)[, nombres_x])
X_test_intercept <- cbind(1, X_test_mat) # Añadimos el intercepto
WX_test_mat <- lag.listw(lw_test, X_test_mat, zero.policy = TRUE)

# 3. Unir todo en una matriz de diseño: [Intercepto, X, WX]
X_final_test <- cbind(X_test_intercept, WX_test_mat)

# 4. Predicción Directa: Y_hat = X_final * Coeficientes
# (Nota: En SDEM la predicción se basa en la parte sistemática Xb + WXt)
pred_log_sdem <- X_final_test %*% m_sdem$coefficients

# 5. Transformación a dólares reales
reales_usd <- exp(test_df$logprice) * 1000
pred_usd_sdem <- exp(as.numeric(pred_log_sdem)) * 1000

# 6. Cálculo del RMSE
rmse_sdem <- sqrt(mean((reales_usd - pred_usd_sdem)^2))

cat("--------------------------------------------\n",
    "ERROR DE PREDICCIÓN SDEM (GANADOR AIC)\n",
    "RMSE Final: $", round(rmse_sdem, 2), "\n",
    "--------------------------------------------\n")



# 1. Ejecutar la predicción para el set de prueba
# Es vital pasarle la matriz lw_test que creamos para el set de validación
pred_sar <- predict(m_sar, newdata = test_df, listw = lw_test, zero.policy = TRUE)

# 2. Extraer los valores reales y predichos en dólares
# (Recordamos: exp() para revertir el logaritmo y *1000 por la escala de Boston)
reales_usd <- exp(test_df$logprice) * 1000
predichos_sar_usd <- exp(as.numeric(pred_sar)) * 1000

# 3. Calcular el RMSE
rmse_sar <- sqrt(mean((reales_usd - predichos_sar_usd)^2))

cat("--------------------------------------------\n",
    "RESULTADO MODELO SAR (Spatial Lag)\n",
    "RMSE Final: $", round(rmse_sar, 2), "\n",
    "--------------------------------------------\n")




# Modelo Lineal Simple (MCO) para comparar
m_mco <- lm(formula_boston, data = train_df)
pred_mco <- predict(m_mco, newdata = test_df)
rmse_mco <- sqrt(mean((reales_usd - (exp(pred_mco)*1000))^2))

cat("RMSE Regresión Lineal (MCO): $", round(rmse_mco, 2), "\n")
cat("Mejora del SEM respecto al MCO: $", round(rmse_mco - rmse_sem, 2), "\n")


# 1. Ajustamos el modelo MCO (Regresión Lineal Simple)
m_mco <- lm(formula_boston, data = train_df)

# 2. Obtenemos el AIC
aic_mco <- AIC(m_mco)

# 3. Obtenemos el Log-Likelihood
ll_mco <- logLik(m_mco)

# Mostramos los resultados
cat("AIC de MCO:", aic_mco, "\n")
cat("Log-Likelihood de MCO:", ll_mco, "\n")

# 1. Calculamos el Test de Moran para los residuos del MCO
# Usamos la matriz de pesos del entrenamiento
residuos_mco <- residuals(m_mco)
test_moran_mco <- lm.morantest(m_mco, lw_train)

print(test_moran_mco)

# Cargamos la librería correcta
library(spdep)

# 1. Extraemos los residuos del modelo ganador
residuos_sem <- residuals(m_sem)

# 2. Ejecutamos el Test de Moran (Monte Carlo es lo más robusto)
# Usamos nsim = 999 para una alta precisión en el p-valor
test_moran_sem <- moran.mc(residuos_sem, lw_train, nsim = 999, zero.policy = TRUE)
# 3. Vemos el resultado
print(test_moran_sem)

# 1. Obtenemos los residuos reales del set de prueba
# (Asegúrate de que 'pred_sem' y 'test_df' estén en la misma escala)
residuos_test <- test_df$logprice - as.numeric(pred_sem)

# 2. Ejecutamos el Test de Moran sobre esos residuos nuevos
test_moran_validacion <- moran.mc(residuos_test, lw_test, nsim = 999, zero.policy = TRUE)

# 3. Ver el resultado
print(test_moran_validacion)

# Función auxiliar para el cálculo (asegúrate de que esté cargada)
calc_rmse_usd <- function(actual_log, pred_obj) {
  reales <- exp(actual_log) * 1000
  predichos <- exp(as.numeric(pred_obj)) * 1000
  sqrt(mean((reales - predichos)^2))
}

# 1. Predicción para el modelo SAR (Spatial Lag Model)
pred_sar <- predict(m_sar, newdata = test_df, listw = lw_test, zero.policy = TRUE)
rmse_sar <- calc_rmse_usd(test_df$logprice, pred_sar)

# 2. Predicción para el modelo SDM (Spatial Durbin Model)
# Si predict() sigue fallando, usamos una alternativa para extraer la tendencia
pred_sdm_fixed <- predict(m_sdm, newdata = test_df, listw = lw_test, zero.policy = TRUE)
rmse_sdm <- calc_rmse_usd(test_df$logprice, pred_sdm_fixed)

# Mostramos los resultados
cat("RMSE Modelo SAR: $", round(rmse_sar, 2), "\n")
cat("RMSE Modelo SDM: $", round(rmse_sdm, 2), "\n")

# 1. Extraemos los residuos del modelo SAR
residuos_sar <- residuals(m_sar)

# 2. Ejecutamos la simulación de Monte Carlo (999 simulaciones)
# Usamos la matriz de pesos del entrenamiento (lw_train)
test_moran_sar <- moran.mc(residuos_sar, lw_train, nsim = 999, zero.policy = TRUE)

# 3. Mostramos el resultado
print(test_moran_sar)