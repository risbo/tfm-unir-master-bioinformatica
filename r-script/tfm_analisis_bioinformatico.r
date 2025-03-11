# https://pubmed.ncbi.nlm.nih.gov/35947121/
# https://faseb.onlinelibrary.wiley.com/doi/epdf/10.1096/fj.202200190R

if (!requireNamespace("rpart", quietly = TRUE)) install.packages("rpart")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("glmnet", quietly = TRUE)) install.packages("glmnet")
if (!requireNamespace("limma", quietly = TRUE)) install.packages("limma")
if (!requireNamespace("Rtsne", quietly = TRUE)) install.packages("Rtsne")
if (!requireNamespace("caret", quietly = TRUE)) install.packages("caret")
if (!requireNamespace("FNN", quietly = TRUE)) install.packages("FNN")
if (!requireNamespace("plotly", quietly = TRUE)) install.packages("plotly")
if (!requireNamespace("uwot", quietly = TRUE)) install.packages("uwot")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("factoextra", quietly = TRUE)) install.packages("factoextra")
if (!requireNamespace("umap", quietly = TRUE)) install.packages("umap")
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
if (!requireNamespace("rpart.plot", quietly = TRUE)) install.packages("rpart.plot")
if (!requireNamespace("naivebayes", quietly = TRUE)) install.packages("naivebayes")
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("randomForest", quietly = TRUE)) install.packages("randomForest")
if (!requireNamespace("xgboost", quietly = TRUE)) install.packages("xgboost")
if (!requireNamespace("nnet", quietly = TRUE)) install.packages("nnet")
if (!requireNamespace("DESeq2", quietly = TRUE)) install.packages("DESeq2")

library(rpart) 
library(tidyr)
library(ggplot2)
library(glmnet)
library(Rtsne)
library(caret)
library(FNN)
library(plotly)
library(uwot) 
library(readr)
library(dplyr)
library(factoextra)
library(ggplot2)
library(naivebayes)

library(umap)
library(pacman)
library(rpart.plot)
library(limma)
library(pROC)
library(pheatmap)
library(randomForest)
library(xgboost)
library(nnet)
library(DESeq2)


###############################################################################
# Script en R: Lectura inicial de GSE277204_matrix_expression.txt
# Objetivo: Identificar cuántas columnas son numéricas vs. texto, y cuántos 
#           valores 0, NA o Inf aparecen.
###############################################################################

# 1) Cargar librerías necesarias (si no las tienes, instala)
# (En este caso, solo usaremos funciones base y 'dplyr' para conveniencia)

if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(dplyr)

# 2) Leer el archivo
file_path <- "GSE277204_matrix_expression.txt"  # Ajusta la ruta si es necesario
data <- read.delim(file_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# 3) Número de columnas totales
total_cols <- ncol(data)
cat("\nTotal de columnas:", total_cols, "\n")

# 4) Identificar columnas numéricas y no numéricas
is_num <- sapply(data, is.numeric)
num_cols <- sum(is_num)
string_cols <- sum(!is_num)

cat("\nColumnas numéricas:", num_cols)
cat("\nColumnas de texto / no-numéricas:", string_cols, "\n")

# 5) Contar valores 0, NA y ±Inf en cada columna numérica
# Creamos un data.frame con estadísticas de cada columna numérica
summary_df <- data.frame(
  Columna = names(data)[is_num],
  Tipo = "numérico",
  CantidadTotal = NA_integer_,
  CantidadNA = NA_integer_,
  CantidadCero = NA_integer_,
  CantidadInf = NA_integer_,
  stringsAsFactors = FALSE
)

# Llenamos la tabla con la información
for (i in seq_along(summary_df$Columna)) {
  colname <- summary_df$Columna[i]
  vec <- data[[colname]]
  
  # Cantidad total de filas en esa columna
  summary_df$CantidadTotal[i] <- length(vec)
  # Cantidad de NA
  summary_df$CantidadNA[i] <- sum(is.na(vec))
  # Cantidad de ceros
  summary_df$CantidadCero[i] <- sum(vec == 0, na.rm=TRUE)
  # Cantidad de Inf o -Inf
  summary_df$CantidadInf[i] <- sum(is.infinite(vec))
}

cat("\nEstadísticas para columnas numéricas:\n")
print(summary_df)

# 6) (Opcional) Contar cuántas filas tienen al menos un NA, 0 o Inf
# Por ejemplo:
rows_with_na <- sum(!complete.cases(data))
rows_with_inf <- sum(apply(data, 1, function(row) any(is.infinite(row))))
# Para contar filas con al menos un cero, habría que definir si 
# lo buscas en columnas numéricas solamente:
rows_with_zero <- sum(apply(data[is_num], 1, function(x) any(x == 0, na.rm=TRUE)))

cat("\nFilas con al menos un NA:", rows_with_na)
cat("\nFilas con al menos un Inf:", rows_with_inf)
cat("\nFilas con al menos un 0 (en columnas numéricas):", rows_with_zero, "\n")

cat("\nLectura y verificación finalizada.\n")


###############################################################################
# Script: Separar filas con NA o 0, y dividir en data_numeric y data_details
###############################################################################

# 1) Identificar columnas que parecen numéricas pero están en formato texto
cols_to_convert <- sapply(data, function(col) {
  if (is.character(col) || is.factor(col)) {  # Solo aplicar a columnas de texto o factores
    all(grepl("^[-0-9.e]+$", col[!is.na(col)]))  # Evitar errores con NAs
  } else {
    FALSE  # Mantener numéricas como están
  }
})

# 2) Convertir esas columnas a numérico
data[, cols_to_convert] <- lapply(data[, cols_to_convert], function(col) as.numeric(as.character(col)))

# 3) Verificar que las conversiones fueron exitosas
str(data)


# 1) Identificar columnas numéricas
is_num <- sapply(data, is.numeric)

# 2) Para cada fila, detectar si existe al menos un valor NA o un 0 en columnas numéricas
rows_na_or_zero <- apply(
  data[, is_num, drop = FALSE], 
  1, 
  function(row) any(is.na(row) | row == 0, na.rm = TRUE)
)

# 3) Crear null_data con esas filas
null_data <- data[rows_na_or_zero, ]

# 4) Remover esas filas del data original
data <- data[!rows_na_or_zero, ]

# 5) Crear data_numeric (solo columnas numéricas) con las filas limpias
data_numeric <- data[, is_num, drop = FALSE]

# 6) Crear data_details (columnas de texto u otros tipos)
data_details <- data[, !is_num, drop = FALSE]

head(data_numeric)
head(data_details)
head(data_scaled)




# Cargar librerías necesarias
library(tidyverse)
library(ggplot2)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(car)


# 1. Estadística Descriptiva
summary(data_numeric)

# 2. Matriz de correlación y visualización
cor_matrix <- cor(data_numeric, use = "pairwise.complete.obs")
corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.7)

# 3. Regresión Lineal Múltiple para evaluar la expresión de ARv7
# Definir la variable dependiente y las predictoras
lm_model <- lm(C42ARV3.fpkm ~ C42NC3.fpkm + C42ARV1.fpkm + RV1NC1.fpkm + RV1sh1.fpkm, data = data_numeric)
summary(lm_model)

# Diagnóstico del modelo
par(mfrow = c(2, 2))
plot(lm_model)

# 4. Análisis de Componentes Principales (PCA)
pca_res <- PCA(data_scaled, graph = FALSE)

# Visualización del PCA
fviz_eig(pca_res) # Scree plot para varianza explicada
fviz_pca_ind(pca_res, geom.ind = "point", col.ind = "cos2", repel = TRUE)

# 5. Visualización de datos
# Histograma de la expresión de ARv7
ggplot(data_numeric, aes(x = C42ARV3.fpkm)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(title = "Distribución de la expresión de ARv7", x = "FPKM", y = "Frecuencia")

# Boxplot de las expresiones de ARv7 en diferentes condiciones
ggplot(data_numeric, aes(x = factor(1), y = C42ARV3.fpkm)) +
  geom_boxplot(fill = "red", alpha = 0.5) +
  labs(title = "Boxplot de la expresión de ARv7", x = "Condición", y = "FPKM")

# Guardar resultados
write.csv(summary(data_numeric), "summary_statistics.csv")
write.csv(cor_matrix, "correlation_matrix.csv")




















# Cargar librerías necesarias
library(tidyverse)
library(ggplot2)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(car)
 

# 1. Estadística Descriptiva
summary(data_numeric)

# 2. Matriz de correlación y visualización
cor_matrix <- cor(data_numeric, use = "pairwise.complete.obs")
corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.7)

# 3. Regresión Lineal Múltiple para evaluar la expresión de ARv7
# Definir la variable dependiente y las predictoras
lm_model <- lm(C42ARV3.fpkm ~ C42NC3.fpkm + C42ARV1.fpkm + RV1NC1.fpkm + RV1sh1.fpkm, data = data_numeric)
summary(lm_model)

# Diagnóstico del modelo
par(mfrow = c(2, 2))
plot(lm_model)

# 4. Análisis de Componentes Principales (PCA)
pca_res <- PCA(data_scaled, graph = FALSE)

# Visualización del PCA
fviz_eig(pca_res) # Scree plot para varianza explicada
fviz_pca_ind(pca_res, geom.ind = "point", col.ind = "cos2", repel = TRUE)

# 5. Visualización de datos
# Histograma de la expresión de ARv7
ggplot(data_numeric, aes(x = C42ARV3.fpkm)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(title = "Distribución de la expresión de ARv7", x = "FPKM", y = "Frecuencia")

# Boxplot de las expresiones de ARv7 en diferentes condiciones
ggplot(data_numeric, aes(x = factor(1), y = C42ARV3.fpkm)) +
  geom_boxplot(fill = "red", alpha = 0.5) +
  labs(title = "Boxplot de la expresión de ARv7", x = "Condición", y = "FPKM")

# Gráfica de coeficientes del modelo de regresión
coeficientes <- coef(lm_model)[-1]  # Excluir el intercepto
variables <- names(coeficientes)
data_coef <- data.frame(variables, coeficientes)

ggplot(data_coef, aes(x = coeficientes, y = reorder(variables, coeficientes), fill = coeficientes > 0)) +
  geom_col() +
  scale_fill_manual(values = c("red", "blue")) +
  geom_text(aes(label = round(coeficientes, 3)), hjust = ifelse(coeficientes > 0, -0.1, 1.1)) +
  labs(title = "Impacto de las Variables en la Expresión de ARv7 (C42ARV3.fpkm)",
       x = "Coeficiente de Regresión", y = "Variables Predictoras") +
  theme_minimal()

# Guardar resultados
write.csv(summary(data_numeric), "summary_statistics.csv")
write.csv(cor_matrix, "correlation_matrix.csv")


# 6. Predicción de respuesta a Olaparib
# Calcular el percentil 75 de C42ARV3.fpkm como umbral para resistencia
umbral_resistencia <- quantile(data_numeric$C42ARV3.fpkm, 0.75)

# Función para predecir respuesta a Olaparib en un nuevo paciente
predecir_respuesta <- function(C42NC3, C42ARV1, RV1NC1, RV1sh1) {
  prediccion <- 0.3399 - 0.0538 * C42NC3 + 1.0133 * C42ARV1 - 0.0948 * RV1NC1 + 0.1144 * RV1sh1
  if (prediccion > umbral_resistencia) {
    return("Resistente a Olaparib")
  } else {
    return("Sensible a Olaparib")
  }
}

# Ejemplo con datos de un paciente
nuevo_paciente <- predecir_respuesta(C42NC3 = 2.5, C42ARV1 = 10, RV1NC1 = 1.5, RV1sh1 = 2)
print(nuevo_paciente)

























###############################################################################
# Resultado:
# - null_data: filas con al menos un NA o un 0 en columnas numéricas
# - data: data original sin esas filas
# - data_numeric: solo columnas numéricas de 'data'
# - data_details: columnas no numéricas de 'data'
###############################################################################



###############################################################################
# 0) Librerías y Datos (ya existen en tu script anterior)
###############################################################################
# Supongamos que ya corriste todo tu bloque de instalación y lectura de datos:
# ...
# data, data_numeric, data_details, null_data, etc.

###############################################################################
# 1) Escalado de data_numeric
###############################################################################
data_scaled <- scale(data_numeric)  # z-score: media=0, sd=1
# Verifica que data_numeric contenga SOLO columnas numéricas.

###############################################################################
# 2) PCA con prcomp
###############################################################################
# center=FALSE y scale.=FALSE en prcomp() porque ya hicimos scale() manualmente
pca_res <- prcomp(data_scaled, center=FALSE, scale.=FALSE)

# Mostrar resumen de varianza explicada
cat("\nResumen de PCA:\n")
print(summary(pca_res))

###############################################################################
# 3) Visualización de varianza explicada
###############################################################################
# A) screeplot base R
screeplot(pca_res, type="lines", main="Varianza Explicada por Componentes")

# B) fviz_eig de factoextra
library(factoextra)
fviz_eig(pca_res, addlabels=TRUE, ylim=c(0,50))

###############################################################################
# 4) Graficar los primeros componentes principales con ggplot (estilo preferido)
###############################################################################
# Crear data frame con coordenadas PCA (scores)
pca_df <- as.data.frame(pca_res$x)
pca_df$Name <- data_details$Name

# Si tienes alguna variable target (clase/grupo), por ejemplo, en data_details,
# podrías unirla a pca_df. Aquí dejo un ejemplo comentado:
# pca_df$Grupo <- data_details$MiColumnaFactor

# A) Plot PC1 vs PC2 con ggplot
# Gráfica PC1 vs PC2 con etiquetas
ggplot(pca_df, aes(x=PC1, y=PC2)) +
  geom_point(alpha=0.5, color="steelblue") +
  geom_text(aes(label=data_details$Name), size=3, vjust=-1) +
  ggtitle("PCA - Primeras dos componentes principales") +
  labs(x="PC1", y="PC2") +
  theme_minimal()

# Gráfica PC1 vs PC3 con etiquetas
ggplot(pca_df, aes(x=PC1, y=PC3)) +
  geom_point(alpha=0.5, color="tomato") +
  geom_text(aes(label=data_details$Name), size=3, vjust=-1) +
  ggtitle("PCA - Componentes 1 y 3") +
  labs(x="PC1", y="PC3") +
  theme_minimal()

cat("\nFin del script PCA con estilo preferido.\n")



# Cargar librerías necesarias
library(ggplot2)
library(cluster)
library(factoextra)

# Realizar clustering jerárquico usando los resultados PCA (primeros dos componentes)
set.seed(123)
dist_matrix <- dist(pca_df[,1:2], method = "euclidean") # Distancia euclídea
hc_model <- hclust(dist_matrix, method = "ward.D2") # Clustering jerárquico

# Elegir número óptimo de clusters (ej. 3 clusters)
k <- 3
clusters <- cutree(hc_model, k = k)
pca_df$Cluster <- factor(clusters)

# Gráfica visualizando los clusters usando PCA con etiquetas de nombres de genes
fviz_cluster(list(data = pca_df[,1:2], cluster = clusters),
             geom = "point",
             ellipse.type = "convex",
             main = "Clustering Jerárquico sobre PCA",
             xlab = "PC1",
             ylab = "PC2") +
  geom_text(aes(label = data_details$Name), size = 3, vjust = -1) +
  theme_minimal()




# Cargar librerías necesarias
library(ggplot2)
library(cluster)
library(factoextra)

# Realizar K-Means usando los resultados PCA (primeros dos componentes)
set.seed(123)
k <- 3
kmeans_model <- kmeans(pca_df[,1:2], centers = k, nstart = 25)
pca_df$Cluster <- factor(kmeans_model$cluster)

# Gráfica visualizando los clusters usando PCA con etiquetas de nombres de genes
fviz_cluster(kmeans_model,
             data = pca_df[,1:2],
             geom = "point",
             ellipse.type = "convex",
             main = "Clustering K-Means sobre PCA",
             xlab = "PC1",
             ylab = "PC2") +
  geom_text(aes(label = data_details$Name), size = 3, vjust = -1) +
  theme_minimal()














##  Paso A: Definir Variable Objetivo Supervisada
# Ejemplo práctico basado en nombres de columnas
library(dplyr)

# Extraer nombres de columnas para identificar condiciones
condiciones <- colnames(data_scaled)

# Ejemplo para crear etiquetas supervisadas basadas en nombre de columna
cond_labels <- ifelse(grepl("NC", condiciones), "Control",
                      ifelse(grepl("ARV", condiciones), "ARV7_treated", "Other"))

# Revisa las etiquetas generadas
print(cond_labels)



data_supervised <- as.data.frame(t(data_scaled))
colnames(data_supervised) <- paste0("gene_", seq_len(ncol(data_supervised)))
data_supervised$Condition <- factor(cond_labels)


##Paso B: División de datos (entrenamiento y validación)
library(caret)

set.seed(123)
indices <- createDataPartition(data_supervised$Condition, p = 0.7, list = FALSE)
train_data <- data_supervised[indices,]
test_data <- data_supervised[-indices,]

## Paso C: Entrenar un modelo supervisado (ejemplo con Random Forest)
library(randomForest)

# Entrenamiento
rf_model <- randomForest(Condition ~ ., data = train_data, importance = TRUE)

# Predicción
predicted <- predict(rf_model, test_data)

# Evaluación
confusionMatrix(predicted, test_data$Condition)












library(ggplot2)

# Seleccionar solo las columnas .fpkm
fpkm_cols <- grep(".fpkm$", colnames(data_numeric), value = TRUE)

# Convertir a formato largo para ggplot
data_long <- reshape2::melt(data_numeric[, fpkm_cols])

# Histograma de distribución de expresión génica
ggplot(data_long, aes(x = value)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.6) +
  scale_x_log10() +  # Transformación logarítmica para normalizar
  ggtitle("Distribución de Expresión Génica (FPKM)") +
  xlab("FPKM (log10)") +
  ylab("Frecuencia") +
  theme_minimal()





library(corrplot)

# Calcular matriz de correlación
cor_matrix <- cor(data_numeric, use = "pairwise.complete.obs")

# Graficar
corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.7)




library(dplyr)

# Agrupar por Name y calcular estadísticas sobre las columnas numéricas (FPKM y Read Count)
data_grouped <- data_numeric %>%
  group_by(data_details$Name) %>%
  summarise(across(where(is.numeric), list(
    mean = ~mean(.x, na.rm = TRUE),
    median = ~median(.x, na.rm = TRUE),
    sd = ~sd(.x, na.rm = TRUE)
  )))

# Ver los primeros resultados
head(data_grouped)




# Boxplot comparativo entre dos condiciones específicas
boxplot(data_numeric[, c("C42NC3.fpkm", "C42ARV1.fpkm")],
        main = "Comparación de FPKM: C42NC3 vs C42ARV1",
        col = c("#66C2A5", "#FC8D62"),
        ylab = "FPKM")


# 2. Prueba estadística (Wilcoxon) entre dos condiciones
wilcox_result <- wilcox.test(data_numeric$C42NC3.fpkm, data_numeric$C42ARV1.fpkm)
print(wilcox_result)




# 3. PCA para visualizar estructura general en los datos escalados
pca_result <- prcomp(data_scaled, scale. = TRUE)
pca_df <- as.data.frame(pca_result$x)
pca_df$gene <- rownames(pca_df)

# Gráfica PCA
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.7, color = "steelblue") +
  theme_minimal() +
  labs(title = "PCA de datos de expresión génica",
       x = "PC1",
       y = "PC2")










# 5. Análisis funcional: conteo de términos GO más frecuentes
library(stringr)

# Extraer y contar términos GO
go_terms <- unlist(str_split(data_details$GO, ";"))
go_counts <- sort(table(go_terms), decreasing = TRUE)
go_df <- as.data.frame(go_counts)
colnames(go_df) <- c("GO_Term", "Count")

# Visualizar top 15 términos GO más frecuentes
ggplot(head(go_df, 15), aes(x = reorder(GO_Term, Count), y = Count)) +
  geom_bar(stat = "identity", fill = "#8DA0CB") +
  coord_flip() +
  theme_light() +
  labs(title = "Top 15 términos GO más frecuentes",
       x = "Término GO",
       y = "Número de genes asociados")






# 6. Chi-cuadrado
# Supongamos que deseas evaluar asociación categórica entre cromosoma y dirección del gen
chi_data <- table(data_details$Chromosome, data_details$Direction)
chi_result <- chisq.test(chi_data)
print(chi_result)

# 7. Regresión lineal simple
# Ejemplo: Relación entre C42NC3.fpkm y C42ARV1.fpkm
lm_result <- lm(C42ARV1.fpkm ~ C42NC3.fpkm, data = data_numeric)
summary(lm_result)

# Visualizar regresión lineal
ggplot(data_numeric, aes(x = C42NC3.fpkm, y = C42ARV1.fpkm)) +
  geom_point(color="darkgreen") +
  geom_smooth(method="lm", color="red") +
  theme_minimal() +
  labs(title = "Regresión Lineal: C42NC3 vs C42ARV1",
       x = "C42NC3.fpkm",
       y = "C42ARV1.fpkm")

