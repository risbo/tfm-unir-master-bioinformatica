# TFM - Farmacogenómica: Una Aliada Estratégica en la Salud Pública y la Medicina de Precisión

Este repositorio contiene el código en R utilizado en mi Trabajo Fin de Máster en Bioinformática, cuyo objetivo es analizar la relación entre la expresión de la variante de splicing ARv7 y la respuesta terapéutica en cáncer de próstata. El estudio integra técnicas de análisis estadístico, machine learning, PCA y clustering para identificar biomarcadores predictivos que optimicen la selección de tratamientos en medicina personalizada.

Además, se ha desarrollado una plataforma web (Binary DNA) que centraliza las bases de datos clave en farmacogenómica, facilitando el acceso a información crítica para el personal sanitario. La URL de la plataforma es: [http://binaricode.com](http://binaricode.com).

## Contenido

- **Código en R:**  
  Scripts para la lectura y preprocesamiento de datos, análisis descriptivo, regresión lineal, PCA, clustering, modelos supervisados (Random Forest) y visualización de resultados.
  
- **Análisis Estadístico y Machine Learning:**  
  Implementación de técnicas para explorar la expresión de ARv7 y evaluar la influencia de variables predictoras en la respuesta terapéutica.
  
- **Visualizaciones:**  
  Generación de gráficos (histogramas, boxplots, diagramas de dispersión, scree plots, visualización de PCA, clustering, etc.) para apoyar la interpretación de los datos.
  
- **Plataforma Web Binary DNA:**  
  La web centraliza bases de datos relevantes (PharmGKB, DrugBank, dbSNP, ClinVar, 1000 Genomes Project, GTR, HGMD, gnomAD) y ofrece información sobre empresas colombianas especializadas en genómica y diagnóstico genético.

## Cómo Ejecutar el Código

### Requisitos
- R (versión 4.0 o superior)
- RStudio (opcional, pero recomendado)
- Conexión a Internet para la instalación de paquetes

### Instalación de Dependencias
El script inicia verificando la existencia de las librerías necesarias y las instala si es necesario. Para ejecutar el código, simplemente abre el archivo `r-script/tfm_analisis_bioinformatico.r` (o el nombre que le hayas asignado) en R o RStudio y ejecútalo línea por línea o en conjunto.

### Ejecución
1. Clona el repositorio:
   ```bash
   gh repo clone risbo/tfm-unir-master-bioinformatica
