

## Required packages

library(ggplot2)
library(car)
library(dplyr)
library(tidyr)
library(ggsignif)

# -------------------------------
# Function for C16/C16:1t Ratio Analysis
# -------------------------------
analizar_ratio <- function(datos) {
  # 1. Data Processing
  datos_ratio <- datos %>%
    mutate(Genotype = factor(Genotype, levels = c("WT", "ntrc", "NTRC_OE")),ratio = .data[["Ratio.C16.C16.1t"]]) %>%
    select(Genotype, ratio) %>%
    drop_na()
  
  # 2. Descriptive Analysis
  cat("\n----------------------------------\n")
  cat("Análisis para Ratio C16/C16:1t\n")
  cat("----------------------------------\n\n")
  print(summary(datos_ratio))
  
  # 3. Exploratory Plots
  p_jitter <- ggplot(datos_ratio, aes(x = Genotype, y = ratio, color = Genotype)) +
    geom_jitter(width = 0.2, size = 3, alpha = 0.7) +
    labs(title = "Ratio C16/C16:1t per Genotype", y = "Ratio C16/C16:1t") +
    theme_minimal() +
    theme(legend.position = "none")
  
  p_boxplot <- ggplot(datos_ratio, aes(x = Genotype, y = ratio, fill = Genotype)) +
    geom_boxplot(alpha = 0.7) +
    geom_point(position = position_jitter(width = 0.1), size = 2) +
    scale_fill_manual(values = c("WT" = "blue", "ntrc" = "red", "NTRC_OE" = "green")) +
    labs(y = "Ratio C16/C16:1t") +
    theme_minimal()
  
  # 4. Normality test per group
  normalidad <- datos_ratio %>%
    group_by(Genotype) %>%
    summarise(p_valor = shapiro.test(ratio)$p.value, .groups = "drop")
  cat("\nTest de normalidad (Shapiro-Wilk):\n"); print(normalidad)
  
  # 5. Homoscedasticity (Levene)
  levene <- leveneTest(ratio ~ Genotype, data = datos_ratio)
  cat("\nTest de homocedasticidad (Levene):\n"); print(levene)
  
  # 6. ANOVA or Kruskal-Wallis
  tukey <- NULL
  if (all(normalidad$p_valor > 0.05)) {
    cat("\nTodos los grupos son normales -> ANOVA\n")
    modelo <- aov(ratio ~ Genotype, data = datos_ratio)
    print(summary(modelo))
    if (summary(modelo)[[1]]$Pr[1] < 0.05) {
      cat("\nANOVA significativo -> Tukey HSD\n")
      tukey <- TukeyHSD(modelo)
      print(tukey)
    }
  } else {
    cat("\nDatos no normales -> Kruskal-Wallis\n")
    kruskal <- kruskal.test(ratio ~ Genotype, data = datos_ratio)
    print(kruskal)
    if (kruskal$p.value < 0.05) {
      cat("\nKruskal significativo -> Dunn's test\n")
      if (!require(FSA)) install.packages("FSA")
      dunn <- dunnTest(ratio ~ Genotype, data = datos_ratio, method = "bh")
      print(dunn)
    }
  }
  
  # 7. Residuals Normality & QQ plot
  p_qq <- NULL
  if (exists("modelo")) {
    cat("\nNormalidad de residuos:\n")
    residuos <- residuals(modelo)
    print(shapiro.test(residuos))
    p_qq <- ggplot(data.frame(residuos), aes(sample = residuos)) +
      stat_qq() + stat_qq_line() +
      ggtitle("QQ Plot de Residuos") + theme_minimal()
  }
  
  # 8. Adjusted p-values Table (Tukey)
  pval_df <- NULL
  if (!is.null(tukey)) {
    tukey_df <- as.data.frame(tukey$Genotype)
    tukey_df$comparacion <- rownames(tukey$Genotype)
    tukey_df$p_symbol <- cut(
      tukey_df$`p adj`,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("***", "**", "*", "ns")
    )
    pval_df <- tukey_df[, c("comparacion", "p adj", "p_symbol")]
  }
  
  # Final graphs
  if (!is.null(p_qq)) {
    print((p_jitter + p_boxplot) / p_qq)
  } else {
    print(p_jitter); print(p_boxplot)
  }
  
  # Result
  list(grafico_jitter    = p_jitter,
    grafico_boxplot   = p_boxplot,
    normalidad        = normalidad,
    homocedasticidad  = levene,
    modelo_anova      = if (exists("modelo")) modelo else NULL,
    test_posthoc      = tukey,
    p_values_posthoc  = pval_df)
}

# -------------------------------
# Auxiliary function to extract comparisons
# -------------------------------
extraer_comparaciones_significativas <- function(pval_df, alpha = 0.05) {
  if (is.null(pval_df) || nrow(pval_df) == 0) return(NULL)
  
  sig <- subset(pval_df, `p adj` < alpha)
  if (nrow(sig) == 0) return(NULL)
  
  # Split comparisons of the form "group1-group2" into pairs
  comps <- strsplit(as.character(sig$comparacion), "-")
  
  # Verify that all pairs are valid (two elements)
  comps_validas <- comps[sapply(comps, length) == 2]
  etiquetas <- sig$`p adj`[sapply(comps, length) == 2]
  
  if (length(comps_validas) == 0) return(NULL)
  
  # Transform p-values into significance symbols
  simbolos <- sapply(etiquetas, function(p) {
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    if (p < 0.1) return(".")
    return("ns")
  })
  
  list(pairs = lapply(comps_validas, function(x) c(x[1], x[2])),labels = simbolos)
}

# -------------------------------
# Function to plot with significance
# -------------------------------
graficar_con_lineas_significancia <- function(datos, pval_df, variable, titulo, subtitulo) {
  datos$Genotype <- factor(datos$Genotype, levels = c("WT", "ntrc", "NTRC_OE"))
  
  comps <- extraer_comparaciones_significativas(pval_df)
  print(comps)  
  
  p <- ggplot(datos, aes(x = Genotype, y = .data[[variable]], fill = Genotype)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.5) +
    scale_fill_manual(values = c("WT" = "blue", "ntrc" = "#FF6B6B", "NTRC_OE" = "green")) +
    labs(
      title = titulo,
      subtitle = subtitulo,
      x = "",
      y = variable) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text = element_text(size = 11),
      axis.title.y = element_text(size = 12))
  
  if (!is.null(comps)) {
    p <- p + geom_signif(
      comparisons = comps$pairs,
      annotations = comps$labels,
      tip_length = 0.01,
      vjust = 0.5,
      textsize = 4)
  }
  print(p)
}

# -------------------------------
# Real example execution
# -------------------------------
ratio <- read.table("Transpose_Ratio_C16_A.txt", header = TRUE, sep = "\t", dec = ".")

ratio$Genotype <- factor(ratio$Genotype, levels = c("WT", "ntrc", "NTRC_OE"))
ratio$Ratio.C16.C16.1t <- as.numeric(ratio$Ratio.C16.C16.1t)

# Analysis execution
resultados_ratio <- analizar_ratio(ratio)

# Plot with significance lines
graficar_con_lineas_significancia(
  datos     = ratio,
  pval_df = resultados_ratio$p_values_posthoc,
  variable  = "Ratio.C16.C16.1t",
  titulo    = "Significant differences in C16/C16:1t ratio",
  subtitulo = "Post-hoc Tukey HSD (one-way ANOVA)") +
  labs(y = "C16/C16:1t ratio")


###############################
###############################
###############################

# --------------------------------------------
# Fatty acid composition data loading and preprocessing
# --------------------------------------------

# Data loading
datos <- read.table(file = "Transpose_Ag_Data_Con.txt", header = TRUE, dec = ",") %>%
  mutate(Genotype = as.factor(Genotype),across(2:10, as.numeric)  )

# -------------------------------
# Function for acid-specific analysis
# -------------------------------
analizar_acido <- function(acido) {
  # 1. Data preprocessing
  datos_acido <- datos %>%
    dplyr::select(Genotype, all_of(acido)) %>%
    drop_na()
  
  colnames(datos_acido)[2] <- "concentracion"  # Estandarizar nombre
  
  # 2. Descriptive Analysis
  cat("\n----------------------------------\n")
  cat("Análisis para:", acido, "\n")
  cat("----------------------------------\n\n")
  
  print(summary(datos_acido))
  
  # 3. Exploratory graphs
  p_jitter <- ggplot(datos_acido, aes(x = Genotype, y = concentracion, color = Genotype)) +
    geom_jitter(width = 0.2, size = 3, alpha = 0.7) +
    labs(title = paste("Concentración de", acido), 
         y = paste(acido, "(mg/g DW)")) +
    theme_minimal() +
    theme(legend.position = "none")
  
  p_boxplot <- ggplot(datos_acido, aes(x = Genotype, y = concentracion, fill = Genotype)) +
    geom_boxplot(alpha = 0.7) +
    scale_fill_manual(values = c("WT" = "blue", "ntrc" = "red", "NTRC_OE" = "green")) +
    labs(y = paste(acido, "(mg/g DW)")) +
    theme_minimal()
  
  # 4. Normality test per group
  normalidad <- datos_acido %>%
    group_by(Genotype) %>%
    summarise(
      p_valor = shapiro.test(concentracion)$p.value,
      .groups = "drop")
  
  cat("\nTest de normalidad (Shapiro-Wilk):\n")
  print(normalidad)
  
  # 5. Homoscedasticity (Levene)
  levene <- leveneTest(concentracion ~ Genotype, data = datos_acido)
  cat("\nTest de homocedasticidad (Levene):\n")
  print(levene)
  
  # 6. ANOVA o Kruskal-Wallis 
  if (all(normalidad$p_valor > 0.05)) {
    cat("\nTodos los grupos son normales -> ANOVA\n")
    modelo <- aov(concentracion ~ Genotype, data = datos_acido)
    print(summary(modelo))
    
    # Post-hoc Tukey HSD test if ANOVA es significative
    if (summary(modelo)[[1]]$Pr[1] < 0.05) {
      cat("\nANOVA significativo -> Tukey HSD\n")
      tukey <- TukeyHSD(modelo)
      print(tukey)
    }
  } else {
    cat("\nDatos no normales -> Kruskal-Wallis\n")
    kruskal <- kruskal.test(concentracion ~ Genotype, data = datos_acido)
    print(kruskal)
    
    if (kruskal$p.value < 0.05) {
      cat("\nKruskal significativo -> Dunn's test\n")
      library(FSA)
      dunn <- dunnTest(concentracion ~ Genotype, data = datos_acido, method = "bh")
      print(dunn)
    }
  }
  
  pval_df <- NULL
  if (exists("tukey")) {
    # Convertir resultados de Tukey a dataframe
    tukey_df <- as.data.frame(tukey$Genotype)
    tukey_df$comparacion <- rownames(tukey$Genotype)
    tukey_df$p_symbol <- cut(
      tukey_df$`p adj`,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("***", "**", "*", "ns"))
    pval_df <- tukey_df[, c("comparacion", "p adj", "p_symbol")]
  } else if (exists("dunn")) {
    # Adaptar para Dunn's test (si usas Kruskal)
    pval_df <- data.frame(
      comparacion = paste0(dunn$res$Comparison),
      `p adj` = dunn$res$P.adj,
      p_symbol = cut(
        dunn$res$P.adj,
        breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
        labels = c("***", "**", "*", "ns")))
  }
  
  # 7. Residuals Normality
  if (exists("modelo")) {
    cat("\nNormalidad de residuos:\n")
    print(shapiro.test(residuals(modelo)))
  }
  
  # Return final plots & results
  return(list(
    grafico_jitter = p_jitter,
    grafico_boxplot = p_boxplot,
    normalidad = normalidad,
    homocedasticidad = levene,
    p_values_posthoc = pval_df))
}

# -------------------------------
# Run analysis for all acids
# -------------------------------
acidos_grasos <- colnames(datos)[2:10]  # Fatty acids as column names

resultados <- lapply(acidos_grasos, analizar_acido)
names(resultados) <- acidos_grasos

# Visualize plots (example for C16)
resultado_c16 <- analizar_acido("C16")
# Plot with significance lines
graficar_con_lineas_significancia(
  datos = datos %>% drop_na(C16), 
  pval_df = resultado_c16$p_values_posthoc,
  variable = "C16",
  titulo = "Significant differences in C16",
  subtitulo = "Post-hoc Tukey HSD (one-way ANOVA)") +
  labs(y = "C16 (mg/g DW)") 
