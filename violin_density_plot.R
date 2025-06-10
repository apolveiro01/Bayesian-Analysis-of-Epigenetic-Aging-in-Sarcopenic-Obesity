# Required packages
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")

library(ggplot2)
library(dplyr)
library(tidyr)

# Load data
dados <- read.csv2("https://raw.githubusercontent.com/apolveiro01/Bayesian-Analysis-of-Epigenetic-Aging-in-Sarcopenic-Obesity/main/dados_epi.csv", header = TRUE)

# Define factors and order of groups
dados$gp <- factor(dados$gp, levels = c("IV", "I", "II", "III"))

# List of original variables
variaveis <- c(
  "EEAA","IEAA_Hannum", "DNAmTL", "DNAmFitAge",
  "DNAmAgeHannum","DNAmAgeHorvath", 
  "DNAmGrimAge", "DNAmGrimAge2",
  "Epigenetic_Age_Zhang", "DNAmPhenoAge"
)

# Rename groups for more user-friendly labels
labels_grupos <- c(
  "IV" = "Normal weight",
  "I" = "Sarcopenic obesity",
  "II" = "Obesity",
  "III" = "Sarcopenia"
)

# Rename variables for more readable titles
labels_variaveis <- c(
  "EEAA" = "EEAA",
  "IEAA_Hannum" = "IEAA",
  "DNAmTL" = "DNAmTL",
  "DNAmFitAge" = "DNAm FitAge",
  "DNAmAgeHannum" = "Hannum clock",
  "DNAmAgeHorvath" = "DNAm Age Horvath",
  "DNAmGrimAge" = "DNAm GrimAge",
  "DNAmGrimAge2" = "DNAm GrimAge 2",
  "Epigenetic_Age_Zhang" = "DNAm Zhang",
  "DNAmPhenoAge" = "DNAm PhenoAge"
)

# Order variables in the desired sequence
ordem_variaveis <- names(labels_variaveis)

# Prepare data for group analysis (long format)
dados_gp <- dados %>%
  select(gp, all_of(variaveis)) %>%
  pivot_longer(cols = -gp, names_to = "Variable", values_to = "Value") # Translated column names

# Prepare data for overall analysis (long format)
dados_geral <- dados %>%
  select(gp, all_of(variaveis)) %>%
  pivot_longer(cols = -gp, names_to = "Variable", values_to = "Value") # Translated column names

dados_completo <- bind_rows(dados_gp, dados_geral)

# Apply ordered factor for variables and rename with letter labels
letras <- LETTERS[1:length(ordem_variaveis)]
labels_com_letra <- paste0(letras, ") ", labels_variaveis[ordem_variaveis])

dados_completo$Variavel <- factor(
  dados_completo$Variavel,
  levels = ordem_variaveis,
  labels = labels_com_letra
)

# Apply factor for group with user-friendly labels
dados_completo$gp <- factor(
  dados_completo$gp,
  levels = c(names(labels_grupos)),
  labels = labels_grupos
)

# Colors for each group
cores_grupos <- c(
  "Normal weight" = "green",
  "Sarcopenic obesity" = "red",
  "Obesity" = "deepskyblue",
  "Sarcopenia" = "gold"
)

# Create plot
p <- ggplot(dados_completo, aes(x = gp, y = Valor, fill = gp)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.1, size = 0.8, alpha = 0.3) +
  facet_wrap(~Variavel, scales = "free", ncol = 2) +
  scale_fill_manual(values = cores_grupos) +
  theme_minimal() +
  labs(
    title = "Distribution of Epigenetic, Metabolic, and Chronological Age Measures",
    x = "",
    y = "Observed value", # Translated y-axis label
    fill = "Groups"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    strip.text = element_text(size = 10, face = "bold", hjust = 0)
  )

# Show the plot
print(p)

# Save as PNG (900 dpi)
ggsave(
  filename = "violin_plot.png", # Translated filename
  plot = p,
  device = "png",
  dpi = 900,
  width = 21, height = 29.7, units = "cm"
)

# Save as vector PDF
ggsave(
  filename = "violin_plot_for_publication.pdf", # Translated filename
  plot = p,
  device = "pdf",
  width = 21, height = 29.7, units = "cm"
)