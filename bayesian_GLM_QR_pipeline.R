# ----------------------------------------------------------------------#
# Bayesian Analysis of Epigenetic Aging in Sarcopenic Obesity         #
# ----------------------------------------------------------------------#

# ---------------------------------#
# 1. IMPORT DATA AND LOAD PACKAGES #
# ---------------------------------#
# Ensure cmdstanr is installed and configured for brms
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cmdstanr::install_cmdstan()
cmdstanr::set_cmdstan_path(
  file.path(path.expand("~"), ".cmdstan", "cmdstan-2.36.0"))
cmdstanr::cmdstan_version()

# Load required packages
library(brms)
library(ISLR)
library(sjPlot)
library(quantreg)
library(bayesQR)
library(loo)
library(ggplot2)
library(dplyr)
library(effectsize)
library(patchwork)
library(tidyr)
library(emmeans)
library(gtsummary)
library(flextable)
library(cmdstanr)

# Import the dataset
dados <- read.csv2(
"https://raw.githubusercontent.com/apolveiro01/Bayesian-Analysis-of-Epigenetic-Aging-in-Sarcopenic-Obesity/main/dados_epi.csv", 
header = TRUE)

# --------------------------#
# 2. DATA PREPROCESSING     #
# --------------------------#

# Convert variables to factors
dados$ID <- as.factor(dados$ID)
dados$gp_i <- as.factor(dados$gp_i)
dados$gp_n <- as.factor(dados$gp_n)
dados$gp <- factor(dados$gp, levels = c("IV", "I", "II", "III"))
# IV = Normal weight; E = Sarcopenic obesity; II = Obesity; E = Sarcopenia

# --------------------------------------------------------------------------
# 3. MODELING
# --------------------------------------------------------------------------

# Define priors for Bayesian linear models
priors_glm <- c(
  set_prior("normal(0, 10)", class = "b"),# priors suaves para o coeficiente
  set_prior("normal(0, 10)", class = "Intercept"), # para o intercepto
  set_prior("student_t(3, 0, 5)", class = "sigma")  # para a variabilidade
)

# Bayesian Linear Models

modelo_glm_bayes_DNAmFitAge <- brm(
  formula = DNAmFitAge ~ gp,
  family = gaussian(),
  data = dados,
  prior = priors_glm,
  chains = 4, iter = 4000, warmup = 1000, seed = 123,
  backend = "cmdstanr"
)

modelo_glm_bayes_DNAmFitAge_aj <- brm(
  formula = DNAmFitAge ~ gp + Age + HOMA_IR,
  family = gaussian(),
  data = dados,
  prior = priors_glm,
  chains = 4, iter = 4000, warmup = 1000, seed = 123,
  backend = "cmdstanr"
)

modelo_glm_bayes_EEAA <- brm(
  formula = EEAA ~ gp,
  family = gaussian(),
  data = dados,
  prior = priors_glm,
  chains = 4, iter = 4000, warmup = 1000, seed = 123,
  backend = "cmdstanr"
)

modelo_glm_bayes_EEAA_aj <- brm(
  formula = EEAA ~ gp + Age + HOMA_IR,
  family = gaussian(),
  data = dados,
  prior = priors_glm,
  chains = 4, iter = 4000, warmup = 1000, seed = 123,
  backend = "cmdstanr"
)

modelo_glm_bayes_DNAmTL <- brm(
  formula = DNAmTL ~ gp,
  family = gaussian(),
  data = dados,
  prior = priors_glm,
  chains = 4, iter = 4000, warmup = 1000, seed = 123,
  backend = "cmdstanr"
)

modelo_glm_bayes_DNAmTL_aj <- brm(
  formula = DNAmTL ~ gp + Age + HOMA_IR,
  family = gaussian(),
  data = dados,
  prior = priors_glm,
  chains = 4, iter = 4000, warmup = 1000, seed = 123,
  backend = "cmdstanr"
)

modelo_glm_bayes_DNAmAgeHannum <- brm(
  formula = DNAmAgeHannum ~ gp,
  family = gaussian(),
  data = dados,
  prior = priors_glm,
  chains = 4, iter = 4000, warmup = 1000, seed = 123,
  backend = "cmdstanr"
)

modelo_glm_bayes_DNAmAgeHannum_aj <- brm(
  formula = DNAmAgeHannum ~ gp + Age + HOMA_IR,
  family = gaussian(),
  data = dados,
  prior = priors_glm,
  chains = 4, iter = 4000, warmup = 1000, seed = 123,
  backend = "cmdstanr"
)

modelo_glm_bayes_DNAmAgeHorvath <- brm(
  formula = DNAmAgeHorvath ~ gp,
  family = gaussian(),
  data = dados,
  prior = priors_glm,
  chains = 4, iter = 4000, warmup = 1000, seed = 123,
  backend = "cmdstanr"
)

modelo_glm_bayes_DNAmAgeHorvath_aj <- brm(
  formula = DNAmAgeHorvath ~ gp + Age + HOMA_IR,
  family = gaussian(),
  data = dados,
  prior = priors_glm,
  chains = 4, iter = 4000, warmup = 1000, seed = 123,
  backend = "cmdstanr"
)

modelo_glm_bayes_DNAmPhenoAge <- brm(
  formula = DNAmPhenoAge ~ gp,
  family = gaussian(),
  data = dados,
  prior = priors_glm,
  chains = 4, iter = 4000, warmup = 1000, seed = 123,
  backend = "cmdstanr"
)

modelo_glm_bayes_DNAmPhenoAge_aj <- brm(
  formula = DNAmPhenoAge ~ gp + Age + HOMA_IR,
  family = gaussian(),
  data = dados,
  prior = priors_glm,
  chains = 4, iter = 4000, warmup = 1000, seed = 123,
  backend = "cmdstanr"
)

modelo_glm_bayes_DNAmGrimAge <- brm(
  formula = DNAmGrimAge ~ gp,
  family = gaussian(),
  data = dados,
  prior = priors_glm,
  chains = 4, iter = 4000, warmup = 1000, seed = 123,
  backend = "cmdstanr"
)

modelo_glm_bayes_DNAmGrimAge_aj <- brm(
  formula = DNAmGrimAge ~ gp + Age + HOMA_IR,
  family = gaussian(),
  data = dados,
  prior = priors_glm,
  chains = 4, iter = 4000, warmup = 1000, seed = 123,
  backend = "cmdstanr"
)

modelo_glm_bayes_DNAmGrimAge2 <- brm(
  formula = DNAmGrimAge2~ gp,
  family = gaussian(),
  data = dados,
  prior = priors_glm,
  chains = 4, iter = 4000, warmup = 1000, seed = 123,
  backend = "cmdstanr"
)

modelo_glm_bayes_DNAmGrimAge2_aj <- brm(
  formula = DNAmGrimAge2 ~ gp + Age + HOMA_IR,
  family = gaussian(),
  data = dados,
  prior = priors_glm,
  chains = 4, iter = 4000, warmup = 1000, seed = 123,
  backend = "cmdstanr"
)

modelo_glm_bayes_DNAmGrimAge2 <- brm(
  formula = DNAmGrimAge2~ gp + Age + HOMA_IR,
  family = gaussian(),
  data = dados,
  prior = priors_glm,
  chains = 4, iter = 4000, warmup = 1000, seed = 123,
  backend = "cmdstanr"
)

modelo_glm_bayes_DNAmGrimAge2_aj <- brm(
  formula = DNAmGrimAge2 ~ gp + Age + HOMA_IR,
  family = gaussian(),
  data = dados,
  prior = priors_glm,
  chains = 4, iter = 4000, warmup = 1000, seed = 123,
  backend = "cmdstanr"
)

modelo_glm_bayes_Epigenetic_Age_Zhang <- brm(
  formula = Epigenetic_Age_Zhang ~ gp,
  family = gaussian(),
  data = dados,
  prior = priors_glm,
  chains = 4, iter = 4000, warmup = 1000, seed = 123,
  backend = "cmdstanr"
)

modelo_glm_bayes_Epigenetic_Age_Zhang_aj <- brm(
  formula = Epigenetic_Age_Zhang ~ gp + Age + HOMA_IR,
  family = gaussian(),
  data = dados,
  prior = priors_glm,
  chains = 4, iter = 4000, warmup = 1000, seed = 123,
  backend = "cmdstanr"
)

modelo_glm_bayes_IEAA_Hannum <- brm(
  formula = IEAA_Hannum ~ gp,
  family = gaussian(),
  data = dados,
  prior = priors_glm,
  chains = 4, iter = 4000, warmup = 1000, seed = 123,
  backend = "cmdstanr"
)

modelo_glm_bayes_IEAA_Hannum_aj <- brm(
  formula = IEAA_Hannum ~ gp + Age + HOMA_IR,
  family = gaussian(),
  data = dados,
  prior = priors_glm,
  chains = 4, iter = 4000, warmup = 1000, seed = 123,
  backend = "cmdstanr"
)

# Define priors for quantile regression models
priors_qr <- function(taus) {
  if (taus == 0.5) {
    return(c(
      set_prior("normal(0, 10)", class = "b"),
      set_prior("normal(0, 10)", class = "Intercept"),
      set_prior("student_t(3, 0, 5)", class = "sigma")
    ))
  } else if (taus %in% c(0.2, 0.8)) {
    return(c(
      set_prior("normal(0, 3)", class = "b"),
      set_prior("normal(0, 3)", class = "Intercept"),
      set_prior("student_t(3, 0, 3)", class = "sigma")
    ))
  } else {
    return(c(
      set_prior("normal(0, 5)", class = "b"),
      set_prior("normal(0, 5)", class = "Intercept"),
      set_prior("student_t(3, 0, 4)", class = "sigma")
    ))
  }
}

# Define desired quantiles
taus <- seq(0.20, 0.80, by = 0.05)

# Bayesian quantile regression models

modelos_bqr_brms_EEAA <- lapply(taus, function(t) {
  brm(
    formula = bf(EEAA ~ gp, quantile = t),
    data = dados,
    family = asym_laplace(),
    prior = priors_qr(t),
    chains = 4, iter = 4000, warmup = 1000, seed = 1234,
    backend = "cmdstanr"
  )
})

# Name the models based on their respective quantile (tau) values
names(modelos_bqr_brms_EEAA) <- paste0("tau_", taus)
 

modelos_bqr_brms_EEAA_aj <- lapply(taus, function(t) {
  brm(
    formula = bf(EEAA ~ gp + Age + HOMA_IR, quantile = t),
    data = dados,
    family = asym_laplace(),
    prior = priors_qr(t),
    chains = 4, iter = 4000, warmup = 1000, seed = 1234,
    backend = "cmdstanr"
  )
})

# Name the models based on their respective quantile (tau) values
names(modelos_bqr_brms_EEAA_aj) <- paste0("tau_", taus)


modelos_bqr_brms_DNAmTL <- lapply(taus, function(t) {
  brm(
    formula = bf(DNAmTL ~ gp, quantile = t),
    data = dados,
    family = asym_laplace(),
    prior = priors_qr(t),
    chains = 4, iter = 4000, warmup = 1000, seed = 1234,
    backend = "cmdstanr"
  )
})

# Name the models based on their respective quantile (tau) values
names(modelos_bqr_brms_DNAmTL) <- paste0("tau_", taus)


modelos_bqr_brms_DNAmTL_aj <- lapply(taus, function(t) {
  brm(
    formula = bf(DNAmTL ~ gp + Age + HOMA_IR, quantile = t),
    data = dados,
    family = asym_laplace(),
    prior = priors_qr(t),
    chains = 4, iter = 4000, warmup = 1000, seed = 1234,
    backend = "cmdstanr"
  )
})

# Name the models based on their respective quantile (tau) values
names(modelos_bqr_brms_DNAmTL_aj) <- paste0("tau_", taus)


modelos_bqr_brms_DNAmAgeHorvath <- lapply(taus, function(t) {
  brm(
    formula = bf(DNAmAgeHorvath ~ gp, quantile = t),
    data = dados,
    family = asym_laplace(),
    prior = priors_qr(t),
    chains = 4, iter = 4000, warmup = 1000, seed = 1234,
    backend = "cmdstanr"
  )
})

# Name the models based on their respective quantile (tau) values
names(modelos_bqr_brms_DNAmAgeHorvath) <- paste0("tau_", taus)


modelos_bqr_brms_DNAmAgeHorvath_aj <- lapply(taus, function(t) {
  brm(
    formula = bf(DNAmAgeHorvath ~ gp + Age + HOMA_IR, quantile = t),
    data = dados,
    family = asym_laplace(),
    prior = priors_qr(t),
    chains = 4, iter = 4000, warmup = 1000, seed = 1234,
    backend = "cmdstanr"
  )
})

# Name the models based on their respective quantile (tau) values
names(modelos_bqr_brms_DNAmAgeHorvath_aj) <- paste0("tau_", taus)


modelos_bqr_brms_DNAmAgeHannum <- lapply(taus, function(t) {
  brm(
    formula = bf(DNAmAgeHannum ~ gp, quantile = t),
    data = dados,
    family = asym_laplace(),
    prior = priors_qr(t),
    chains = 4, iter = 4000, warmup = 1000, seed = 1234,
    backend = "cmdstanr"
  )
})

# Name the models based on their respective quantile (tau) values
names(modelos_bqr_brms_DNAmAgeHannum) <- paste0("tau_", taus)


modelos_bqr_brms_DNAmAgeHannum_aj <- lapply(taus, function(t) {
  brm(
    formula = bf(DNAmAgeHannum ~ gp + Age + HOMA_IR, quantile = t),
    data = dados,
    family = asym_laplace(),
    prior = priors_qr(t),
    chains = 4, iter = 4000, warmup = 1000, seed = 1234,
    backend = "cmdstanr"
  )
})

# Name the models based on their respective quantile (tau) values
names(modelos_bqr_brms_DNAmAgeHannum_aj) <- paste0("tau_", taus)


modelos_bqr_brms_DNAmPhenoAge <- lapply(taus, function(t) {
  brm(
    formula = bf(DNAmPhenoAge ~ gp, quantile = t),
    data = dados,
    family = asym_laplace(),
    prior = priors_qr(t),
    chains = 4, iter = 4000, warmup = 1000, seed = 1234,
    backend = "cmdstanr"
  )
})

# Name the models based on their respective quantile (tau) values
names(modelos_bqr_brms_DNAmPhenoAge) <- paste0("tau_", taus)


modelos_bqr_brms_DNAmPhenoAge_aj <- lapply(taus, function(t) {
  brm(
    formula = bf(DNAmPhenoAge ~ gp + Age + HOMA_IR, quantile = t),
    data = dados,
    family = asym_laplace(),
    prior = priors_qr(t),
    chains = 4, iter = 4000, warmup = 1000, seed = 1234,
    backend = "cmdstanr"
  )
})

# Name the models based on their respective quantile (tau) values
names(modelos_bqr_brms_DNAmPhenoAge_aj) <- paste0("tau_", taus)


modelos_bqr_brms_IEAA_Hannum <- lapply(taus, function(t) {
  brm(
    formula = bf(IEAA_Hannum ~ gp, quantile = t),
    data = dados,
    family = asym_laplace(),
    prior = priors_qr(t),
    chains = 4, iter = 4000, warmup = 1000, seed = 1234,
    backend = "cmdstanr"
  )
})

# Name the models based on their respective quantile (tau) values
names(modelos_bqr_brms_IEAA_Hannum) <- paste0("tau_", taus)


modelos_bqr_brms_IEAA_Hannum_aj <- lapply(taus, function(t) {
  brm(
    formula = bf(IEAA_Hannum ~ gp + Age + HOMA_IR, quantile = t),
    data = dados,
    family = asym_laplace(),
    prior = priors_qr(t),
    chains = 4, iter = 4000, warmup = 1000, seed = 1234,
    backend = "cmdstanr"
  )
})

# Name the models based on their respective quantile (tau) values
names(modelos_bqr_brms_IEAA_Hannum_aj) <- paste0("tau_", taus)


modelos_bqr_brms_Epigenetic_Age_Zhang <- lapply(taus, function(t) {
  brm(
    formula = bf(Epigenetic_Age_Zhang ~ gp, quantile = t),
    data = dados,
    family = asym_laplace(),
    prior = priors_qr(t),
    chains = 4, iter = 4000, warmup = 1000, seed = 1234,
    backend = "cmdstanr"
  )
})

# Name the models based on their respective quantile (tau) values
names(modelos_bqr_brms_Epigenetic_Age_Zhang) <- paste0("tau_", taus)


modelos_bqr_brms_Epigenetic_Age_Zhang_aj <- lapply(taus, function(t) {
  brm(
    formula = bf(Epigenetic_Age_Zhang ~ gp + Age + HOMA_IR, quantile = t),
    data = dados,
    family = asym_laplace(),
    prior = priors_qr(t),
    chains = 4, iter = 4000, warmup = 1000, seed = 1234,
    backend = "cmdstanr"
  )
})

# Name the models based on their respective quantile (tau) values
names(modelos_bqr_brms_Epigenetic_Age_Zhang_aj) <- paste0("tau_", taus)


modelos_bqr_brms_DNAmFitAge <- lapply(taus, function(t) {
  brm(
    formula = bf(DNAmFitAge ~ gp, quantile = t),
    data = dados,
    family = asym_laplace(),
    prior = priors_qr(t),
    chains = 4, iter = 4000, warmup = 1000, seed = 1234,
    backend = "cmdstanr"
  )
})

# Name the models based on their respective quantile (tau) values
names(modelos_bqr_brms_DNAmFitAge) <- paste0("tau_", taus)


modelos_bqr_brms_DNAmFitAge_aj <- lapply(taus, function(t) {
  brm(
    formula = bf(DNAmFitAge ~ gp + Age + HOMA_IR, quantile = t),
    data = dados,
    family = asym_laplace(),
    prior = priors_qr(t),
    chains = 4, iter = 4000, warmup = 1000, seed = 1234,
    backend = "cmdstanr"
  )
})

# Name the models based on their respective quantile (tau) values
names(modelos_bqr_brms_DNAmFitAge_aj) <- paste0("tau_", taus)


modelos_bqr_brms_DNAmGrimAge2<- lapply(taus, function(t) {
  brm(
    formula = bf(DNAmGrimAge2 ~ gp, quantile = t),
    data = dados,
    family = asym_laplace(),
    prior = priors_qr(t),
    chains = 4, iter = 4000, warmup = 1000, seed = 1234,
    backend = "cmdstanr"
  )
})

# Name the models based on their respective quantile (tau) values
names(modelos_bqr_brms_DNAmGrimAge2) <- paste0("tau_", taus)


modelos_bqr_brms_DNAmGrimAge2_aj <- lapply(taus, function(t) {
  brm(
    formula = bf(EEAA ~ gp + Age + HOMA_IR, quantile = t),
    data = dados,
    family = asym_laplace(),
    prior = priors_qr(t),
    chains = 4, iter = 4000, warmup = 1000, seed = 1234,
    backend = "cmdstanr"
  )
})

# Name the models based on their respective quantile (tau) values
names(modelos_bqr_brms_DNAmGrimAge2_aj) <- paste0("tau_", taus)


modelos_bqr_brms_DNAmGrimAge <- lapply(taus, function(t) {
  brm(
    formula = bf(DNAmGrimAge ~ gp, quantile = t),
    data = dados,
    family = asym_laplace(),
    prior = priors_qr(t),
    chains = 4, iter = 4000, warmup = 1000, seed = 1234,
    backend = "cmdstanr"
  )
})

# Name the models based on their respective quantile (tau) values
names(modelos_bqr_brms_DNAmGrimAge) <- paste0("tau_", taus)


modelos_bqr_brms_DNAmGrimAge_aj <- lapply(taus, function(t) {
  brm(
    formula = bf(DNAmGrimAge ~ gp + Age + HOMA_IR, quantile = t),
    data = dados,
    family = asym_laplace(),
    prior = priors_qr(t),
    chains = 4, iter = 4000, warmup = 1000, seed = 1234,
    backend = "cmdstanr"
  )
})

# Name the models based on their respective quantile (tau) values
names(modelos_bqr_brms_DNAmGrimAge_aj) <- paste0("tau_", taus)

#=====================================#
# LIST OF MODELS - BAYESIAN GLMs      #
#=====================================#
#modelo_glm_bayes_EEAA
#modelo_glm_bayes_EEAA_aj
#modelo_glm_bayes_DNAmTL
#modelo_glm_bayes_DNAmTL_aj
#modelo_glm_bayes_DNAmAgeHorvath
#modelo_glm_bayes_DNAmAgeHorvath_aj
#modelo_glm_bayes_DNAmAgeHannum
#modelo_glm_bayes_DNAmAgeHannum_aj
#modelo_glm_bayes_DNAmPhenoAge
#modelo_glm_bayes_DNAmPhenoAge_aj
#modelo_glm_bayes_IEAA_Hannum
#modelo_glm_bayes_IEAA_Hannum_aj
#modelo_glm_bayes_Epigenetic_Age_Zhang
#modelo_glm_bayes_Epigenetic_Age_Zhang_aj
#modelo_glm_bayes_DNAmFitAge
#modelo_glm_bayes_DNAmFitAge_aj
#modelo_glm_bayes_DNAmGrimAge2
#modelo_glm_bayes_DNAmGrimAge2_aj
#modelo_glm_bayes_DNAmGrimAge
#modelo_glm_bayes_DNAmGrimAge_aj

#=====================================#
# LIST OF MODELS - BAYESIAN QR        #
#=====================================#
#modelos_bqr_brms_EEAA
#modelos_bqr_brms_EEAA_aj
#modelos_bqr_brms_DNAmTL
#modelos_bqr_brms_DNAmTL_aj
#modelos_bqr_brms_DNAmAgeHorvath
#modelos_bqr_brms_DNAmAgeHorvath_aj
#modelos_bqr_brms_DNAmAgeHannum
#modelos_bqr_brms_DNAmAgeHannum_aj
#modelos_bqr_brms_DNAmPhenoAge
#modelos_bqr_brms_DNAmPhenoAge_aj
#modelos_bqr_brms_IEAA_Hannum
#modelos_bqr_brms_IEAA_Hannum_aj
#modelos_bqr_brms_Epigenetic_Age_Zhang
#modelos_bqr_brms_Epigenetic_Age_Zhang_aj
#modelos_bqr_brms_DNAmFitAge
#modelos_bqr_brms_DNAmFitAge_aj
#modelos_bqr_brms_DNAmGrimAge2
#modelos_bqr_brms_DNAmGrimAge2_aj
#modelos_bqr_brms_DNAmGrimAge
#modelos_bqr_brms_DNAmGrimAge_aj

# ----------------------------#
# Models Diagnostics #
# ----------------------------#

# Pacotes necessários
library(brms)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(purrr)
library(stringr)
library(writexl)
library(posterior)

# -----------------------------------------#
# Function to export plots and Bayesian R² #
# -----------------------------------------#
run_diagnostics <- function(model, model_name = "Model", output_dir = "diagnostics/") { 
  cat("\n====== Diagnostics for:", model_name, "======\n") 
  
  model_dir <- file.path(output_dir, model_name)
  if (!dir.exists(model_dir)) dir.create(model_dir, recursive = TRUE)
  
  arr_model <- as.array(model)
  
  # Select valid parameter names containing "b_"
  all_pars <- dimnames(arr_model)$parameters
  b_pars <- grep("^b_", all_pars, value = TRUE)
  if ("b_Intercept" %in% all_pars) b_pars <- c("b_Intercept", b_pars)
  
  plots <- list(
    trace = mcmc_trace(arr_model, pars = b_pars) + ggtitle(paste("Trace plot -", model_name)),
    acf = mcmc_acf(arr_model, pars = b_pars) + ggtitle(paste("Autocorrelation -", model_name)),
    rank = mcmc_rank_overlay(arr_model, pars = b_pars) + ggtitle(paste("Rank plot -", model_name)),
    ppc = pp_check(model) + ggtitle(paste("PPC -", model_name)),
    density = mcmc_dens_overlay(arr_model, pars = b_pars) + 
      ggtitle(paste("Posterior Densities -", model_name))
  )
  
  # Residuals vs. Fitted Values Plot
  # Calculate predicted (fitted) values and residuals
  preds  <- fitted(model)[, "Estimate"]
  resid  <- residuals(model)[, "Estimate"]
  
  # Create a data frame for plotting
  df_resid <- data.frame(Fitted = preds, Residual = resid)
  
  # Generate the residuals vs. fitted plot
  resid_plot <- ggplot(df_resid, aes(x = Fitted, y = Residual)) + 
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    theme_minimal() +
    ggtitle(paste("Residuals vs. Fitted Values -", model_name))
  
  plots$resid <- resid_plot
  
  # Export Plots (PNG and PDF)
  for (gname in names(plots)) {
    g <- plots[[gname]]
    ggsave(filename = file.path(model_dir, paste0(gname, ".png")), plot = g, width = 8, height = 6, dpi = 300)
    ggsave(filename = file.path(model_dir, paste0(gname, ".pdf")), plot = g, width = 8, height = 6)
  }
  
  # Check Rhat > 1.1 (Convergence Diagnostic)
  rhats <- rhat(as_draws_array(model))
  max_rhat <- max(rhats, na.rm = TRUE)
  rhats_high <- rhats[rhats > 1.1]
  
  # Effective Number of Samples (N_eff Ratio)
  neff <- tryCatch({
    neff_ratio(model)
  }, error = function(e) {
    warning("Error calculating neff_ratio: ", conditionMessage(e))
    NULL
  })
  
  neff_low <- if (!is.null(neff)) neff[neff < 0.1] else NA
  
  
  # Bayesian R² (if available)
  r2_file <- file.path(model_dir, "bayes_R2.txt")
  r2_result <- tryCatch({
    r2 <- bayes_R2(model)
    r2_mean <- mean(r2)
    writeLines(paste("Bayesian R² (mean):", round(r2_mean, 4)), r2_file)
    cat("Bayesian R² (mean):", round(r2_mean, 4), "\n")
  }, error = function(e) {
    writeLines("Bayesian R² not available for this model.", r2_file)
    cat("⚠️ Bayesian R² not available for:", model_name, "\n")
  })
  
  
  # Numerical Summary of the Model
  sink(file.path(model_dir, "summary.txt"))
  print(summary(model))
  sink()
  
  cat("Diagnostics saved to:", model_dir, "\n")
}

# -------------------------------------------------------------
# Diagnostics for List of Bayesian Quantile Regression (BQR) Models
# -------------------------------------------------------------
run_diagnostics_bqr_seq <- function(models_list, tau_seq, base_name = "BQR", output_dir = "diagnostics/") {
  for (i in seq_along(tau_seq)) {
    model <- models_list[[i]]
    tau_name <- paste0(base_name, "_tau_", tau_seq[i])
    run_diagnostics(model, model_name = tau_name, output_dir = output_dir)
  }
}

#=======================#
# Function Calls        #
#=======================#

# Generalized Linear Models (GLM)
modelos_bglm <- c(
  "modelo_glm_bayes_EEAA", "modelo_glm_bayes_EEAA_aj", "modelo_glm_bayes_DNAmTL", "modelo_glm_bayes_DNAmTL_aj",
  "modelo_glm_bayes_DNAmAgeHorvath", "modelo_glm_bayes_DNAmAgeHorvath_aj", "modelo_glm_bayes_DNAmAgeHannum",
  "modelo_glm_bayes_DNAmAgeHannum_aj", "modelo_glm_bayes_DNAmPhenoAge", "modelo_glm_bayes_DNAmPhenoAge_aj", 
  "modelo_glm_bayes_IEAA_Hannum","modelo_glm_bayes_IEAA_Hannum_aj","modelo_glm_bayes_Epigenetic_Age_Zhang",
  "modelo_glm_bayes_Epigenetic_Age_Zhang_aj","modelo_glm_bayes_DNAmFitAge","modelo_glm_bayes_DNAmFitAge_aj", 
  "modelo_glm_bayes_DNAmGrimAge2", "modelo_glm_bayes_DNAmGrimAge2_aj","modelo_glm_bayes_DNAmGrimAge",
  "modelo_glm_bayes_DNAmGrimAge_aj"
)

cat("### START: Diagnostics for BGLM models ###\n")
for (mod_name in modelos_bglm) {
  cat("\n>> Running diagnostics for:", mod_name, "\n")
  model <- tryCatch(get(mod_name), error = function(e) NULL)
  if (!is.null(model)) run_diagnostics(model, model_name = mod_name)
}
cat("### END: Diagnostics for BGLM models ###\n\n")

# Quantile Regression (QR)
tau_seq <- seq(0.2, 0.8, by = 0.05)

modelos_bqr <- c(
  "modelos_bqr_brms_EEAA", "modelos_bqr_brms_EEAA_aj", "modelos_bqr_brms_DNAmTL", "modelos_bqr_brms_DNAmTL_aj",
  "modelos_bqr_brms_DNAmAgeHorvath", "modelos_bqr_brms_DNAmAgeHorvath_aj", "modelos_bqr_brms_DNAmAgeHannum",
  "modelos_bqr_brms_DNAmAgeHannum_aj", "modelos_bqr_brms_DNAmPhenoAge", "modelos_bqr_brms_DNAmPhenoAge_aj",
  "modelos_bqr_brms_IEAA_Hannum", "modelos_bqr_brms_IEAA_Hannum_aj",
  "modelos_bqr_brms_Epigenetic_Age_Zhang", "modelos_bqr_brms_Epigenetic_Age_Zhang_aj", "modelos_bqr_brms_DNAmFitAge",
  "modelos_bqr_brms_DNAmFitAge_aj", "modelos_bqr_brms_DNAmGrimAge2", "modelos_bqr_brms_DNAmGrimAge2_aj",
  "modelos_bqr_brms_DNAmGrimAge", "modelos_bqr_brms_DNAmGrimAge_aj"
)

cat("### START: Diagnostics for BQR models ###\n")
for (bqr_name in modelos_bqr) {
  cat("\n>> Running diagnostics for model set:", bqr_name, "\n")
  modelos_tau <- tryCatch(get(bqr_name), error = function(e) NULL)
  if (is.list(modelos_tau)) run_diagnostics_bqr_seq(modelos_tau, tau_seq, base_name = bqr_name)
}
cat("### END: Diagnostics for BQR models ###\n\n")

# ------------------------------------------------------------------#
# Diagnostics for List of Bayesian GLM Models                       #
# ------------------------------------------------------------------#
library(tidyverse)
library(brms)
library(loo)

# List of brms models (not a list of lists!)
todas_listas_glm <- list(
  EEAA = modelo_glm_bayes_EEAA,
  EEAA_aj = modelo_glm_bayes_EEAA_aj,
  DNAmTL = modelo_glm_bayes_DNAmTL,
  DNAmTL_aj = modelo_glm_bayes_DNAmTL_aj,
  DNAmAgeHorvath = modelo_glm_bayes_DNAmAgeHorvath,
  DNAmAgeHorvath_aj = modelo_glm_bayes_DNAmAgeHorvath_aj,
  DNAmAgeHannum = modelo_glm_bayes_DNAmAgeHannum,
  DNAmAgeHannum_aj = modelo_glm_bayes_DNAmAgeHannum_aj,
  DNAmPhenoAge = modelo_glm_bayes_DNAmPhenoAge,
  DNAmPhenoAge_aj = modelo_glm_bayes_DNAmPhenoAge_aj,
  IEAA = modelo_glm_bayes_IEAA_Hannum,
  IEAA_aj = modelo_glm_bayes_IEAA_Hannum_aj,
  DNAmAge_Zhang = modelo_glm_bayes_Epigenetic_Age_Zhang,
  DNAmAge_Zhang_aj = modelo_glm_bayes_Epigenetic_Age_Zhang_aj,
  DNAmFitAge = modelo_glm_bayes_DNAmFitAge,
  DNAmFitAge_aj = modelo_glm_bayes_DNAmFitAge_aj,
  DNAmGrimAge = modelo_glm_bayes_DNAmGrimAge,
  DNAmGrimAge_aj = modelo_glm_bayes_DNAmGrimAge_aj,
  DNAmGrimAge2 = modelo_glm_bayes_DNAmGrimAge2,
  DNAmGrimAge2_aj = modelo_glm_bayes_DNAmGrimAge2_aj
)

# Apply LOO (Leave-One-Out cross-validation) to each model
diagnostico_pareto <- imap_dfr(todas_listas_glm, function(modelo, nome_modelo) {
  loo_result <- loo(modelo, reloo = TRUE)
  pareto_k <- loo_result$diagnostics$pareto_k
  
  tibble(
    model = nome_modelo, 
    outcome = str_remove(nome_modelo, "_aj$"), 
    adjusted = ifelse(str_detect(nome_modelo, "_aj$"), "Yes", "No"),
    elpd_loo = loo_result$estimates["elpd_loo", "Estimate"],
    se_elpd_loo = loo_result$estimates["elpd_loo", "SE"],
    looic = loo_result$estimates["looic", "Estimate"],
    max_pareto_k = max(pareto_k),
    n_obs = length(pareto_k)
  )
})

# Export final table
write_xlsx(diagnostico_pareto, "diagnostics/pareto_k_diagnostics_glm_models.xlsx") 
cat("✅ Pareto k diagnostics exported to 'diagnostics/pareto_k_diagnostics_glm_models.xlsx'\n") #

# -------------------------------------------#
# Diagnostics for List of Bayesian QR Models #
# -------------------------------------------#

# Combine all model lists
todas_listas_QR <- list(
  EEAA = modelos_bqr_brms_EEAA,
  EEAA_aj = modelos_bqr_brms_EEAA_aj,
  DNAmTL = modelos_bqr_brms_DNAmTL,
  DNAmTL_aj = modelos_bqr_brms_DNAmTL_aj,
  DNAmAgeHorvath = modelos_bqr_brms_DNAmAgeHorvath,
  DNAmAgeHorvath_aj = modelos_bqr_brms_DNAmAgeHorvath_aj,
  DNAmAgeHannum = modelos_bqr_brms_DNAmAgeHannum,
  DNAmAgeHannum_aj = modelos_bqr_brms_DNAmAgeHannum_aj,
  DNAmPhenoAge = modelos_bqr_brms_DNAmPhenoAge,
  DNAmPhenoAge_aj = modelos_bqr_brms_DNAmPhenoAge_aj,
  IEAA_Hannum = modelos_bqr_brms_IEAA_Hannum,
  IEAA_Hannum_aj = modelos_bqr_brms_IEAA_Hannum_aj,
  Epigenetic_Age_Zhang = modelos_bqr_brms_Epigenetic_Age_Zhang,
  Epigenetic_Age_Zhang_aj = modelos_bqr_brms_Epigenetic_Age_Zhang_aj,
  DNAmFitAge = modelos_bqr_brms_DNAmFitAge,
  DNAmFitAge_aj = modelos_bqr_brms_DNAmFitAge_aj,
  DNAmGrimAge2 = modelos_bqr_brms_DNAmGrimAge2,
  DNAmGrimAge2_aj = modelos_bqr_brms_DNAmGrimAge2_aj,
  DNAmGrimAge = modelos_bqr_brms_DNAmGrimAge,
  DNAmGrimAge_aj = modelos_bqr_brms_DNAmGrimAge_aj
)

# Flatten the list of lists into a single list of models with metadata
modelos_planos_QR <- imap(todas_listas_QR, function(lista_modelos, nome_lista) {
  imap(lista_modelos, function(modelo, nome_modelo) {
    full_name <- paste0(nome_lista, "_", nome_modelo) # Renamed 'nome_completo' to 'full_name'
    list(
      model = modelo, # Renamed 'modelo' to 'model'
      full_name = full_name,
      outcome = str_remove(nome_lista, "_aj$"), # Renamed 'desfecho' to 'outcome'
      adjusted = ifelse(str_detect(nome_lista, "_aj$"), "Yes", "No"), # Renamed 'ajustado' to 'adjusted'
      quantile = as.numeric(str_extract(nome_modelo, "\\d+")) / 100 # Renamed 'quantil' to 'quantile'
    )
  })
}) %>% flatten()

# Apply LOO (Leave-One-Out cross-validation)
diagnostico_pareto <- map_dfr(modelos_planos_QR, function(m) {
  loo_result <- loo(m$model, reloo = TRUE) # Used m$model
  pareto_k <- loo_result$diagnostics$pareto_k
  tibble(
    model = m$full_name, # Consistent naming
    outcome = m$outcome,
    adjusted = m$adjusted,
    quantile = m$quantile,
    elpd_loo = loo_result$estimates["elpd_loo", "Estimate"],
    se_elpd_loo = loo_result$estimates["elpd_loo", "SE"],
    looic = loo_result$estimates["looic", "Estimate"],
    max_pareto_k = max(pareto_k),
    n_obs = length(pareto_k)
  )
})

# Export final table
write_xlsx(diagnostico_pareto, "diagnostics/pareto_k_diagnostics_bqr_models.xlsx") 
cat("✅ Pareto k diagnostics exported to 'diagnostics/pareto_k_diagnostics_bqr_models.xlsx'\n") #


# Measure execution time
start_time <- Sys.time()
# ... model code ... 
end_time <- Sys.time()
print(end_time - start_time)

# ---------------------------------------------------
# 3. GRAPHICAL VISUALIZATION OF MODELS
# ---------------------------------------------------
# Required packages
library(tidyverse)
library(patchwork)

# Extract coefficients from quantile regression models
coefs_bqr1 <- map_dfr(1:length(taus), function(i) {
  fixef(modelos_bqr_brms_Epigenetic_Age_Zhang[[i]]) %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    filter(term %in% c("gpI", "gpII", "gpIII")) %>% 
    mutate(
      tau = taus[i],
      term = recode(term,
                    "gpI"   = "Sarcopenic obesity",
                    "gpII"  = "Obesity",
                    "gpIII" = "Sarcopenia"),
      term = factor(term, levels = c("Sarcopenic obesity", "Obesity", "Sarcopenia"))
    )
})

# Extract coefficients from quantile regression models (Adjusted)
coefs_bqr2 <- map_dfr(1:length(taus), function(i) {
  fixef(modelos_bqr_brms_Epigenetic_Age_Zhang_aj[[i]]) %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    filter(term %in% c("gpI", "gpII", "gpIII")) %>% 
    mutate(
      tau = taus[i],
      term = recode(term,
                    "gpI"   = "Sarcopenic obesity",
                    "gpII"  = "Obesity",
                    "gpIII" = "Sarcopenia"),
      term = factor(term, levels = c("Sarcopenic obesity", "Obesity", "Sarcopenia"))
    )
})

# Extract coefficients from quantile regression models
coefs_bqr3 <- map_dfr(1:length(taus), function(i) {
  fixef(modelos_bqr_brms_DNAmPhenoAge[[i]]) %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    filter(term != "Intercept") %>%
    mutate(
      tau = taus[i]
    ) %>%
    mutate(
      term = recode(term,
                    "gpI"   = "Sarcopenic obesity",
                    "gpII"  = "Obesity",
                    "gpIII" = "Sarcopenia",
                    "gpIV"  = "Normal weight"),
      term = factor(term, levels = c("Sarcopenic obesity", "Obesity", "Sarcopenia", "Normal weight"))
    )
})

# Extract coefficients from quantile regression models (Adjusted)
coefs_bqr4 <- map_dfr(1:length(taus), function(i) {
  fixef(modelos_bqr_brms_DNAmPhenoAge_aj[[i]]) %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    filter(term %in% c("gpI", "gpII", "gpIII")) %>% 
    mutate(
      tau = taus[i],
      term = recode(term,
                    "gpI"   = "Sarcopenic obesity",
                    "gpII"  = "Obesity",
                    "gpIII" = "Sarcopenia"),
      term = factor(term, levels = c("Sarcopenic obesity", "Obesity", "Sarcopenia"))
    )
})

# ---------------------------------------------------
# Create 4 individual plots
# ---------------------------------------------------

g1 <- ggplot(coefs_bqr1, aes(x = tau, y = Estimate, color = term)) +
  geom_point(size = 3) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, fill = term), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2) +
  labs(x = "Quantile (τ)", y = "Estimates DNAm Zhang") +
  theme_minimal() +
  scale_color_manual(
    name = "Groups",  # Legend title for color
    values = c("Sarcopenic obesity" = "red", "Obesity" = "deepskyblue", "Sarcopenia" = "gold", "Normal weight" = "green"),
    guide = "none"
  ) +
  scale_fill_manual(
    name = "Groups",  # Legend title for fill
    values = c("Sarcopenic obesity" = "red", "Obesity" = "deepskyblue", "Sarcopenia" = "gold", "Normal weight" = "green"),
    guide = "none"  
  )

g2 <- ggplot(coefs_bqr2, aes(x = tau, y = Estimate, color = term)) +
  geom_point(size = 3) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, fill = term), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2) +
  labs(x = "Quantile (τ)", y = "Estimates DNAm Zhang model 2") +
  theme_minimal() +
  scale_color_manual(
    name = "Groups",  # Legend title for color
    values = c("Sarcopenic obesity" = "red", "Obesity" = "deepskyblue", "Sarcopenia" = "gold", "Normal weight" = "green")
  ) +
  scale_fill_manual(
    name = "Groups",  # Legend title for fill
    values = c("Sarcopenic obesity" = "red", "Obesity" = "deepskyblue", "Sarcopenia" = "gold", "Normal weight" = "green")
  )

g3 <- ggplot(coefs_bqr3, aes(x = tau, y = Estimate, color = term)) +
  geom_point(size = 3) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, fill = term), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2) +
  labs(x = "Quantile (τ)", y = "Estimates DNAm PhenoAge") +
  theme_minimal() +
  scale_color_manual(
    name = "Groups",  # Legend title for color
    values = c("Sarcopenic obesity" = "red", "Obesity" = "deepskyblue", "Sarcopenia" = "gold", "Normal weight" = "green"),
    guide = "none"
  ) +
  scale_fill_manual(
    name = "Groups",  # Legend title for fill
    values = c("Sarcopenic obesity" = "red", "Obesity" = "deepskyblue", "Sarcopenia" = "gold", "Normal weight" = "green"),
    guide = "none"
  )

g4 <- ggplot(coefs_bqr4, aes(x = tau, y = Estimate, color = term)) +
  geom_point(size = 3) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, fill = term), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2) +
  labs(x = "Quantile (τ)", y = "Estimates DNAm PhenoAge model 2") +
  theme_minimal() +
  scale_color_manual(
    name = "Groups",  # Legend title for color
    values = c("Sarcopenic obesity" = "red", "Obesity" = "deepskyblue", "Sarcopenia" = "gold", "Normal weight" = "green"),
    guide = "none"
  ) +
  scale_fill_manual(
    name = "Groups",  # Legend title for fill
    values = c("Sarcopenic obesity" = "red", "Obesity" = "deepskyblue", "Sarcopenia" = "gold", "Normal weight" = "green"),
    guide = "none"
  )
# ---------------------------------------------------
# Combine the 4 plots with tags A), B), C), D)...
# ---------------------------------------------------
painel <- (g1 + g2) / (g3 + g4) +
  plot_annotation(tag_levels = list(LETTERS[1:4]))

library(ggplot2)
library(patchwork) 

# Assuming your unified plot is stored in 'painel'
painel <- (g1 + g2) / (g3 + g4) +
  plot_annotation(tag_levels = list(LETTERS[1:4]))

# Save image as PNG with high resolution (900 dpi)
ggsave(
  filename = "unified_plot_models_names_900dpi.png",
  plot = painel,
  width = 14,      # width in inches (adjust as desired)
  height = 10,     # height in inches
  dpi = 900,       # resolution in dpi
  units = "in"
)

ggsave(
  filename = "unified_plot_models_names_900dpi.tiff", # Updated filename
  plot = painel,
  width = 14,
  height = 10,
  dpi = 900,
  units = "in",
  compression = "lzw"
)

# ---------------------------------------------------
# 3. TABLE VISUALIZATION
# ---------------------------------------------------
library(dplyr)
library(purrr)
library(tibble)
library(broom) 
library(tidyr)
library(flextable) 
library(officer)   
library(writexl)

# 1. Extract coefficients from GLM model
coef_glm <- fixef(modelo_glm_bayes_DNAmPhenoAge) %>%
  as.data.frame() %>%
  rownames_to_column("Parameter") %>% 
  filter(grepl("^gp", Parameter)) %>% 
  mutate(
    Model = "GLM (Bayes)",   
    Quantile = "Posterior mean",        
    Est_CI = sprintf("%.2f [%.2f; %.2f]", Estimate, Q2.5, Q97.5)
  )

# 2. Extract coefficients from quantile regression models
coef_bqr <- map2_df(modelos_bqr_brms_DNAmPhenoAge, taus, function(mod, tau) {
  fixef(mod) %>%
    as.data.frame() %>%
    rownames_to_column("Parameter") %>% 
    filter(grepl("^gp", Parameter)) %>%  
    mutate(
      Model = "BQR (Bayes)",
      Quantile = paste0("τ = ", tau),
      Est_CI = sprintf("%.2f [%.2f; %.2f]", Estimate, Q2.5, Q97.5)
    )
})

# 3. Unify the two datasets
tabela_final <- bind_rows(coef_glm, coef_bqr) %>%
  select(Model, Quantile, Parameter, Est_CI) %>% 
  pivot_wider(names_from = Parameter, values_from = Est_CI)

# Visualize in flextable
# Create formatted flextable
ft <- flextable(tabela_final) |>
  theme_booktabs() |> 
  fontsize(size = 10, part = "all") |> 
  align(align = "center", part = "header") |> 
  align(align = "left", part = "body") |> 
  bold(part = "header") |> 
  color(color = "black", part = "all") |> 
  set_header_labels(
    Model = "",
    Quantile = "Model",
    gpIV = "Normal weight",
    gpI = "Sarcopenic obesity",
    gpII = "Obesity",
    gpIII = "Sarcopenia"
    # Add as per your actual parameter names, such as Age and HOMA-IR, but don't 
    # forget to change coefficient extraction as well
  ) |>
  autofit()


# Export to Word
doc <- read_docx() |>
  body_add_flextable(value = ft)

print(doc, target = "table_models_DNAmPhenoAge.docx") 


# Export to Excel (without flextable formatting)
write_xlsx(tabela_final, "table_models_HOMA_IR.xlsx") 