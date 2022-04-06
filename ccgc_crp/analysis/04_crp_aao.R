# MR analysis CRP-AAO --------------------------------------------------
library(tidyverse)
library(MendelianRandomization)
library(TwoSampleMR)

# Prepare the data -----------------------------------------------------
load(here::here("ccgc_crp", "data", "mr_data_aao.rda"))
source(here::here("ccgc_crp", "analysis", "mr_misc_functions.R"))
set.seed(1456)

# R2<0.001 and R2<0.01 have the same instruments -----------------------

# 1) R2<0.001 ----------------------------------------------------------
# Forest plots of Wald estimates ---------------------------------------
dd <- mr_data_r2_0001$mr_dataset

inp <- MendelianRandomization::mr_input(
  bx = dd$beta_x,
  bxse = dd$se_x,
  by = dd$beta_y,
  byse = dd$se_y,
  exposure = "log CRP concentration",
  outcome = "AAO",
  snps = dd$SNP
)

fp <- MendelianRandomization::mr_forest(
  object = inp,
  snp_estimates = TRUE,
  methods = "ivw",
  ordered = TRUE
)

res <- MendelianRandomization::mr_ivw(inp)

# Compute summary statistics of the result -----------------------------
est <- round(res@Estimate, 3L)
lower <- round(res@CILower, 3L)
upper <- round(res@CIUpper, 3L)

results <- tibble(
  SNP = dd$SNP,
  Gene = dd$Closest_gene,
  `R^2` = dd$beta_x^2 * dd$freq_x * (1 - dd$freq_x),
  `F-statistic` = dd$beta_x^2/dd$se_x^2,
  est = est,
  lower = lower,
  upper = upper,
  pval = round(res@Pvalue, 3L)
)

primary_aao_r2_0001 <- list(
  "data" = dd, "mr_input" = inp, "mr_results" = results
)

# 2) R2<0.1 ------------------------------------------------------------
# Forest plots of Wald estimates ---------------------------------------
dd <- mr_data_r2_01$mr_dataset

inp <- MendelianRandomization::mr_input(
  bx = dd$beta_x,
  bxse = dd$se_x,
  by = dd$beta_y,
  byse = dd$se_y,
  exposure = "log CRP concentration",
  outcome = "AAO",
  snps = dd$SNP
)

inp_corr <- MendelianRandomization::mr_input(
  bx = dd$beta_x,
  bxse = dd$se_x,
  by = dd$beta_y,
  byse = dd$se_y,
  exposure = "log CRP concentration",
  outcome = "AAO",
  snps = dd$SNP,
  correlation = mr_data_r2_01$ld_matrix
)

fp <- MendelianRandomization::mr_forest(
  object = inp,
  snp_estimates = TRUE,
  methods = "ivw",
  ordered = TRUE
)

res <- MendelianRandomization::mr_ivw(inp)
res_corr <- MendelianRandomization::mr_ivw(inp_corr, correl = TRUE)

# Compute summary statistics of the result -----------------------------
est <- round(res@Estimate, 3L)
lower <- round(res@CILower, 3L)
upper <- round(res@CIUpper, 3L)

results <- tibble(
  SNP = dd$SNP,
  Gene = dd$Closest_gene,
  `R^2` = dd$beta_x^2 * dd$freq_x * (1 - dd$freq_x),
  `F-statistic` = dd$beta_x^2/dd$se_x^2,
  est = est,
  lower = lower,
  upper = upper,
  pval = round(res@Pvalue, 3L)
)

primary_aao_r2_01 <- list(
  "data" = dd, "mr_input" = inp, "mr_results" = results,
  "mr_results_corr" = res_corr
)

# 4) Save results ------------------------------------------------------
save(
  primary_aao_r2_0001, primary_aao_r2_01,
  file = here::here("ccgc_crp", "analysis", "mr_results_aao.rda")
)
