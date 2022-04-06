# MR analysis IL-6 -> AAO ----------------------------------------------
library(tidyverse)
library(MendelianRandomization)
library(TwoSampleMR)

# 1) Load the datasets -------------------------------------------------
load(here::here("il6r_il6", "data", "mr_data_aao.rda"))
set.seed(1456)

# 2) Primary MR analysis -----------------------------------------------
dd <- mr_data$mr_dataset

inp <- MendelianRandomization::mr_input(
  bx = dd$beta_x,
  bxse = dd$se_x,
  by = dd$beta_y,
  byse = dd$se_y,
  exposure = "log IL-6 concentration",
  outcome = "Age at Onset PD",
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

primary_aao <- list(
  "data" = dd, "mr_input" = inp, "mr_results" = results
)

# 4) Save results ------------------------------------------------------
save(
  primary_aao,
  file = here::here("il6r_il6", "analysis", "results_aao.rda")
)

