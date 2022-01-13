# MR inflammation: TNF->PD ---------------------------------------------
library(tidyverse)
library(MendelianRandomization)
library(TwoSampleMR)

# 1) Load the data -----------------------------------------------------
load(here::here("prins_tnfalpha", "data", "mr_data_pd.rda"))
source(here::here("prins_tnfalpha", "analysis", "mr_misc_functions.R"))
set.seed(1456)

# 2) Forest plots of Wald estimates ------------------------------------
df <- mr_data$mr_dataset |>
  mutate(
    f_stat = beta_x^2/se_x^2,
    r2 = 2 * beta_x^2 * freq_x * (1 - freq_x)
  )

inp <- MendelianRandomization::mr_input(
  bx = df$beta_x,
  bxse = df$se_x,
  by = df$beta_y,
  byse = df$se_y,
  exposure = "log TNF concentration",
  outcome = "log-odds PD",
  snps = df$SNP
)

fp <- MendelianRandomization::mr_forest(
  object = inp,
  snp_estimates = TRUE,
  methods = "ivw",
  ordered = TRUE
)

# 3) Primary MR analysis -----------------------------------------------
res <- MendelianRandomization::mr_ivw(inp, model = "fixed")

# Compute summary statistics of the result -----------------------------
est <- round(res@Estimate, 3L)
lower <- round(res@CILower, 3L)
upper <- round(res@CIUpper, 3L)

results <- tibble(
  `N IVs` = length(df$SNP),
  `R2` = sum(df$r2),
  `F-statistic` = glue::glue(
    "{round(median(df$f_stat), 2L)} [{round(min(df$f_stat), 2L)}; {round(max(df$f_stat), 2L)}]"
  ),
  est = est,
  lower = lower,
  upper = upper,
  pval = round(res@Pvalue, 3L)
)

primary_pd <- list(
  "data" = df, "mr_input" = inp, "mr_results" = results,
  "mr_model" = res
)

# 4) Sensitivity analysis ----------------------------------------------
df <- mr_data$mr_dataset %>%
  mutate(
    f_stat = beta_x^2/se_x^2,
    r2 = 2 * beta_x^2 * freq_x * (1 - freq_x)
  )

q_stat_res <- q_stat_fun(data = df)

df_trim_stat <- q_stat_res$data_trim %>%
  mutate(
    f_stat = beta_x^2/se_x^2,
    r2 = 2 * beta_x^2 * freq_x * (1 - freq_x)
  )

# Exclusion based on Phenoscanner look-up ------------------------------
# Known risk factors for PD:
# - Smoking
# - Rheumatoid arthritis
# - Cognitive performances
# - Educational attainment
# - Iron
# - Total cholesterol
# - BMI
# - Urate in serum
# - LDL cholesterol
# - Alcohol

snps_lookup <- phenoscanner(
  snpquery = df_trim_stat$SNP, pvalue = 5e-08
)$results %>%
  dplyr::filter(ancestry == "European") %>%
  # Group traits
  mutate(
    trait_macro = case_when(
      str_detect(trait, "(BMI)|(Body Mass Index)|(body mass index)|(Body mass index)|(obesity)|(Obesity)|(Overweight)|(overweight)") ~ "BMI",
      str_detect(trait, "(Diabetes)|(diabetest)") ~ "Diabetes",
      str_detect(trait, "(Cholesterol)|(cholesterol)|(LDL)|(HDL)") ~ "Cholesterol",
      str_detect(trait, "(Smoking)|(smoking)|(Smoker)|(smoker)") ~ "Smoking",
      str_detect(trait, "(Alcoho)|(alcoho)|(Drink)|(drink)") ~ "Alcohol",
      str_detect(trait, "(Urate)|(urate)") ~ "Urate",
      str_detect(trait, "(Rheumatoid)|(rheumatoid)|(arthritis)|(Arthritis)") ~ "Rheumatoid arthritis",
      str_detect(trait, "(Cognitive)|(cognitive)|(Educational)|(educational)") ~ "Cognitive/Educational",
      TRUE ~ "other"
    )
  ) %>%
  # Take only SNPs associated with known risk factors
  dplyr::filter(trait_macro != "other")

pleiotropic_snps <- snps_lookup$snp

# Prepare data ---------------------------------------------------------
df_trim <- df_trim_stat %>%
  dplyr::filter(!SNP %in% pleiotropic_snps)

inp <- mr_input(
  bx = df_trim$beta_x,
  bxse = df_trim$se_x,
  by = df_trim$beta_y,
  byse = df_trim$se_y,
  exposure = "log TNF concentration",
  outcome = "log-odds PD",
  snps = df_trim$SNP
)

fp_sens <- MendelianRandomization::mr_forest(
  object = inp,
  snp_estimates = TRUE,
  methods = "ivw",
  ordered = TRUE
)

res <- MendelianRandomization::mr_ivw(inp, model = "fixed")

# Compute R^2, F-stat and I^2 ------------------------------------------
est <- round(res@Estimate, 3L)
lower <- round(res@CILower, 3L)
upper <- round(res@CIUpper, 3L)

results <- tibble(
  `N IVs` = length(df$SNP),
  `R2` = sum(df$r2),
  `F-statistic` = glue::glue(
    "{round(median(df$f_stat), 2L)} [{round(min(df$f_stat), 2L)}; {round(max(df$f_stat), 2L)}]"
  ),
  est = est,
  lower = lower,
  upper = upper,
  pval = round(res@Pvalue, 3L)
)

sensitivity_pd <- list(
  "data" = df_trim, "mr_input" = inp, "forest_plot" = fp_sens,
  "mr_results" = results,
  "mr_model" = res,
  "pleiotropic_snps" = list(
    "stat" = df_trim_stat, "bio" = df_trim, "lookup" = snps_lookup
  )
)

# 3) Save results ------------------------------------------------------
save(
  primary_pd, sensitivity_pd,
  file = here::here("prins_tnfalpha", "analysis", "mr_results_pd.rda")
)

