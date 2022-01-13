# MR analysis IL-6 -> PD -----------------------------------------------
library(tidyverse)
library(MendelianRandomization)
library(TwoSampleMR)
library(MRMix)
library(MRPRESSO)
library(MRPracticals)
library(LDlinkR)

# 1) Load the datasets -------------------------------------------------
load(here::here("ahluwalia_il6", "data", "mr_data_pd.rda"))
source(here::here("ahluwalia_il6", "analysis", "mr_misc_functions.R"))
set.seed(1456)

# 2) Forest plot wald estimates ----------------------------------------
dd <- mr_data$mr_dataset

i <- MendelianRandomization::mr_input(
  bx = dd$beta_x,
  bxse = dd$se_x,
  by = dd$beta_y,
  byse = dd$se_y,
  exposure = "log IL-6 concentration",
  outcome = "log-odds PD",
  snps = dd$SNP
)

fp <- MendelianRandomization::mr_forest(
  object = i,
  snp_estimates = TRUE,
  methods = "ivw",
  ordered = TRUE
)

fp_res <- list(
  "plot" = fp,
  "wald_est" = fp$data %>%
    dplyr::select(snps, estimates) %>%
    dplyr::slice(-nrow(.))
)

# 3) Primary MR analysis -----------------------------------------------
df <- dd %>%
  mutate(
    f_stat = beta_x^2/se_x^2,
    r2 = 2 * beta_x^2 * freq_x * (1 - freq_x)
  )

inp <- mr_input(
  bx = df$beta_x,
  bxse = df$se_x,
  by = df$beta_y,
  byse = df$se_y,
  exposure = "log IL-6",
  outcome = "log-odds PD",
  snps = df$SNP
)

res <- MendelianRandomization::mr_ivw(inp)

# Compute R^2, F-stat and I^2 ------------------------------------------
stats <- tibble(
  SNP = glue::glue("{df$SNP[1]}, {df$SNP[2]}"),
  r2 = sum(df$r2),
  f_stat = glue::glue(
    "{round(df$f_stat[1], 3L)}, {round(df$f_stat[2], 3L)}"
  ),
  q_stat = round(res@Heter.Stat[1], 3L),
  pval_q_stat = if_else(
    res@Heter.Stat[2] < 0.001, "<0.001",
    as.character(round(res@Heter.Stat[2], 3L))
  ),
  i2_stat = i2_fun(df)
)

primary_pd <- list(
  "data" = df, "mr_input" = inp, "mr_results" = res,
  "statistics" = stats
)

# 4) Pleiotropy scenarios ----------------------------------------------
df <- dd %>%
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
  snpquery = df$SNP, pvalue = 5e-08
)$results %>%
  dplyr::filter(ancestry == "European") %>%
  dplyr::select(snp, trait, hg19_coordinates, study, pmid, beta, se, p) %>%
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
df_trim <- df %>%
  dplyr::filter(!SNP %in% pleiotropic_snps)

inp <- mr_input(
  bx = df_trim$beta_x,
  bxse = df_trim$se_x,
  by = df_trim$beta_y,
  byse = df_trim$se_y,
  exposure = "log IL-6",
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
stats <- tibble(
  r2 = sum(df_trim$r2),
  f_stat = glue::glue(
    "{round(median(df_trim$f_stat), 2L)} [{round(min(df_trim$f_stat), 2L)}; {round(max(df_trim$f_stat), 2L)}]"
  )
)

sensitivity_pd <- list(
  "data" = df_trim, "mr_input" = inp, "forest_plot" = fp_sens,
  "mr_results" = res, "statistics" = stats,
  "pleiotropic_snps" = list(
    "stat" = df_trim, "bio" = df_trim, "lookup" = snps_lookup
  )
)

# 5) Save results ------------------------------------------------------
save(
  primary_pd, sensitivity_pd,
  file = here::here(
    "ahluwalia_il6", "analysis", "mr_results_pd.rda"
  )
)

