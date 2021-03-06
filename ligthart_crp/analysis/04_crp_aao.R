# MR analysis CRP-AAO --------------------------------------------------
library(tidyverse)
library(MendelianRandomization)
library(TwoSampleMR)
library(MRMix)
library(MRPRESSO)
library(MRPracticals)

# 1) Prepare the data --------------------------------------------------
load(here::here("ligthart_crp", "data", "mr_data_aao.rda"))
source(here::here("ligthart_crp", "analysis", "mr_misc_functions.R"))
set.seed(1456)

# A) R2<0.001 ----------------------------------------------------------
# 2) Forest plot wald estimates ----------------------------------------
df <- mr_data_r2_0001$mr_dataset

inp <- MendelianRandomization::mr_input(
  bx = df$beta_x,
  bxse = df$se_x,
  by = df$beta_y,
  byse = df$se_y,
  exposure = "log CRP concentration",
  outcome = "AAO",
  snps = df$SNP
)

fp <- MendelianRandomization::mr_forest(
  object = inp,
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
df <- mr_data_r2_0001$mr_dataset %>%
  mutate(
    f_stat = beta_x^2/se_x^2,
    r2 = 2 * beta_x^2 * freq_x * (1 - freq_x)
  )

res <- mr_fun(data = df, input = inp)

# Compute R^2, F-stat and I^2 ------------------------------------------
stats <- tibble(
  r2 = sum(df$r2),
  f_stat = glue::glue(
    "{round(median(df$f_stat), 2L)} [{round(min(df$f_stat), 2L)}; {round(max(df$f_stat), 2L)}]"
  ),
  q_stat = res$ivw$het@Heter.Stat[1],
  pval_q_stat = if_else(
    res$ivw$het@Heter.Stat[2] < 0.001, "<0.001",
    as.character(round(res$ivw$het@Heter.Stat[2], 3L))
  ),
  i2_stat = i2_fun(df),
  pval_outlier_presso = res$presso$models$`MR-PRESSO results`$`Global Test`$Pvalue
)

primary_aao_r2_0001 <- list(
  "data" = df, "mr_input" = inp, "mr_results" = res,
  "statistics" = stats
)

# 4) Sensitivity analyses ----------------------------------------------
df <- mr_data_r2_0001$mr_dataset %>%
  mutate(
    f_stat = beta_x^2/se_x^2,
    r2 = 2 * beta_x^2 * freq_x * (1 - freq_x)
  )

# 4A) Exclude pleitropic IVs -------------------------------------------
# Exclusion based on extreme Q-stat (Bonferroni) -----------------------
q_stat_res <- q_stat_fun(data = df)

df_trim_stat <- q_stat_res$data_trim %>%
  mutate(
    f_stat = beta_x^2/se_x^2,
    r2 = 2 * beta_x^2 * freq_x * (1 - freq_x)
  )

# Exclusion based on Phenoscanner look-up ------------------------------
# Known risk factors for aao:
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
df_trim <- df_trim_stat %>%
  dplyr::filter(!SNP %in% pleiotropic_snps)

inp <- mr_input(
  bx = df_trim$beta_x,
  bxse = df_trim$se_x,
  by = df_trim$beta_y,
  byse = df_trim$se_y,
  exposure = "log CRP concentration",
  outcome = "AAO",
  snps = df_trim$SNP
)

fp_sens <- MendelianRandomization::mr_forest(
  object = inp,
  snp_estimates = TRUE,
  methods = "ivw",
  ordered = TRUE
)

res <- mr_fun(data = df_trim, input = inp)

# Compute R^2, F-stat and I^2 ------------------------------------------
stats <- tibble(
  r2 = sum(df_trim$r2),
  f_stat = glue::glue(
    "{round(median(df_trim$f_stat), 2L)} [{round(min(df_trim$f_stat), 2L)}; {round(max(df_trim$f_stat), 2L)}]"
  ),
  q_stat = res$ivw$het@Heter.Stat[1],
  pval_q_stat = if_else(
    res$ivw$het@Heter.Stat[2] < 0.001, "<0.001",
    as.character(
      round(res$ivw$het@Heter.Stat[2], 3L)
    )
  ),
  i2_stat = i2_fun(df_trim),
  pval_outlier_presso = res$presso$models$`MR-PRESSO results`$`Global Test`$Pvalue
)

sensitivity_aao_r2_0001 <- list(
  "data" = df_trim, "mr_input" = inp, "forest_plot" = fp_sens,
  "mr_results" = res, "statistics" = stats,
  "q_stat_plot" = q_stat_res,
  "pleiotropic_snps" = list(
    "stat" = df_trim_stat, "bio" = df_trim, "lookup" = snps_lookup
  )
)

# B) R2<0.01 -----------------------------------------------------------
# 2) Forest plot wald estimates ----------------------------------------
df <- mr_data_r2_001$mr_dataset

inp <- MendelianRandomization::mr_input(
  bx = df$beta_x,
  bxse = df$se_x,
  by = df$beta_y,
  byse = df$se_y,
  exposure = "log CRP concentration",
  outcome = "AAO",
  snps = df$SNP
)

inp_corr <- MendelianRandomization::mr_input(
  bx = df$beta_x,
  bxse = df$se_x,
  by = df$beta_y,
  byse = df$se_y,
  exposure = "log CRP concentration",
  outcome = "AAO",
  snps = df$SNP,
  correlation = mr_data_r2_001$ld_mat
)

fp <- MendelianRandomization::mr_forest(
  object = inp,
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
df <- mr_data_r2_001$mr_dataset %>%
  mutate(
    f_stat = beta_x^2/se_x^2,
    r2 = 2 * beta_x^2 * freq_x * (1 - freq_x)
  )

res <- mr_fun(data = df, input = inp)
res_corr <- MendelianRandomization::mr_ivw(
  inp_corr, model = "random", correl = TRUE
)

# Compute R^2, F-stat and I^2 ------------------------------------------
stats <- tibble(
  r2 = sum(df$r2),
  f_stat = glue::glue(
    "{round(median(df$f_stat), 2L)} [{round(min(df$f_stat), 2L)}; {round(max(df$f_stat), 2L)}]"
  ),
  q_stat = res$ivw$het@Heter.Stat[1],
  pval_q_stat = if_else(
    res$ivw$het@Heter.Stat[2] < 0.001, "<0.001",
    as.character(round(res$ivw$het@Heter.Stat[2], 3L))
  ),
  i2_stat = i2_fun(df),
  pval_outlier_presso = res$presso$models$`MR-PRESSO results`$`Global Test`$Pvalue
)

primary_aao_r2_001 <- list(
  "data" = df, "mr_input" = inp, "mr_results" = res,
  "statistics" = stats, "mr_results_corr" = res_corr
)

# 4) Sensitivity analyses ----------------------------------------------
df <- mr_data_r2_001$mr_dataset %>%
  mutate(
    f_stat = beta_x^2/se_x^2,
    r2 = 2 * beta_x^2 * freq_x * (1 - freq_x)
  )

# 4A) Exclude pleitropic IVs -------------------------------------------
# Exclusion based on extreme Q-stat (Bonferroni) -----------------------
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
df_trim <- df_trim_stat %>%
  dplyr::filter(!SNP %in% pleiotropic_snps)

inp <- mr_input(
  bx = df_trim$beta_x,
  bxse = df_trim$se_x,
  by = df_trim$beta_y,
  byse = df_trim$se_y,
  exposure = "log CRP concentration",
  outcome = "AAO",
  snps = df_trim$SNP
)

inp <- mr_input(
  bx = df_trim$beta_x,
  bxse = df_trim$se_x,
  by = df_trim$beta_y,
  byse = df_trim$se_y,
  exposure = "log CRP concentration",
  outcome = "AAO",
  snps = df_trim$SNP,
  correlation = mr_data_r2_001$ld_mat[df_trim$SNP, df_trim$SNP]
)

fp_sens <- MendelianRandomization::mr_forest(
  object = inp,
  snp_estimates = TRUE,
  methods = "ivw",
  ordered = TRUE
)

res <- mr_fun(data = df_trim, input = inp)
res_corr <- MendelianRandomization::mr_ivw(
  inp_corr, model = "random", correl = TRUE
)

# Compute R^2, F-stat and I^2 ------------------------------------------
stats <- tibble(
  r2 = sum(df_trim$r2),
  f_stat = glue::glue(
    "{round(median(df_trim$f_stat), 2L)} [{round(min(df_trim$f_stat), 2L)}; {round(max(df_trim$f_stat), 2L)}]"
  ),
  q_stat = res$ivw$het@Heter.Stat[1],
  pval_q_stat = if_else(
    res$ivw$het@Heter.Stat[2] < 0.001, "<0.001",
    as.character(
      round(res$ivw$het@Heter.Stat[2], 3L)
    )
  ),
  i2_stat = i2_fun(df_trim),
  pval_outlier_presso = res$presso$models$`MR-PRESSO results`$`Global Test`$Pvalue
)

sensitivity_aao_r2_001 <- list(
  "data" = df_trim, "mr_input" = inp, "forest_plot" = fp_sens,
  "mr_results" = res, "statistics" = stats,
  "q_stat_plot" = q_stat_res,
  "pleiotropic_snps" = list(
    "stat" = df_trim_stat, "bio" = df_trim, "lookup" = snps_lookup
  ),
  "mr_results_corr" = res_corr
)

# C) R2<0.1 ------------------------------------------------------------
# 2) Forest plot wald estimates ----------------------------------------
df <- mr_data_r2_01$mr_dataset

inp <- MendelianRandomization::mr_input(
  bx = df$beta_x,
  bxse = df$se_x,
  by = df$beta_y,
  byse = df$se_y,
  exposure = "log CRP concentration",
  outcome = "AAO",
  snps = df$SNP
)

inp_corr <- MendelianRandomization::mr_input(
  bx = df$beta_x,
  bxse = df$se_x,
  by = df$beta_y,
  byse = df$se_y,
  exposure = "log CRP concentration",
  outcome = "AAO",
  snps = df$SNP,
  correlation = mr_data_r2_01$ld_mat
)

fp <- MendelianRandomization::mr_forest(
  object = inp,
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
df <- mr_data_r2_01$mr_dataset %>%
  mutate(
    f_stat = beta_x^2/se_x^2,
    r2 = 2 * beta_x^2 * freq_x * (1 - freq_x)
  )

res <- mr_fun(data = df, input = inp)
res_corr <- MendelianRandomization::mr_ivw(
  inp_corr, model = "random", correl = TRUE
)

# Compute R^2, F-stat and I^2 ------------------------------------------
stats <- tibble(
  r2 = sum(df$r2),
  f_stat = glue::glue(
    "{round(median(df$f_stat), 2L)} [{round(min(df$f_stat), 2L)}; {round(max(df$f_stat), 2L)}]"
  ),
  q_stat = res$ivw$het@Heter.Stat[1],
  pval_q_stat = if_else(
    res$ivw$het@Heter.Stat[2] < 0.001, "<0.001",
    as.character(round(res$ivw$het@Heter.Stat[2], 3L))
  ),
  i2_stat = i2_fun(df),
  pval_outlier_presso = res$presso$models$`MR-PRESSO results`$`Global Test`$Pvalue
)

primary_aao_r2_01 <- list(
  "data" = df, "mr_input" = inp, "mr_results" = res,
  "statistics" = stats, "mr_results_corr" = res_corr
)

# 4) Sensitivity analyses ----------------------------------------------
df <- mr_data_r2_01$mr_dataset %>%
  mutate(
    f_stat = beta_x^2/se_x^2,
    r2 = 2 * beta_x^2 * freq_x * (1 - freq_x)
  )

# 4A) Exclude pleitropic IVs -------------------------------------------
# Exclusion based on extreme Q-stat (Bonferroni) -----------------------
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
df_trim <- df_trim_stat %>%
  dplyr::filter(!SNP %in% pleiotropic_snps)

inp <- mr_input(
  bx = df_trim$beta_x,
  bxse = df_trim$se_x,
  by = df_trim$beta_y,
  byse = df_trim$se_y,
  exposure = "log CRP concentration",
  outcome = "AAO",
  snps = df_trim$SNP
)

inp_corr <- mr_input(
  bx = df_trim$beta_x,
  bxse = df_trim$se_x,
  by = df_trim$beta_y,
  byse = df_trim$se_y,
  exposure = "log CRP concentration",
  outcome = "AAO",
  snps = df_trim$SNP,
  correlation = mr_data_r2_01$ld_mat[df_trim$SNP, df_trim$SNP]
)

fp_sens <- MendelianRandomization::mr_forest(
  object = inp,
  snp_estimates = TRUE,
  methods = "ivw",
  ordered = TRUE
)

res <- mr_fun(data = df_trim, input = inp)
res_corr <- MendelianRandomization::mr_ivw(
  inp_corr, model = "random", correl = TRUE
)

# Compute R^2, F-stat and I^2 ------------------------------------------
stats <- tibble(
  r2 = sum(df_trim$r2),
  f_stat = glue::glue(
    "{round(median(df_trim$f_stat), 2L)} [{round(min(df_trim$f_stat), 2L)}; {round(max(df_trim$f_stat), 2L)}]"
  ),
  q_stat = res$ivw$het@Heter.Stat[1],
  pval_q_stat = if_else(
    res$ivw$het@Heter.Stat[2] < 0.001, "<0.001",
    as.character(
      round(res$ivw$het@Heter.Stat[2], 3L)
    )
  ),
  i2_stat = i2_fun(df_trim),
  pval_outlier_presso = res$presso$models$`MR-PRESSO results`$`Global Test`$Pvalue
)

sensitivity_aao_r2_01 <- list(
  "data" = df_trim, "mr_input" = inp, "forest_plot" = fp_sens,
  "mr_results" = res, "statistics" = stats,
  "q_stat_plot" = q_stat_res,
  "pleiotropic_snps" = list(
    "stat" = df_trim_stat, "bio" = df_trim, "lookup" = snps_lookup
  ),
  "mr_results_corr" = res_corr
)

# 3) Save results ------------------------------------------------------
save(
  primary_aao_r2_0001, sensitivity_aao_r2_0001,
  primary_aao_r2_001, sensitivity_aao_r2_001,
  primary_aao_r2_01, sensitivity_aao_r2_01,
  file = here::here("ligthart_crp", "analysis", "mr_results_aao.rda")
)

