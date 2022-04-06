# MR analysis IL-6 -> PD -----------------------------------------------
library(tidyverse)
library(MendelianRandomization)
library(TwoSampleMR)
library(MRMix)
library(MRPRESSO)
library(MRPracticals)
library(LDlinkR)

# 1) Load the datasets -------------------------------------------------
load(here::here("folkersen_il6", "data", "mr_data_pd.rda"))
source(here::here("folkersen_il6", "analysis", "mr_misc_functions.R"))
set.seed(1456)

# Same instruments for r2<0.001 and r2<0.01

# A) R2<0.001 and R2<0.01 ----------------------------------------------
# 2) Forest plot wald estimates ----------------------------------------
dd <- mr_data_r2_0001$mr_dataset

i <- MendelianRandomization::mr_input(
  bx = dd$beta_x,
  bxse = dd$se_x,
  by = dd$beta_y,
  byse = dd$se_y,
  exposure = "IL-6 concentration",
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

primary_pd_r2_0001 <- list(
  "data" = df, "mr_input" = inp, "mr_results" = res,
  "statistics" = stats
)

# 4) Pleiotropy scenarios ----------------------------------------------
df <- dd %>%
  mutate(
    f_stat = beta_x^2/se_x^2,
    r2 = 2 * beta_x^2 * freq_x * (1 - freq_x)
  )

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
  # Remove rs4959106 in HLA-DQA1 locus (highly pleiotropic SNP)
  dplyr::filter(SNP != "rs4959106")

inp <- mr_input(
  bx = df_trim$beta_x,
  bxse = df_trim$se_x,
  by = df_trim$beta_y,
  byse = df_trim$se_y,
  exposure = "log IL-6",
  outcome = "log-odds",
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

sensitivity_pd_r2_0001 <- list(
  "data" = df_trim, "mr_input" = inp, "forest_plot" = fp_sens,
  "mr_results" = res, "statistics" = stats,
  "pleiotropic_snps" = snps_lookup
)

# B) R2<0.1 ------------------------------------------------------------
# 2) Forest plot wald estimates ----------------------------------------
dd <- mr_data_r2_01$mr_dataset

i <- MendelianRandomization::mr_input(
  bx = dd$beta_x,
  bxse = dd$se_x,
  by = dd$beta_y,
  byse = dd$se_y,
  exposure = "IL-6 concentration",
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

# 3) Primary and secondary MR analysis ---------------------------------
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

res <- mr_fun(data = df, input = inp)

# Analysis with correlated variants
inp_corr <- mr_input(
  bx = df$beta_x,
  bxse = df$se_x,
  by = df$beta_y,
  byse = df$se_y,
  exposure = "log IL-6",
  outcome = "log-odds PD",
  snps = df$SNP,
  correlation = mr_data_r2_01$ld_mat
)

res_corr <- MendelianRandomization::mr_ivw(
  inp_corr, model = "random", correl = TRUE
)

# Compute R^2, F-stat and I^2 ------------------------------------------
stats <- tibble(
  SNP = nrow(df),
  r2 = sum(df$r2),
  f_stat = glue::glue(
    "{round(median(df$f_stat), 3L)} [{round(min(df$f_stat), 3L)}; {round(max(df$f_stat), 3L)}]"
  ),
  q_stat = round(res$ivw$het@Heter.Stat[1], 3L),
  pval_q_stat = if_else(
    res$ivw$het@Heter.Stat[2] < 0.001, "<0.001",
    as.character(round(res$ivw$het@Heter.Stat[2], 3L))
  ),
  i2_stat = i2_fun(df)
)

primary_pd_r2_01 <- list(
  "data" = df, "mr_input" = inp, "mr_results" = res,
  "statistics" = stats, "mr_results_corr" = res_corr
)

# 4) Pleiotropy scenarios ----------------------------------------------
df <- dd %>%
  mutate(
    f_stat = beta_x^2/se_x^2,
    r2 = 2 * beta_x^2 * freq_x * (1 - freq_x)
  )

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
  # Remove rs4959106 in HLA-DQA1 locus (highly pleiotropic SNP) and
  # rs2228145
  dplyr::filter(!SNP %in% c("rs4959106", "rs2228145"))

inp <- mr_input(
  bx = df_trim$beta_x,
  bxse = df_trim$se_x,
  by = df_trim$beta_y,
  byse = df_trim$se_y,
  exposure = "log IL-6",
  outcome = "log-odds",
  snps = df_trim$SNP
)

fp_sens <- MendelianRandomization::mr_forest(
  object = inp,
  snp_estimates = TRUE,
  methods = "ivw",
  ordered = TRUE
)

res <- MendelianRandomization::mr_ivw(inp, model = "fixed")

inp_corr <- mr_input(
  bx = df_trim$beta_x,
  bxse = df_trim$se_x,
  by = df_trim$beta_y,
  byse = df_trim$se_y,
  exposure = "log IL-6",
  outcome = "log-odds",
  snps = df_trim$SNP,
  correlation = mr_data_r2_01$ld_mat[df_trim$SNP, df_trim$SNP]
)

res_corr <- MendelianRandomization::mr_ivw(
  inp_corr, model = "fixed", correl = TRUE
)

# Compute R^2, F-stat and I^2 ------------------------------------------
stats <- tibble(
  SNP = nrow(df_trim),
  r2 = sum(df_trim$r2),
  f_stat = glue::glue(
    "{round(median(df_trim$f_stat), 3L)} [{round(min(df_trim$f_stat), 3L)}; {round(max(df_trim$f_stat), 3L)}]"
  )
)

sensitivity_pd_r2_01 <- list(
  "data" = df_trim, "mr_input" = inp, "forest_plot" = fp_sens,
  "mr_results" = res, "statistics" = stats,
  "pleiotropic_snps" = snps_lookup, "mr_results_corr" = res_corr
)

# C) Save results ------------------------------------------------------
save(
  primary_pd_r2_0001, sensitivity_pd_r2_0001,
  primary_pd_r2_01, sensitivity_pd_r2_01,
  file = here::here(
    "folkersen_il6", "analysis", "mr_results_pd.rda"
  )
)

