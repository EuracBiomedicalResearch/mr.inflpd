# Organize results for the paper ---------------------------------------

# The aim of this script is to organized the results of primary and
# sensitivity MR analysis for tables and figures of manuscript and
# supplementary

# 1) CRP ---------------------------------------------------------------
## 1A) Ligthart et al. -------------------------------------------------
load(here::here("ligthart_crp", "analysis", "mr_results_pd.rda"))
load(here::here("ligthart_crp", "analysis", "mr_results_aao.rda"))
ligthart_results <- list(
  "pd" = list("primary" = primary_pd, "sens" = sensitivity_pd),
  "aao" = list(
    "primary" = primary_aao, "sens" = sensitivity_aao
  )
)

rm(primary_pd, primary_aao, sensitivity_pd, sensitivity_aao)

# 1B) CCGC -------------------------------------------------------------
load(here::here("ccgc_crp", "analysis", "mr_results_pd.rda"))
load(here::here("ccgc_crp", "analysis", "mr_results_aao.rda"))

ccgc_results <- list(
  "pd" = list("primary" = primary_pd), "aao" = list("primary" = primary_aao)
)

rm(primary_pd, primary_aao)

# 2) IL-6 --------------------------------------------------------------
# 2A) Ahluwalia et al. -------------------------------------------------
load(here::here("ahluwalia_il6", "analysis", "mr_results_pd.rda"))
load(here::here("ahluwalia_il6", "analysis", "mr_results_aao.rda"))

ahluwalia_results <- list(
  "pd" = list("primary" = primary_pd, "sens" = sensitivity_pd),
  "aao" = list("primary" = primary_aao, "sens" = sensitivity_aao)
)

rm(primary_pd, primary_aao, sensitivity_pd, sensitivity_aao)

## 2B) IL6R Genetic consortium -----------------------------------------
load(here::here("il6r_il6", "analysis", "results_pd.rda"))
load(here::here("il6r_il6", "analysis", "results_aao.rda"))

il6r_results <- list(
  "pd" = list("primary" = primary_pd), "aao" = list("primary" = primary_aao)
)

rm(primary_pd, primary_aao)

## 2C) Folkersen et al. ------------------------------------------------
load(here::here("folkersen_il6", "analysis", "mr_results_pd.rda"))
load(here::here("folkersen_il6", "analysis", "mr_results_aao.rda"))

folkersen_results <- list(
  "pd" = list("primary" = primary_pd, "sens" = sensitivity_pd),
  "aao" = list("primary" = primary_aao, "sens" = sensitivity_aao)
)

rm(primary_pd, primary_aao, sensitivity_pd, sensitivity_aao)

# 3) IL-1ra ------------------------------------------------------------
## 3A) Folkersen -------------------------------------------------------
load(here::here("folkersen_il1ra", "analysis", "mr_results_pd.rda"))
load(here::here("folkersen_il1ra", "analysis", "mr_results_aao.rda"))

il1ra_results <- list(
  "pd" = list("primary" = primary_pd, "sens" = sensitivity_pd),
  "aao" = list("primary" = primary_aao, "sens" = sensitivity_aao)
)

rm(primary_pd, primary_aao, sensitivity_pd,sensitivity_aao)

## 3B) Herder ----------------------------------------------------------
load(here::here("herder_il1ra", "analysis", "mr_results_pd.rda"))
load(here::here("herder_il1ra", "analysis", "mr_results_aao.rda"))

herder_results <- list(
  "pd" = list("primary" = primary_pd),
  "aao" = list("primary" = primary_aao)
)

rm(primary_pd, primary_aao)

# 4) TNF-alpha --------------------------------------------------------
load(here::here("prins_tnfalpha", "analysis", "mr_results_pd.rda"))
load(here::here("prins_tnfalpha", "analysis", "mr_results_aao.rda"))

tnf_results <- list(
  "pd" = list("primary" = primary_pd, "sens" = sensitivity_pd),
  "aao" = list("primary" = primary_aao, "sens" = sensitivity_aao)
)

rm(primary_pd, primary_aao, sensitivity_pd, sensitivity_aao)

# 5) Save results ------------------------------------------------------
save(
  ligthart_results,
  ccgc_results,
  ahluwalia_results,
  il6r_results,
  folkersen_results,
  il1ra_results,
  herder_results,
  tnf_results,
  file = here::here("all", "paper_results.rda")
)

