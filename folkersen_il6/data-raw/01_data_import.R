# Import GWAS results from Folkersen et al. (2020) ---------------------
library(tidyverse)
library(readr)
library(janitor)

# 1) Read summary statistics -------------------------------------------
snps <- read_table2(
  gzfile(here::here("folkersen_il6", "data-raw", "IL-6.txt.gz"))
) |>
  # Select only genome-wide significant data
  dplyr::filter(`P-value` < 5e-08)

# 2) Save the data -----------------------------------------------------
save(
  snps,
  file = here::here(
    "folkersen_il6", "data-raw", "gwas_snps_folkersen.rda"
  )
)
