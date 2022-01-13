# Import raw data ------------------------------------------------------
library(tidyverse)
library(readr)

ma <- read_table2(
  file = here::here(
    "aao", "data-raw", "IPDGC_AAO_GWAS_sumstats_april_2018.txt"
  )
)

# Save as rda in "data" ------------------------------------------------
save(
  ma,
  file = here::here("aao", "data", "aao2019_data.rda")
)
