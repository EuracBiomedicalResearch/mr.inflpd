# Import raw data ------------------------------------------------------
library(tidyverse)
library(readr)

ma <- read_table2(
  file = here::here(
    "pd", "data-raw", "nallsEtAl2019_excluding23andMe_allVariants.tab"
  )
)

# Save as rda in "data" ------------------------------------------------
save(
  ma,
  file = here::here("pd", "data", "nalls2019_data.rda")
)
