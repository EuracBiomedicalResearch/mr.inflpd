# Import raw data for instruments --------------------------------------
library(tidyverse)
library(readxl)

# 1) Import and clean data ---------------------------------------------
instr <- read_xlsx(
  path = here::here("prins_tnfalpha", "data-raw", "instruments_tnf.xlsx")
) |>
  mutate(
    # Compute the F-statistic
    F_statistic = Beta^2/SE^2,
    # Compute the R^2
    R_2 = 2 * Beta^2 * Effect_allele_frequency * (1 - Effect_allele_frequency)
  )

# 2) Save the data -----------------------------------------------------
save(
  instr,
  file = here::here("prins_tnfalpha", "data", "iv_data.rda")
)

