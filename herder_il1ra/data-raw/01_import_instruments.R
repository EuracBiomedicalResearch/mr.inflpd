# Import raw data for instruments --------------------------------------
library(tidyverse)
library(readxl)

# 1) Import the instruments --------------------------------------------
instr <- read_xlsx(
  path = here::here("herder_il1ra", "data-raw", "instruments_herder.xlsx")
) %>%
  mutate(
    # Compute the F-statistic
    F_statistic = Beta^2/SE^2,
    # Compute the R^2
    R_2 = 2 * Beta^2 * Effect_allele_frequency * (1 - Effect_allele_frequency)
  )

# 2) Save the data -----------------------------------------------------
save(
  instr,
  file = here::here("herder_il1ra", "data", "iv_data.rda")
)


