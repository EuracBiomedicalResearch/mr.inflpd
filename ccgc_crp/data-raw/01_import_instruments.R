# Import raw data for instruments --------------------------------------
library(tidyverse)
library(readxl)

# 1) Import and prepare data -------------------------------------------
instr <- read_xlsx(
  path = here::here("ccgc_crp", "data-raw", "instruments_ccgc.xlsx")
) %>%
  mutate_at(
    vars(Effect_allele_frequency, Beta, Lower, Upper),
    ~ as.double(.)
  ) %>%
  # From beta and CI compute SE and p-value
  mutate(
    SE = (Upper - Lower)/(1.96 * 2),
    z = Beta/SE,
    P_value = 2 * pnorm(abs(z), lower.tail = FALSE)
  ) %>%
  dplyr::select(SNP:Beta, SE, P_value, Closest_gene) %>%
  mutate(
    # Compute the F-statistic
    F_statistic = Beta^2/SE^2,
    # Compute the R^2
    R_2 = 2 * Beta^2 * Effect_allele_frequency * (1 - Effect_allele_frequency)
  )

# 2) Save the data -----------------------------------------------------
save(
  instr,
  file = here::here("ccgc_crp", "data", "iv_data.rda")
)


