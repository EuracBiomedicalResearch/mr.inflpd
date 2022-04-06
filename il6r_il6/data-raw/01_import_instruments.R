# Import raw data for instruments --------------------------------------
library(tidyverse)
library(readxl)

instr <- read_xlsx(
  path = here::here("il6r_il6", "data-raw", "instruments_il6r.xlsx")
) %>%
  mutate_at(
    vars(Effect_allele_frequency, Beta, Lower, Upper),
    ~ as.double(.)
  ) %>%
  # From beta and CI compute SE and p-value
  mutate(
    SE = (Upper - Lower)/(1.96 * 2),
    z = Beta/SE
  ) %>%
  dplyr::select(SNP:Beta, SE, P_value, Closest_gene) |>
  mutate(
    # Compute the F-statistic
    F_statistic = Beta^2/SE^2,
    # Compute the R^2
    R_2 = 2 * Beta^2 * Effect_allele_frequency * (1 - Effect_allele_frequency)
  )

# Save instruments data ------------------------------------------------
save(
  instr,
  file = here::here("il6r_il6", "data", "iv_data.rda")
)


