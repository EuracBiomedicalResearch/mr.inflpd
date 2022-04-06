# Import raw data for instruments --------------------------------------
library(tidyverse)
library(readxl)

# 1) Import data -------------------------------------------------------
instr <- read_xlsx(
  path = here::here(
    "ligthart_crp", "data-raw", "instruments_Ligthart.xlsx"
  )
) %>%
  mutate_at(
    vars(Coded_allele_frequency, Beta, SE, P_value),
    ~ as.double(.)
  ) %>%
  mutate(
    # Retrieve the reference allele
    Reference_allele = case_when(
      Effect_allele == "A" ~ "G",
      Effect_allele == "G" ~ "A",
      Effect_allele == "T" ~ "C",
      Effect_allele == "C" ~ "T",
    ),
    # Compute the F-statistic
    F_statistic = Beta^2/SE^2,
    # Compute the R^2
    R_2 = 2 * Beta^2 * Coded_allele_frequency * (1 - Coded_allele_frequency)
  ) %>%
  # Group SNPs by Chromosome and sort them in ascending order by P-values
  mutate(
    # Compute F-statisti with other formula
    F_stat_2 = R_2 * ((204402 - 2)/(1 - R_2))
  ) %>%
  group_by(Chromosome) %>%
  arrange(P_value, .by_group = TRUE) %>%
  ungroup()

# 2) Save data ---------------------------------------------------------
save(
  instr,
  file = here::here("ligthart_crp", "data", "iv_data.rda")
)


