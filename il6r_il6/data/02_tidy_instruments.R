# Tidy IVs from IL6R consortium for IL-6 -------------------------------
library(tidyverse)
library(readxl)

# 1) Load imported data on IVs -----------------------------------------
load(here::here("il6r_il6", "data", "iv_data.rda"))

# 2) Tidy instruments data ---------------------------------------------
r2_0001 <- instr %>%
  mutate(
    # Swap effect and reference alleles if Beta < 0
    Effect_allele_1 = if_else(
      Beta < 0, Reference_allele, Effect_allele
    ),
    Reference_allele = if_else(
      Beta < 0, Effect_allele, Reference_allele
    ),
    # Compute the complementary Coded allele frequency
    Effect_allele_frequency = if_else(
      Beta < 0, 1 - Effect_allele_frequency, Effect_allele_frequency
    ),
    # Switch the sign of Beta
    Beta = if_else(Beta < 0, -Beta, Beta)
  ) %>%
  # Select all the columns
  dplyr::select(SNP:Position, Effect_allele_1, Reference_allele:R_2) %>%
  rename(Effect_allele = Effect_allele_1) %>%
  mutate(chr_pos = glue::glue("{Chromosome}_{Position}")) %>%
  arrange(P_value)

# 2) Save into rda -----------------------------------------------------
save(
  r2_0001,
  file = here::here("il6r_il6", "data", "tidied_instruments.rda")
)
