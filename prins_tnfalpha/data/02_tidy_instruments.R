# Prepare TNF-alpha instruments data -----------------------------------
library(tidyverse)
library(LDlinkR)

# 1) Load the data -----------------------------------------------------
load(here::here("prins_tnfalpha", "data", "iv_data.rda"))

# 2) Prepare instruments data ------------------------------------------
tab_1 <- instr %>%
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
  arrange(Pvalue)

# 3) LD clumping -------------------------------------------------------
ld <- SNPclip(
  snps = tab_1$SNP[c(2, 3)],
  pop = "CEU",
  r2_threshold = "0.001",
  maf_threshold = "0.01",
  token = Sys.getenv("LDLINK_TOKEN")
)

r2_0001 <- tab_1 %>%
  filter(!SNP %in% ld[ld$Details != "Variant kept.", 1])

# 4) Save the data -----------------------------------------------------
save(
  r2_0001,
  file = here::here(
    "prins_tnfalpha", "data", "tidied_instruments.rda"
  )
)
