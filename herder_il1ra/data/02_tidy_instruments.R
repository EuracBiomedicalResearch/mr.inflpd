# Tidy IL-1ra instruments ----------------------------------------------
library(tidyverse)
library(LDlinkR)

# 1) Load the data -----------------------------------------------------
load(here::here("herder_il1ra", "data", "iv_data.rda"))

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
# 3A) R2<0.001 ---------------------------------------------------------
ld <- SNPclip(
  snps = tab_1$SNP,
  pop = "CEU",
  r2_threshold = "0.001",
  maf_threshold = "0.01",
  token = Sys.getenv("LDLINK_TOKEN")
)

r2_0001 <- tab_1 %>%
  filter(!SNP %in% ld[ld$Details != "Variant kept.", 1])

# 3B) R2<0.01 ----------------------------------------------------------
ld <- SNPclip(
  snps = tab_1$SNP,
  pop = "CEU",
  r2_threshold = "0.01",
  maf_threshold = "0.01",
  token = Sys.getenv("LDLINK_TOKEN")
)

r2_001 <- tab_1 %>%
  filter(!SNP %in% ld[ld$Details != "Variant kept.", 1])

# 3C) R2<0.1 -----------------------------------------------------------
ld <- SNPclip(
  snps = tab_1$SNP,
  pop = "CEU",
  r2_threshold = "0.1",
  maf_threshold = "0.01",
  token = Sys.getenv("LDLINK_TOKEN")
)

r2_01 <- tab_1 %>%
  filter(!SNP %in% ld[ld$Details != "Variant kept.", 1])

# Same set of IVs for each LD clumping threshold.

# 3) Save into the data ------------------------------------------------
save(
  r2_0001,
  file = here::here("herder_il1ra", "data", "tidied_instruments.rda")
)
