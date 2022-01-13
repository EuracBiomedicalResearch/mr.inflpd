# Tidy IVs from Ahluwalia for IL-6 -------------------------------------
library(tidyverse)
library(readxl)
library(openxlsx)
library(LDlinkR)

# 1) Load imported data on IVs -----------------------------------------
load(
  here::here("ahluwalia_il6", "data", "iv_data.rda")
)

# 2) LD clumping at r2 < 0.001 -----------------------------------------
chr1 <- SNPclip(
  snps = instr %>%
    dplyr::filter(Chromosome_x == 1) %>%
    .[["SNP"]],
  pop = "CEU",
  r2_threshold = "0.001",
  maf_threshold = "0.01",
  token = Sys.getenv("LDLINK_TOKEN"),
) %>%
  dplyr::filter(Details == "Variant kept.") %>%
  .[["RS_Number"]]

chr6 <- SNPclip(
  snps = instr %>%
    dplyr::filter(Chromosome_x == 6) %>%
    .[["SNP"]],
  pop = "CEU",
  r2_threshold = "0.001",
  maf_threshold = "0.01",
  token = Sys.getenv("LDLINK_TOKEN"),
) %>%
  dplyr::filter(Details == "Variant kept.") %>%
  .[["RS_Number"]]

# 3) Ensure positive genetic associations ------------------------------
r2_0001 <- instr %>%
  dplyr::filter(SNP %in% c(chr1, chr6)) %>%
  mutate(
    # Swap effect and reference alleles if Beta < 0
    Effect_allele_1 = if_else(
      Beta < 0, Reference_allele, Effect_allele
    ),
    Reference_allele = if_else(
      Beta < 0, Effect_allele, Reference_allele
    ),
    # Compute the complementary Coded allele frequency
    EAF = if_else(Beta < 0, 1 - EAF, EAF),
    # Switch the sign of Beta
    Beta = if_else(Beta < 0, -Beta, Beta)
  ) %>%
  # Select all the columns
  dplyr::select(
    SNP:Position_x, Effect_allele_1, Reference_allele:Gene
  ) %>%
  rename(Effect_allele = Effect_allele_1) %>%
  mutate(chr_pos = glue::glue("{Chromosome_x}_{Position_x}"))

# 4) Save tidied instruments data --------------------------------------
save(
  r2_0001,
  file = here::here(
    "ahluwalia_il6", "data", "tidied_instruments.rda"
  )
)
