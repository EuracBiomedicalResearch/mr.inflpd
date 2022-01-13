# Prepare IVs using Folkersen et al. summary data ----------------------
library(tidyverse)

# 1) Load summary stats ------------------------------------------------
load(here::here("folkersen_il1ra", "data-raw", "gwas_snps_folkersen.rda"))

instr <- snps %>%
  # Remove non biallelic variants
  dplyr::mutate(
    biallelic = str_length(Allele1 == 1) & str_length(Allele2) == 1
  ) %>%
  dplyr::filter(biallelic) %>%
  dplyr::select(-biallelic) %>%
  # Effect allele and Reference allele to upper case
  dplyr::mutate(dplyr::across(c(Allele1, Allele2), ~ str_to_upper(.))) %>%
  # Store Chromosome and position into separated columns
  dplyr::mutate(
    Chromosome = as.double(str_extract(MarkerName, "[^:]+")),
    Position = str_extract(str_extract(MarkerName, ":(.*?):"), "[0-9]+"),
    chr_pos = glue::glue("{Chromosome}_{Position}")
  ) %>%
  # Rename some variables
  dplyr::rename(
    Effect_allele = Allele1,
    Reference_allele = Allele2,
    EAF = Freq1,
    Beta = Effect,
    SE = StdErr,
    P = `P-value`
  ) %>%
  # Compute F-statistic and R^2 for each SNP
  dplyr::mutate(
    F_statistic = Beta^2/SE^2,
    R_2 = 2 * Beta^2 * EAF * (1 - EAF)
  ) %>%
  # Reorder the columns
  dplyr::select(
    Chromosome, Position, chr_pos, Effect_allele, Reference_allele,
    EAF, Beta, SE, P, F_statistic, R_2, Direction, TotalSampleSize
  )

# Save results ---------------------------------------------------------
save(instr, file = here::here("folkersen_il1ra", "data", "iv_data.rda"))
