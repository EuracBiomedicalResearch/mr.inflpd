# Tidy CRP instruments data --------------------------------------------
library(tidyverse)
library(LDlinkR)

# 1) Load the data -----------------------------------------------------
load(here::here("ligthart_crp", "data", "iv_data.rda"))

# 2) Prepare the data --------------------------------------------------
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
    Coded_allele_frequency = if_else(
      Beta < 0, 1 - Coded_allele_frequency, Coded_allele_frequency
    ),
    # Switch the sign of Beta
    Beta = if_else(Beta < 0, -Beta, Beta)
  ) %>%
  # Select all the columns
  dplyr::select(SNP:Position, Effect_allele_1, Reference_allele:R_2) %>%
  rename(Effect_allele = Effect_allele_1) %>%
  mutate(chr_pos = glue::glue("{Chromosome}_{Position}"))

# 3) LD clumping (r2<0.001) --------------------------------------------
list_chrs <- unique(tab_1$Chromosome)

ld_clump <- map_dfr(
  .x = list_chrs,
  ~ {

    # Select SNPs in the Chromosome
    seq_snps <- tab_1 |>
      dplyr::filter(Chromosome == .x)

    if(nrow(seq_snps) > 1) {

      # Perform LD clumping
      ld <- SNPclip(
        snps = seq_snps$SNP,
        pop = "CEU",
        r2_threshold = "0.001",
        maf_threshold = "0.01",
        token = Sys.getenv("LDLINK_TOKEN")
      )

      retained_snps <- ld |>
        dplyr::filter(Details == "Variant kept.") |>
        pull(RS_Number)

      # Final tibble with retained SNPs
      res <- seq_snps |>
        dplyr::filter(SNP %in% retained_snps)

    } else {

      res <- seq_snps

    }

    message(glue::glue("Chromosome {.x} completed!"))

    res

  }
)

r2_0001 <- ld_clump %>%
  mutate(
    # Swap effect and reference alleles if Beta < 0
    Effect_allele_1 = if_else(
      Beta < 0, Reference_allele, Effect_allele
    ),
    Reference_allele = if_else(
      Beta < 0, Effect_allele, Reference_allele
    ),
    # Compute the complementary Coded allele frequency
    Coded_allele_frequency = if_else(
      Beta < 0, 1 - Coded_allele_frequency, Coded_allele_frequency
    ),
    # Switch the sign of Beta
    Beta = if_else(Beta < 0, -Beta, Beta)
  ) %>%
  # Select all the columns
  dplyr::select(SNP:Position, Effect_allele_1, Reference_allele:R_2) %>%
  rename(Effect_allele = Effect_allele_1) %>%
  mutate(chr_pos = glue::glue("{Chromosome}_{Position}"))

# 3) Save the data -----------------------------------------------------
save(
  r2_0001,
  file = here::here("ligthart_crp", "data", "tidied_instruments.rda")
)
