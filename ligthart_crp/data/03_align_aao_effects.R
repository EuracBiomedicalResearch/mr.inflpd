# Harmonize CRP-AAO data for MR analysis -------------------------------
library(tidyverse)
library(readr)

# 1) Load the data -----------------------------------------------------
load(here::here("aao", "data", "aao2019_data.rda"))
load(here::here("ligthart_crp", "data", "tidied_instruments.rda"))

# 2) Prepare AAO data --------------------------------------------------
aao_data <- ma %>%
  mutate(
    Chromosome = parse_number(str_sub(MarkerName, 1, 5)),
    Position = sub(".*:", "", MarkerName),
    chr_pos = glue::glue("{Chromosome}_{Position}")
  ) %>%
  rename(
    Chromosome_y = Chromosome,
    Position_y = Position,
    A1_y = Allele1,
    A2_y = Allele2,
    freq_y = Freq1,
    beta_y = Effect,
    se_y = StdErr,
    p_y = `P-value`
  ) %>%
  dplyr::mutate(across(c(A1_y, A2_y), ~ str_to_upper(.))) %>%
  dplyr::select(
    chr_pos,
    Chromosome_y,
    Position_y,
    A1_y,
    A2_y,
    freq_y,
    beta_y,
    se_y,
    p_y
  )

# 3) Consistent column names -------------------------------------------
instr_mr <- r2_0001 %>%
  rename(
    Chromosome_x = Chromosome,
    Position_x = Position,
    A1_x = Effect_allele,
    A2_x = Reference_allele,
    freq_x = Coded_allele_frequency,
    beta_x = Beta,
    se_x = SE,
    p_x = P_value
  ) %>%
  dplyr::select(
    SNP,
    Closest_gene,
    chr_pos,
    Chromosome_x,
    Position_x,
    A1_x,
    A2_x,
    freq_x,
    beta_x,
    se_x,
    p_x
  )

# 4) Select instruments from PD data -----------------------------------
instr_nm <- instr_mr$chr_pos

aao_df <- aao_data %>%
  filter(chr_pos %in% instr_nm)

# Join CRP and PD data -------------------------------------------------
lj <- instr_mr %>%
  left_join(aao_df, by = "chr_pos")

# 5) Align PD estimates ------------------------------------------------
dd <- lj %>%
  # Here I must perform the check of the alleles. The pairs are:
  # A - G and C - T. If some of the pairs are not found, the
  # correspondence is A = T and C = G.
  mutate(
    A1_y = case_when(
      (A1_x == "C" & A2_x == "T") & (A1_y == "A" & A2_y == "C") ~ "T",
      (A1_x == "G" & A2_x == "A") & (A1_y == "T" & A2_y == "G") ~ "A",
      (A1_x == "A" & A2_x == "G") & (A1_y == "T" & A2_y == "G") ~ "A",
      (A1_x == "T" & A2_x == "C") & (A1_y == "A" & A2_y == "T") ~ "C",
      (A1_x == "A" & A2_x == "G") & (A1_y == "C" & A2_y == "G") ~ "A",
      TRUE ~ A1_y
    )
  ) %>%
  # Align the PD estimates given the exposure alleles
  mutate(
    beta_y = if_else((A1_x == A1_y), beta_y, -beta_y),
    freq_y = if_else((A1_x == A1_y), freq_y, 1 - freq_y)
  )

# 6) Prepare the data for the analysis ---------------------------------
df <- dd %>%
  dplyr::select(
    SNP, Closest_gene, Chromosome_x, Position_x, freq_x, freq_y,
    beta_x, se_x, p_x, beta_y, se_y, p_y, A1_x, A2_x, A1_y, A2_y
  )

# 7) Save the data -----------------------------------------------------
mr_data <- list(
  "instruments_CRP" = instr_mr,
  "instruments_AAO" = aao_df,
  "joined_instruments" = lj,
  "aligned_estimates" = dd,
  "mr_dataset" = df
)

save(
  mr_data,
  file = here::here("ligthart_crp", "data", "mr_data_aao.rda")
)

