# Harmonized CRP-PD data for MR analysis -------------------------------
library(tidyverse)
library(readr)

# 1) Load the data -----------------------------------------------------
load(here::here("pd", "data", "nalls2019_data.rda"))
load(here::here("ccgc_crp", "data", "tidied_instruments.rda"))

# 2) Prepare PD data ---------------------------------------------------
pd_data <- ma %>%
  mutate(
    Chromosome = parse_number(str_sub(SNP, 1, 5)),
    Position = sub(".*:", "", SNP),
    chr_pos = glue::glue("{Chromosome}_{Position}")
  ) %>%
  rename(
    Chromosome_y = Chromosome,
    Position_y = Position,
    A1_y = A1,
    A2_y = A2,
    freq_y = freq,
    beta_y = b,
    se_y = se,
    p_y = p
  ) %>%
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
    freq_x = Effect_allele_frequency,
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

pd_df <- pd_data %>%
  filter(chr_pos %in% instr_nm)

# Join CRP and PD data -------------------------------------------------
lj <- instr_mr %>%
  left_join(pd_df, by = "chr_pos")

# 5) Align PD estimates ------------------------------------------------
dd <- lj %>%
  # Here I must perform the check of the alleles. The pairs are:
  # A - G and C - T. If some of the pairs are not found, the
  # correspondence is A = T and C = G.
  mutate(
    A1_y = case_when(
      (A1_x == "T" & A2_x == "C") & (A1_y == "A" & A2_y == "G") ~ "T",
      (A1_x == "G" & A2_x == "A") & (A1_y == "C" & A2_y == "G") ~ "G",
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
  "instruments_PD" = pd_df,
  "joined_instruments" = lj,
  "aligned_estimates" = dd,
  "mr_dataset" = df
)

save(
  mr_data,
  file = here::here("ccgc_crp", "data", "mr_data_pd.rda")
)
