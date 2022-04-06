# Harmonized CRP-AAO data for MR analysis ------------------------------
library(tidyverse)
library(readr)
library(TwoSampleMR)

# 1) Load the data -----------------------------------------------------
load(here::here("aao", "data", "aao2019_data.rda"))
load(here::here("ccgc_crp", "data", "tidied_instruments.rda"))

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
instr_mr_r2_0001 <- r2_0001 %>%
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

instr_mr_r2_001 <- r2_001 %>%
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

instr_mr_r2_01 <- r2_01 %>%
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

# 4) Select instruments from AAO data ----------------------------------
instr_nm_r2_0001 <- instr_mr_r2_0001$chr_pos
instr_nm_r2_001 <- instr_mr_r2_001$chr_pos
instr_nm_r2_01 <- instr_mr_r2_01$chr_pos

aao_df_r2_0001 <- aao_data %>%
  filter(chr_pos %in% instr_nm_r2_0001)

aao_df_r2_001 <- aao_data %>%
  filter(chr_pos %in% instr_nm_r2_001)

aao_df_r2_01 <- aao_data %>%
  filter(chr_pos %in% instr_nm_r2_01)

# Join CRP and AAO data ------------------------------------------------
lj_r2_0001 <- instr_mr_r2_0001 %>%
  left_join(aao_df_r2_0001, by = "chr_pos")

lj_r2_001 <- instr_mr_r2_001 %>%
  left_join(aao_df_r2_001, by = "chr_pos")

lj_r2_01 <- instr_mr_r2_01 %>%
  left_join(aao_df_r2_01, by = "chr_pos")

# 5) Align PD estimates ------------------------------------------------
dd_r2_0001 <- lj_r2_0001 %>%
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

dd_r2_001 <- lj_r2_001 %>%
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

dd_r2_01 <- lj_r2_01 %>%
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
df_r2_0001 <- dd_r2_0001 %>%
  dplyr::select(
    SNP, Closest_gene, Chromosome_x, Position_x, freq_x, freq_y,
    beta_x, se_x, p_x, beta_y, se_y, p_y, A1_x, A2_x, A1_y, A2_y
  )

df_r2_001 <- dd_r2_001 %>%
  dplyr::select(
    SNP, Closest_gene, Chromosome_x, Position_x, freq_x, freq_y,
    beta_x, se_x, p_x, beta_y, se_y, p_y, A1_x, A2_x, A1_y, A2_y
  )

df_r2_01 <- dd_r2_01 %>%
  dplyr::select(
    SNP, Closest_gene, Chromosome_x, Position_x, freq_x, freq_y,
    beta_x, se_x, p_x, beta_y, se_y, p_y, A1_x, A2_x, A1_y, A2_y
  )

# 7) Get LD matrix when R2<0.1 -----------------------------------------
ld_r2_01 <- ld_matrix(df_r2_01$SNP, with_alleles = FALSE)

# 8) Save the data -----------------------------------------------------
mr_data_r2_0001 <- list(
  "instruments_CRP" = instr_mr_r2_0001,
  "instruments_AAO" = aao_df_r2_0001,
  "joined_instruments" = lj_r2_0001,
  "aligned_estimates" = dd_r2_0001,
  "mr_dataset" = df_r2_0001
)

mr_data_r2_001 <- list(
  "instruments_CRP" = instr_mr_r2_001,
  "instruments_AAO" = aao_df_r2_001,
  "joined_instruments" = lj_r2_001,
  "aligned_estimates" = dd_r2_001,
  "mr_dataset" = df_r2_001
)

mr_data_r2_01 <- list(
  "instruments_CRP" = instr_mr_r2_01,
  "instruments_AAO" = aao_df_r2_01,
  "joined_instruments" = lj_r2_01,
  "aligned_estimates" = dd_r2_01,
  "mr_dataset" = df_r2_01,
  "ld_matrix" = ld_r2_01
)

save(
  mr_data_r2_0001, mr_data_r2_001, mr_data_r2_01,
  file = here::here("ccgc_crp", "data", "mr_data_aao.rda")
)
