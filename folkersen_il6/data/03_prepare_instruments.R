# Prepare instruments using Folkersen et al. (2020) data ---------------
library(tidyverse)
library(readxl)
library(openxlsx)
library(LDlinkR)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(rsnps)

load(here::here("folkersen_il6", "data", "iv_data.rda"))

# 1) Prepare SNPs data -------------------------------------------------
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
    EAF = if_else(Beta < 0, 1 - EAF, EAF),
    # Switch the sign of Beta
    Beta = if_else(Beta < 0, -Beta, Beta)
  ) %>%
  # Select all the columns
  dplyr::select(
    Chromosome:chr_pos, Effect_allele_1, Reference_allele:TotalSampleSize
  ) %>%
  dplyr::rename(Effect_allele = Effect_allele_1) %>%
  group_by(Chromosome) %>%
  arrange(P, .by_group = TRUE) %>%
  ungroup()

# 2) Get SNP number based on position ----------------------------------
ss <- SNPlocs.Hsapiens.dbSNP144.GRCh37

chr_1 <- snpsBySeqname(ss, "1")
pos <- as.double(tab_1[tab_1$Chromosome == 1, ]$Position)
idx <- match(pos, start(chr_1))
snp_chr_1 <- mcols(chr_1)$RefSNP_id[idx]
snps_ch1 <- tibble(
  Position = tab_1[tab_1$Chromosome == 1, ]$Position,
  SNP = snp_chr_1
) %>%
  add_column(Chromosome = 1, .before = 1)

# chr_6 <- snpsBySeqname(ss, "6")
# pos <- as.double(tab_1[tab_1$Chromosome == 6, ]$Position)
# idx <- match(pos, start(chr_6))
# snp_chr_6 <- mcols(chr_6)$RefSNP_id[idx]

# For chromosome 6 no position was found. The RS number will be set
# manually using dbSNP

snps_ch6 <- tibble(
  Chromosome = rep(6, length(tab_1[tab_1$Chromosome == 6, ]$Position)),
  Position = tab_1[tab_1$Chromosome == 6, ]$Position
) %>%
  mutate(rg = glue::glue("{Chromosome}:{Position}")) %>%
  # Manually add RS number
  mutate(
    SNP = case_when(
      rg == "6:32583159" ~ "rs4959106",
      rg == "6:32583682" ~ "rs34136174",
      rg == "6:32583426" ~ "rs35029150",
      rg == "6:32583557" ~ "rs9269909",
      rg == "6:32583677" ~ "rs36124427",
      rg == "6:32582627" ~ "rs3129752",
      rg == "6:32583046" ~ "rs3129754",
      rg == "6:32609094" ~ "rs1129737",
      rg == "6:32583299" ~ "rs34850435",
      rg == "6:32582406" ~ "rs5004279",
      rg == "6:32582401" ~ "rs5004278",
      rg == "6:32584355" ~ "rs1136758",
      rg == "6:32582932" ~ "rs9271340",
      rg == "6:32582548" ~ "rs9271328"
    )
  )

snps <- tab_1 %>%
  dplyr::left_join(
    bind_rows(snps_ch1, snps_ch6), by = c("Chromosome", "Position")
  ) %>%
  dplyr::select(SNP, everything()) %>%
  # Add the Gene
  dplyr::mutate(Gene = ncbi_snp_query(snps = .$SNP)$gene)

# 3) Perform LD clumping (r2<0.001)-------------------------------------
ld_ch1 <- SNPclip(
  snps = snps$SNP[snps$Chromosome == 1],
  pop = "CEU",
  r2_threshold = "0.001",
  maf_threshold = "0.01",
  token = Sys.getenv("LDLINK_TOKEN")
)

ch1 <- snps %>%
  filter(SNP %in% ld_ch1[ld_ch1$Details == "Variant kept.", 1])

ld_ch6 <- SNPclip(
  snps = snps$SNP[snps$Chromosome == 6],
  pop = "CEU",
  r2_threshold = "0.001",
  maf_threshold = "0.01",
  token = Sys.getenv("LDLINK_TOKEN")
)

ch6 <- snps %>%
  filter(SNP %in% ld_ch6[ld_ch6$Details == "Variant kept.", 1])

r2_0001 <- bind_rows(ch1, ch6) %>% dplyr::select(-rg)

# 3) Save the data -----------------------------------------------------
save(
  r2_0001,
  file = here::here("folkersen_il6", "data", "tidied_instruments.rda")
)
