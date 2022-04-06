# Prepare instruments using Folkersen et al. (2020) data ---------------
library(tidyverse)
library(readxl)
library(openxlsx)
library(LDlinkR)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(rsnps)

load(here::here("folkersen_il1ra", "data", "iv_data.rda"))

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

chr_2 <- snpsBySeqname(ss, "2")
pos <- as.double(tab_1[tab_1$Chromosome == 2, ]$Position)
idx <- match(pos, start(chr_2))
snp_chr_2 <- mcols(chr_2)$RefSNP_id[idx]
snps_ch2 <- tibble(
  Position = tab_1[tab_1$Chromosome == 2, ]$Position,
  SNP = snp_chr_2
) %>%
  add_column(Chromosome = 2, .before = 1)

chr_17 <- snpsBySeqname(ss, "17")
pos <- as.double(tab_1[tab_1$Chromosome == 17, ]$Position)
idx <- match(pos, start(chr_17))
snp_chr_17 <- mcols(chr_17)$RefSNP_id[idx]
snps_ch17 <- tibble(
  Position = tab_1[tab_1$Chromosome == 17, ]$Position,
  SNP = snp_chr_17
) %>%
  add_column(Chromosome = 17, .before = 1)

chr_19 <- snpsBySeqname(ss, "19")
pos <- as.double(tab_1[tab_1$Chromosome == 19, ]$Position)
idx <- match(pos, start(chr_19))
snp_chr_19 <- mcols(chr_19)$RefSNP_id[idx]
snps_ch19 <- tibble(
  Position = tab_1[tab_1$Chromosome == 19, ]$Position,
  SNP = snp_chr_19
) %>%
  add_column(Chromosome = 19, .before = 1)

snps <- tab_1 %>%
  dplyr::left_join(
    bind_rows(snps_ch2, snps_ch17, snps_ch19),
    by = c("Chromosome", "Position")
  ) %>%
  # For now remove the SNPs that were not found in the assembly 37
  filter(!is.na(SNP)) %>%
  dplyr::select(SNP, everything()) %>%
  # Add the gene
  mutate(Gene = ncbi_snp_query(snps = .$SNP)$gene)

# 3) LD clumping -------------------------------------------------------
# R2<0.001 -------------------------------------------------------------
ld_ch2 <- SNPclip(
  snps = snps$SNP[snps$Chromosome == 2],
  pop = "CEU",
  r2_threshold = "0.001",
  maf_threshold = "0.01",
  token = Sys.getenv("LDLINK_TOKEN")
)

ch2 <- snps %>%
  filter(SNP %in% ld_ch2[ld_ch2$Details == "Variant kept.", 1])

ld_ch17 <- SNPclip(
  snps = snps$SNP[snps$Chromosome == 17],
  pop = "CEU",
  r2_threshold = "0.001",
  maf_threshold = "0.01",
  token = Sys.getenv("LDLINK_TOKEN")
)

ch17 <- snps %>%
  filter(SNP %in% ld_ch17[ld_ch17$Details == "Variant kept.", 1])

ld_ch19 <- SNPclip(
  snps = snps$SNP[snps$Chromosome == 19],
  pop = "CEU",
  r2_threshold = "0.001",
  maf_threshold = "0.01",
  token = Sys.getenv("LDLINK_TOKEN")
)

ch19 <- snps %>%
  filter(SNP %in% ld_ch19[ld_ch19$Details == "Variant kept.", 1])

r2_0001 <- bind_rows(ch2, ch17, ch19)

# R2<0.01 --------------------------------------------------------------
ld_ch2 <- SNPclip(
  snps = snps$SNP[snps$Chromosome == 2],
  pop = "CEU",
  r2_threshold = "0.01",
  maf_threshold = "0.01",
  token = Sys.getenv("LDLINK_TOKEN")
)

ch2 <- snps %>%
  filter(SNP %in% ld_ch2[ld_ch2$Details == "Variant kept.", 1])

ld_ch17 <- SNPclip(
  snps = snps$SNP[snps$Chromosome == 17],
  pop = "CEU",
  r2_threshold = "0.01",
  maf_threshold = "0.01",
  token = Sys.getenv("LDLINK_TOKEN")
)

ch17 <- snps %>%
  filter(SNP %in% ld_ch17[ld_ch17$Details == "Variant kept.", 1])

ld_ch19 <- SNPclip(
  snps = snps$SNP[snps$Chromosome == 19],
  pop = "CEU",
  r2_threshold = "0.01",
  maf_threshold = "0.01",
  token = Sys.getenv("LDLINK_TOKEN")
)

ch19 <- snps %>%
  filter(SNP %in% ld_ch19[ld_ch19$Details == "Variant kept.", 1])

r2_001 <- bind_rows(ch2, ch17, ch19)

# R2<0.1 --------------------------------------------------------------
ld_ch2 <- SNPclip(
  snps = snps$SNP[snps$Chromosome == 2],
  pop = "CEU",
  r2_threshold = "0.1",
  maf_threshold = "0.01",
  token = Sys.getenv("LDLINK_TOKEN")
)

ch2 <- snps %>%
  filter(SNP %in% ld_ch2[ld_ch2$Details == "Variant kept.", 1])

ld_ch17 <- SNPclip(
  snps = snps$SNP[snps$Chromosome == 17],
  pop = "CEU",
  r2_threshold = "0.1",
  maf_threshold = "0.01",
  token = Sys.getenv("LDLINK_TOKEN")
)

ch17 <- snps %>%
  filter(SNP %in% ld_ch17[ld_ch17$Details == "Variant kept.", 1])

ld_ch19 <- SNPclip(
  snps = snps$SNP[snps$Chromosome == 19],
  pop = "CEU",
  r2_threshold = "0.1",
  maf_threshold = "0.01",
  token = Sys.getenv("LDLINK_TOKEN")
)

ch19 <- snps %>%
  filter(SNP %in% ld_ch19[ld_ch19$Details == "Variant kept.", 1])

r2_01 <- bind_rows(ch2, ch17, ch19)

# 3) Save into the data ------------------------------------------------
save(
  r2_0001, r2_001, r2_01,
  file = here::here(
    "folkersen_il1ra", "data", "tidied_instruments.rda"
  )
)
