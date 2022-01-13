# Import IVs from Ahluwalia for IL-6 -----------------------------------
library(tidyverse)
library(readxl)
library(janitor)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(rsnps)
library(openxlsx)

# The raw data obtained from Ahluwalia et al. can be accessed here:
# https://academic.oup.com/hmg/article/30/5/393/6124523#supplementary-data.
# The genetic associations are shown in Table S1, file
# "IL6_GWAS_Supp_Tables_HMG_11012021_ddab023".

# 1) Data import -------------------------------------------------------
ivs <- read_xlsx(
  here::here("ahluwalia_il6", "data-raw", "instruments_ahluwalia.xlsx")
) %>%
  dplyr::rename(
    SNP = rsID,
    Effect_allele = effect_allele,
    Reference_allele = other_allele,
    EAF = Effect_allele_frequency,
    Beta = beta
  )

# 2) Get the gene of the SNPs ------------------------------------------
id_snps <- ivs$SNP

dbsnps_query <- ncbi_snp_query(id_snps)

# 3) Get SNPs position using 37 build ----------------------------------
snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37

# The following SNPs were not found: rs8192284, rs8192282, rs660895,
# rs6910071, rs2395175. For them the position will be set manually
snps_query_37 <- snpsById(snps, ids = id_snps, ifnotfound = "drop")
ref_snps <- snps_query_37@elementMetadata@listData$RefSNP_id
pos_snps <- snps_query_37@ranges@pos %>%
  set_names(nm = ref_snps)

instr <- ivs %>%
  mutate(
    # Add position of build 37
    position_37 = if_else(
      SNP %in% ref_snps,
      unname(pos_snps[match(.data$SNP, names(pos_snps))]),
      NA_integer_
    ),
    # The position for 5 SNPs was not retrieved by the snpsById function.
    # For those SNP the position will be manually retrieved from dbSNP
    # (https://www.ncbi.nlm.nih.gov/snp/)
    position_37 = case_when(
      SNP == "rs8192284" ~ as.integer(154426970),
      SNP == "rs8192282" ~ as.integer(154401679),
      SNP == "rs660895" ~ as.integer(32577380),
      SNP == "rs6910071" ~ as.integer(32282854),
      SNP == "rs2395175" ~ as.integer(32405026),
      TRUE ~ position_37
    )
  ) %>%
  dplyr::select(-Position) %>%
  dplyr::rename(Position = position_37) %>%
  dplyr::select(
    SNP, Chromosome, Position,
    Effect_allele, Reference_allele, EAF, Beta, SE, P_value
  ) %>%
  dplyr::rename(Chromosome_x = Chromosome, Position_x = Position) %>%
  mutate(
    # Compute the F-statistic
    F_statistic = Beta^2/SE^2,
    # Compute the R^2
    R_2 = 2 * Beta^2 * EAF * (1 - EAF),
    # Add the gene
    Gene = dbsnps_query$gene
  )

# 4) Save instruments data ---------------------------------------------
save(instr, file = here::here("ahluwalia_il6", "data", "iv_data.rda"))
