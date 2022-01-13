find_proxy <- function(snp, exp_data, out_data) {

  assertive::assert_is_character(snp)
  assertive::assert_is_data.frame(exp_data)
  assertive::assert_is_data.frame(out_data)

  # # Retrieve the chromosome and position of the SNP
  # cp <- exp_data %>%
  #   dplyr::filter(SNP == snp) %>%
  #   .[["chr_pos"]]

  if(!snp %in% exp_data$SNP) {
    usethis::ui_stop("'snp' must be contained in the exposure data")
  }

  # The safe version of LDproxy is created for not biallelic variant.
  # Those variants will not be used in the analysis
  safe_LDproxy <- purrr::safely(
    LDlinkR::LDproxy, otherwise = NA_real_
  )

  # Find the proxy
  query <- safe_LDproxy(
    snp = snp,
    pop = "CEU",
    r2d = "r2",
    token = Sys.getenv("LDLINK_TOKEN")
  )

  if (
    query$result[1, 1] == "  error: rs645692 is not a biallelic variant.,"
  ) {

    res <- exp_data %>%
      dplyr::filter(SNP == snp)


  } else {

    # Retrieve the seqname of the proxy
    proxy_seqname <- query$result %>%
      # Take only the SNPs which are present in the outcome data
      dplyr::filter(Coord %in% out_data$seqname) %>%
      # Remove the SNPs that are already present in the exposure data
      dplyr::filter(!RS_Number %in% exp_data$SNP) %>%
      # Take as proxy the one with the highest r^2 and the minimum
      # distance, which is in the first row
      dplyr::slice(1) %>%
      dplyr::select(Coord, R2)

    # Take the proxy estimate from the outcome dataset
    res <- out_data %>%
      dplyr::filter(seqname == proxy_seqname$Coord) %>%
      dplyr::mutate(r2 = proxy_seqname$R2) %>%
      dplyr::bind_cols(
        exp_data %>%
          dplyr::filter(SNP == snp) %>%
          dplyr::select(SNP:R_2)
      ) %>%
      dplyr::rename(chr_pos = chr_pos...2) %>%
      dplyr::select(-chr_pos...13)

  }

  message(glue::glue("{snp} is finished!"))

  res

}

