# Computing I^2 statistic ----------------------------------------------
i2_fun <- function(mr_df) {

  assertive::assert_is_data.frame(mr_df)

  # Wald estimates
  est <- mr_df$beta_y/mr_df$beta_x
  se <- abs(mr_df$se_y/mr_df$beta_x)

  ma <- meta::metagen(est, se)

  est <- round(ma$I2 * 100, 2L)
  lower <- round(ma$lower.I2 * 100, 2L)
  upper <- round(ma$upper.I2 * 100, 2L)
  glue::glue("{est} [{lower}; {upper}]")

}

# Function for Q-statistics for IVs ------------------------------------
q_stat_fun <- function(data) {

  assertive::assert_is_data.frame(data)

  # Retrieve object for computing Q-stat
  b_x <- data$beta_x
  b_y <- data$beta_y

  se_x <- data$se_x
  se_y <- data$se_y

  ratio <- b_y/b_x
  se_ratio <- abs(se_y/b_x)

  w <- 1/se_ratio^2  # 1st order weights

  # Compute Q-stat
  mu_hat <- sum(ratio * w)/sum(w)
  q_ivw <- w * (ratio - mu_hat)^2

  # Compute thresholds for plotting
  l1 <- qchisq(1 - 0.05, df = 1)
  l2 <- qchisq(1 - 0.01, df = 1)
  l3 <- qchisq(1 - 0.0019, df = 1)

  # Create the dataframe for the plot
  plot_df <- tibble::tibble(
    snps = data$SNP,
    q_stat = q_ivw
  )

  # Plot Q-statistics
  q_plot <- ggplot2::ggplot(
    data = plot_df,
    mapping = aes(x = snps, y = q_stat)
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(
      aes(yintercept = l1, linetype = "5th percentile"),
      colour = "forestgreen", size = 1
    ) +
    ggplot2::geom_hline(
      aes(yintercept = l2, linetype = "1st percentile"),
      colour = "dodgerblue2", size = 1
    ) +
    ggplot2::geom_hline(
      aes(yintercept = l3, linetype = "0.19th percentile"),
      colour = "firebrick", size = 1
    ) +
    ggplot2::scale_linetype_manual(
      name = "Thresholds",
      values = c("dashed", "dashed", "dashed"),
      guide = guide_legend(
        override.aes = list(
          color = c("firebrick", "dodgerblue2", "forestgreen")
        )
      )
    ) +
    ggplot2::xlab("SNP") +
    ggplot2::ylab("Q-statistic") +
    ggplot2::theme_light() +
    ggplot2::theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  # Delete IVs above 0.19th percentile and 5th percentile
  df_trim_l3 <- data %>%
    dplyr::mutate(q_ivw = q_ivw) %>%
    dplyr::filter(q_ivw < l3) %>%
    dplyr::select(-q_ivw)

  df_trim_l1 <- data %>%
    dplyr::mutate(q_ivw = q_ivw) %>%
    dplyr::filter(q_ivw < l1) %>%
    dplyr::select(-q_ivw)

  # Output within a list
  list("plot" = q_plot, "data_l3" = df_trim_l3, "data_l1" = df_trim_l1)

}

# Function for Radial Plot method --------------------------------------
radial_fun <- function(data) {

  assertive::assert_is_data.frame(data)

  # Prepare the data for radial MR
  fr <- RadialMR::format_radial(
    BXG = data$beta_x,
    BYG = data$beta_y,
    seBXG = data$se_x,
    seBYG = data$se_x,
    RSID = data$SNP
  )

  # Select outliers identified with radial method. Moreover, retrieve
  # a dataframe with IVs excluded by the radial plot method
  mr_r <- MRPracticals::ivw_radial(fr)
  fr_to_remove <- mr_r$data$Outliers

  df_incl <- data %>%
    dplyr::mutate(outlier = fr_to_remove) %>%
    dplyr::filter(outlier == "Variant") %>%
    dplyr::select(-outlier)

  df_excl <- data %>%
    dplyr::mutate(outlier = fr_to_remove) %>%
    dplyr::filter(outlier != "Variant") %>%
    dplyr::select(-outlier)

  # Plot the Radial results
  pr <- MRPracticals::plot_radial(
    r_object = mr_r,
    radial_scale = TRUE,
    show_outliers = TRUE,
    scale_match = TRUE
  )

  # List with results of the function
  list("data_incl" = df_incl, "data_excl" = df_excl, "radial_plot" = pr)

}

# Main function for MR -------------------------------------------------
mr_fun <- function(data, input) {

  assertive::assert_is_data.frame(data)
  if(class(input) != "MRInput") {
    usethis::ui_stop("'input' must be of class 'MRInput'")
  }

  # IVW ----------------------------------------------------------------
  ivw_fit <- MendelianRandomization::mr_ivw(
    input, model = "random", robust = TRUE, penalized = TRUE
  )

  ivw_het <- MendelianRandomization::mr_ivw(input, model = "fixed")

  ivw_res <- tibble(
    m = "IVW-RE-Robust",
    type = "Estimate",
    est = ivw_fit@Estimate,
    se = ivw_fit@StdError,
    lower = ivw_fit@CILower,
    upper = ivw_fit@CIUpper,
    pval = ivw_fit@Pvalue
  )

  # MR-Egger -----------------------------------------------------------
  egger_fit <- MendelianRandomization::mr_egger(input)

  egger_res <- tibble(
    m = "MR-Egger",
    type = "Estimate",
    est = egger_fit@Estimate,
    se = egger_fit@StdError.Est,
    lower = egger_fit@CILower.Est,
    upper = egger_fit@CIUpper.Est,
    pval = egger_fit@Pvalue.Est
  )

  # Median -------------------------------------------------------------
  median_fit <- MendelianRandomization::mr_median(
    input, weighting = "weighted"
  )

  median_res <- tibble(
    m = "Weighted-Median",
    type = "Estimate",
    est = median_fit@Estimate,
    se = median_fit@StdError,
    lower = median_fit@CILower,
    upper = median_fit@CIUpper,
    pval = median_fit@Pvalue
  )

  # Mode-based ---------------------------------------------------------
  mbe_fit <- MendelianRandomization::mr_mbe(
    input, weighting = "weighted"
  )

  mbe_res <- tibble(
    m = "Weighted-Mode",
    type = "Estimate",
    est = mbe_fit@Estimate,
    se = mbe_fit@StdError,
    lower = mbe_fit@CILower,
    upper = mbe_fit@CIUpper,
    pval = mbe_fit@Pvalue
  )

  # MR-PRESSO ----------------------------------------------------------
  presso <- MRPRESSO::mr_presso(
    BetaOutcome = "beta_y",
    BetaExposure = "beta_x",
    SdOutcome = "se_y",
    SdExposure = "se_x",
    data = as.data.frame(data),
    OUTLIERtest = TRUE,
    DISTORTIONtest = TRUE
  )

  presso_res <- tibble::tibble(
    m = rep("MR PRESSO", 2L),
    method = c("Raw", "Outlier-corrected"),
    type = rep("Estimate", 2L),
    est = presso$`Main MR results`$`Causal Estimate`,
    se = presso$`Main MR results`$Sd,
    lower = presso$`Main MR results`$`Causal Estimate` - 1.96 *
      presso$`Main MR results`$Sd,
    upper = presso$`Main MR results`$`Causal Estimate` + 1.96 *
      presso$`Main MR results`$Sd,
    pval = presso$`Main MR results`$`P-value`
  )

  # Results into a final list ------------------------------------------
  list(
    "ivw" = list(
      "model" = ivw_fit, "estimate" = ivw_res,
      "het" = ivw_het
    ),
    "egger" = list("model" = egger_fit, "estimate" = egger_res),
    "mr_median" = list("model" = median_fit, "estimate" = median_res),
    "mr_mode" = list("model" = mbe_fit, "estimate" = mbe_res),
    "presso" = list("model" = presso, "estimate" = presso_res)
  )

}
