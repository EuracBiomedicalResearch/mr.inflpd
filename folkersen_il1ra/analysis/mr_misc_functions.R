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

  # Compute Bonferroni threshold
  alpha <- 0.05/length(q_ivw)
  q_thr <- stats::qchisq(1 - alpha, df = 1)

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
    ggplot2::geom_hline(aes(yintercept = q_thr), linetype = "dashed") +
    ggplot2::xlab("SNP") +
    ggplot2::ylab("Q-statistic") +
    ggplot2::theme_light() +
    ggplot2::theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  # Delete IVs above 0.19th percentile and 5th percentile
  df_trim <- data %>%
    dplyr::mutate(q_ivw = q_ivw) %>%
    dplyr::filter(q_ivw < q_thr) %>%
    dplyr::select(-q_ivw)

  # Output within a list
  list("plot" = q_plot, "data_trim" = df_trim)


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
