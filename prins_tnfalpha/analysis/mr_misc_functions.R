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
  # All the combinations of the inputs will be checked
  list_inp <- tidyr::crossing(
    model = c("fixed", "random"),
    robust = c(FALSE, TRUE),
    penalized = c(FALSE, TRUE),
    weights = c("simple", "delta")
  ) %>%
    tibble::add_column(object = map(1:nrow(.), ~ input), .before = 1)

  method <- ifelse(list_inp$model == "fixed", "Fixed", "Random")
  robust <- ifelse(list_inp$robust, "Robust", "")
  pen <- ifelse(list_inp$penalized, "Penalized", "")
  wts <- ifelse(list_inp$weights == "simple", "1st order", "2nd order")

  nm <- stringr::str_trim(glue::glue("{method} {wts} {robust} {pen}"))

  ivw <- purrr::pmap(list_inp, MendelianRandomization::mr_ivw) %>%
    purrr::set_names(x = ., nm = nm)

  ivw_res <- ivw %>%
    purrr::imap_dfr(
      .x = .,
      ~ tibble(
        m = "IVW",
        method = .y,
        type = "Estimate",
        est = .x@Estimate,
        se = .x@StdError,
        lower = .x@CILower,
        upper = .x@CIUpper,
        pval = .x@Pvalue
      )
    )

  # MR-Egger -----------------------------------------------------------
  # All the combinations of the inputs will be checked
  list_inp <- tidyr::crossing(
    robust = c(FALSE, TRUE),
    penalized = c(FALSE, TRUE),
  ) %>%
    tibble::add_column(object = map(1:nrow(.), ~ input), .before = 1)

  robust <- ifelse(list_inp$robust, "Robust", "")
  pen <- ifelse(list_inp$penalized, "Penalized", "")

  nm <- stringr::str_trim(glue::glue("{robust} {pen}"))

  egger <- purrr::pmap(list_inp, MendelianRandomization::mr_egger) %>%
    purrr::set_names(x = ., nm = nm)

  egger_res <- egger %>%
    purrr::imap_dfr(
      .x = .,
      ~ tibble(
        m = rep("MR Egger", 2L),
        method = rep(.y, 2L),
        type = c("Estimate", "Intercept"),
        est = c(.x@Estimate, .x@Intercept),
        se = c(.x@StdError.Est, .x@StdError.Int),
        lower = c(.x@CILower.Est, .x@CILower.Int),
        upper = c(.x@CIUpper.Est, .x@CIUpper.Int),
        pval = c(.x@Pvalue.Est, .x@Pvalue.Int)
      )
    )

  # Median -------------------------------------------------------------
  # All the combinations of the inputs will be checked
  list_inp <- tibble::tibble(
    weighting = c("simple", "weighted", "penalized")
  ) %>%
    tibble::add_column(object = map(1:nrow(.), ~ input), .before = 1)

  nm <- c("Simple", "Weighted", "Penalized")

  med <- purrr::pmap(list_inp, MendelianRandomization::mr_median) %>%
    purrr::set_names(x = ., nm = nm)

  med_res <- med %>%
    purrr::imap_dfr(
      .x = .,
      ~ tibble(
        m = "MR Median",
        method = .y,
        type = "Estimate",
        est = .x@Estimate,
        se = .x@StdError,
        lower = .x@CILower,
        upper = .x@CIUpper,
        pval = .x@Pvalue
      )
    )

  # Mode-based ---------------------------------------------------------
  # Both unweighted and weighted estimators will be used
  list_inp <- tibble::tibble(
    weighting = c("unweighted", "weighted")
  ) %>%
    tibble::add_column(object = map(1:nrow(.), ~ input), .before = 1)

  nm <- c("Unweighted", "Weighted")

  mbe <- purrr::pmap(list_inp, MendelianRandomization::mr_mbe) %>%
    purrr::set_names(x = ., nm = nm)

  mbe_res <- mbe %>%
    purrr::imap_dfr(
      .x = .,
      ~ tibble(
        m = "MR Mode",
        method = .y,
        type = "Estimate",
        est = .x@Estimate,
        se = .x@StdError,
        lower = .x@CILower,
        upper = .x@CIUpper,
        pval = .x@Pvalue
      )
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

  # Lasso --------------------------------------------------------------
  # Default parameters values are used: lambda parameter calculated
  # using heterogeneity stopping rule

  # If MR Lasso throws an error use the safe function that returns NA
  # instead of stopping
  safe_mr_lasso <- purrr::safely(
    MendelianRandomization::mr_lasso, otherwise = NA_real_
  )

  lasso_model <- safe_mr_lasso(object = inp)

  if (class(lasso_model$result) != "MRLasso") {

    lasso_res <- tibble::tibble(
      m = "MR Lasso",
      method = "Het stopping rule",
      type = "Estimate",
      est = NA_real_,
      se = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      pval = NA_real_
    )

  } else {

    lasso_res <- tibble::tibble(
      m = "MR Lasso",
      method = "Het stopping rule",
      type = "Estimate",
      est = lasso_model$result@Estimate,
      se = lasso_model$result@StdError,
      lower = lasso_model$result@CILower,
      upper = lasso_model$result@CIUpper,
      pval = lasso_model$result@Pvalue
    )

  }

  # MR-RAPS ------------------------------------------------------------
  # Default parameters (overdispersion = TRUE, loss.function = Tukey)
  raps <- mr.raps::mr.raps(
    b_exp = df$beta_x,
    b_out = df$beta_y,
    se_exp = df$se_x,
    se_out = df$se_y,
    over.dispersion = TRUE,
    loss.function = "tukey"
  )

  raps_res <- tibble::tibble(
    m = "MR RAPS",
    method = "Overdispersion Tukey",
    type = "Estimate",
    est = raps$b,
    se = raps$se,
    lower = raps$b - 1.96 * raps$se,
    upper = raps$b + 1.96 * raps$se,
    pval = raps$pval
  )

  # Maximum-Likelihood -------------------------------------------------
  # Both fixed and random effects will be used
  list_inp <- tibble::tibble(
    model = c("fixed", "random")
  ) %>%
    tibble::add_column(object = map(1:nrow(.), ~ input), .before = 1)

  nm <- c("Fixed", "Random")

  maxlik <- purrr::pmap(list_inp, MendelianRandomization::mr_maxlik) %>%
    purrr::set_names(x = ., nm = nm)

  maxlik_res <- mbe %>%
    purrr::imap_dfr(
      .x = .,
      ~ tibble(
        m = "MR MaxLik",
        method = .y,
        type = "Estimate",
        est = .x@Estimate,
        se = .x@StdError,
        lower = .x@CILower,
        upper = .x@CIUpper,
        pval = .x@Pvalue
      )
    )

  # Contamination mixture ----------------------------------------------
  # Default values will be used. However, a deeper exploration of the
  # values should be conducted, since default values may not be
  # suitable for many applications (as suggested by the authors)
  conmix <- MendelianRandomization::mr_conmix(inp)

  conmix_res <- tibble::tibble(
    m = "MR ConMix",
    method = "Default",
    type = "Estimate",
    est = conmix@Estimate,
    se = NA_real_,
    lower = conmix@CILower,
    upper = conmix@CIUpper,
    pval = conmix@Pvalue
  )

  # MR-Mixture ---------------------------------------------------------
  # Default parameters values will be used except for the grid over
  # theta values (from -2 to 2 by 0.01)
  mix <- MRMix::MRMix(
    betahat_x = df$beta_x,
    betahat_y = df$beta_y,
    sx = df$se_x,
    sy = df$se_y,
    theta_temp_vec = seq(from = -2, to = 2, by = 0.01)
  )

  mix_res <- tibble::tibble(
    m = "MR Mixture",
    method = "Default",
    type = "Estimate",
    est = mix$theta,
    se = mix$SE_theta,
    lower = mix$theta - 1.96 * mix$SE_theta,
    upper = mix$theta + 1.96 * mix$SE_theta,
    pval = mix$pvalue_theta
  )

  # Results into a final list ------------------------------------------
  list(
    "ivw" = list("models" = ivw, "estimate" = ivw_res),
    "egger" = list("models" = egger, "estimate" = egger_res),
    "mr_median" = list("models" = med, "estimate" = med_res),
    "mr_mode" = list("models" = mbe, "estimate" = mbe_res),
    "presso" = list("models" = presso, "estimate" = presso_res),
    "lasso" = list("models" = lasso_model, "estimate" = lasso_res),
    "raps" = list("models" = raps, "estimate" = raps_res),
    "conmix" = list("models" = conmix, "estimate" = conmix_res),
    "mix" = list("models" = mix, "estimate" = mix_res)
  )

}
