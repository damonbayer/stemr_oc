`%notin%` <- Negate(`%in%`)

# to_human_scale ----------------------------------------------------------
to_human_scale = function(params_est) {
  c(R0 = exp(params_est[["R0_est"]]),
    dur_latent = exp(-params_est[["dur_latent_est"]]),
    dur_early = exp(-params_est[["dur_early_est"]]),
    dur_progress = exp(-params_est[["dur_progress_est"]]),
    ifr = expit(params_est[["ifr_est"]]),
    rho_death = expit(params_est[["rho_death_est"]]),
    phi_death = exp(-2 * params_est[["phi_death_est"]]),
    alpha0 = exp(params_est[["alpha0_est"]]),
    alpha1 = expit(params_est[["alpha1_est"]]),
    kappa = exp(-2 * params_est[["kappa_est"]]))
}


# extract_stem_parameter_posterior ----------------------------------------
extract_stem_parameter_posterior <- function(multi_chain_stem_fit, transform = "estimation") {
  map(multi_chain_stem_fit$stem_fit_list,
      function(stem_fit) {
        if (is.character(transform) && transform == "estimation") {
          parameter_samples <- stem_fit$results$posterior$parameter_samples_est
        } else if (is.character(transform) && transform == "natural") {
          parameter_samples <- stem_fit$results$posterior$parameter_samples_nat
        } else {
          parameter_samples <- t(apply(stem_fit$results$posterior$parameter_samples_est, 1, transform))
        }
        mcmc(parameter_samples, start = 100, end = 200000, thin = 100)}) %>%
    as.mcmc.list()
}


# extract_epi_curves ------------------------------------------------------
# Probably a faster way to write this
extract_epi_curves <- function(multi_chain_stem_fit, curve_type = "prevalence", tidy = F) {
  if (curve_type == "prevalence") curve_type <- "p"
  if (curve_type == "incidence") curve_type <- "i"
  if (curve_type %notin% c("i", "p")) stop('curve_type must be one of "i", "incidence", "p", or "prevalence"')
  if (!is.logical(tidy)) stop('tidy myust be TRUE or FALSE')

  epi_curves <- imap_dfr(multi_chain_stem_fit$stem_fit_list, function(stem_fit, chain) {

    map_dfr(1:multi_chain_stem_fit$n_iterations, function(i) {
      if (curve_type == "p") {
        if (stem_fit$dynamics$fixed_inits) {
          init_state <- stem_fit$dynamics$initdist_params
        } else {
          init_state <- stem_fit$result$posterior$initdist_samples[i,]
        }
        path <- incidence2prevalence(path = stem_fit$result$posterior$latent_paths[,,i],
                                     flow_matrix = stem_fit$dynamics$flow_matrix_ode,
                                     init_state = init_state)
      } else {
        path <- stem_fit$result$posterior$latent_paths[,,i]
      }

      path %>%
        as_tibble() %>%
        mutate(.iteration = i)
    }) %>%
      mutate(.chain = chain)
  }) %>%
    group_by(time) %>%
    mutate(.draw = row_number()) %>%
    ungroup()

  if (tidy == T) {
    epi_curves <- epi_curves %>%
      pivot_longer(-c(time, .iteration, .chain, .draw)) %>%
      mutate(name = fct_inorder(name))
  }
  epi_curves
}


# plot_epi_curves ---------------------------------------------------------
# possible add an option to separate by chain
plot_epi_curves <- function(multi_chain_stem_fit, curve_type = "p") {
  get_epi_curves(multi_chain_stem_fit = multi_chain_stem_fit, curve_type = curve_type, tidy = T) %>%
    ggplot(aes(time, value)) +
    stat_lineribbon() +
    facet_wrap(. ~ name, scales = "free_y") +
    scale_fill_brewer()
}
