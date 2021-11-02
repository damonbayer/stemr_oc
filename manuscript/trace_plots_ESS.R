library(tidyverse)
library(tidybayes)
library(scales)
library(cowplot)
library(fs)
library(stemr)
library(coda)
library(MCMCvis)
library(xtable)
source('stemr_functions.R')
theme_set(cowplot::theme_minimal_hgrid())

main_results <- read_rds("final_models/results/2021-02-28.rds")
simulated_results <- read_rds("final_models/results/simulated_data.rds")

extract_log_posterior <- function(multi_chain_stem_fit) {
  names_to_extract <-
    c("data_log_lik",
      "params_log_prior",
      "initdist_log_lik",
      "tparam_log_lik") %>%
    set_names(., .)

  posterior_list <-
    multi_chain_stem_fit$stem_fit_list %>%
    map("results") %>%
    map("posterior")

  map(names_to_extract, ~map(posterior_list, .)) %>%
    enframe() %>%
    unnest(value) %>%
    group_by(name) %>%
    mutate(.chain = row_number()) %>%
    unnest(value) %>%
    group_by(name, .chain) %>%
    mutate(.iteration = row_number()) %>%
    ungroup() %>%
    pivot_wider(-value) %>%
    mutate(log_lik = data_log_lik + initdist_log_lik + tparam_log_lik) %>%
    mutate(log_post = log_lik + params_log_prior)
}

extract_mcmc_list <- function(multi_chain_stem_fit) {
  map(multi_chain_stem_fit$stem_fit_list, ~{
    initdist_draws <-
      .$results$posterior$initdist_draws[[1]] %>%
      t() %>%
      `colnames<-`(
        c("$\\frac{\\exp(S_0)}{\\exp(S_0) + 1}$",
          "$\\frac{\\exp(\\tilde{I}_0)}{\\exp(\\tilde{I}_0) + 1}$",
          "$\\frac{\\exp(\\tilde{R}_0)}{\\exp(\\tilde{R}_0) + 1}$",
          "$\\frac{\\exp(\\tilde{D}_0)}{\\exp(\\tilde{D}_0) + 1}$"))
    parameter_samples <- .$results$posterior$parameter_samples_est

    mcmc(cbind(parameter_samples, initdist_draws))
  }) %>%
    as.mcmc.list()
}

main_log_lik <- extract_log_posterior(main_results)

simulated_log_lik <- extract_log_posterior(simulated_results)

main_log_lik_plot <-
  main_log_lik %>%
  rename(Chain = .chain) %>%
  ggplot(aes(.iteration, log_lik)) +
  facet_wrap(. ~ Chain,
             labeller = label_both,
             nrow = 5) +
  geom_line() +
  scale_y_continuous(name = "Log-Likelihood", label = comma) +
  scale_x_continuous(name = "Iteration", label = comma) +
  ggtitle("Orange County Data Traceplots")

simulated_log_lik_plot <-
  simulated_log_lik %>%
  rename(Chain = .chain) %>%
  ggplot(aes(.iteration, log_lik)) +
  facet_wrap(. ~ Chain,
             labeller = label_both,
             nrow = 5) +
  geom_line() +
  scale_y_continuous(name = "Log-Likelihood", label = comma) +
  scale_x_continuous(name = "Iteration", label = comma) +
  ggtitle("Simulated Data Traceplots")

main_log_post_plot <-
  main_log_lik %>%
  rename(Chain = .chain) %>%
  ggplot(aes(.iteration, log_post)) +
  facet_wrap(. ~ Chain,
             labeller = label_both,
             nrow = 5) +
  geom_line() +
  scale_y_continuous(name = "Log-Posterior", label = comma) +
  scale_x_continuous(name = "Iteration", label = comma) +
  ggtitle("Orange County Data Traceplots")

simulated_log_post_plot <-
  simulated_log_lik %>%
  rename(Chain = .chain) %>%
  ggplot(aes(.iteration, log_post)) +
  facet_wrap(. ~ Chain,
             labeller = label_both,
             nrow = 5) +
  geom_line() +
  scale_y_continuous(name = "Log-Posterior", label = comma) +
  scale_x_continuous(name = "Iteration", label = comma) +
  ggtitle("Simulated Data Traceplots")

main_results_mcmc_list <- extract_mcmc_list(main_results)
simulated_results_mcmc_list <- extract_mcmc_list(simulated_results)

variable_latex_converter <-
  c(R0_init_est = "$\\tilde{R}_{01}$",
    dur_latent_est = "$\\log(\\gamma)$",
    dur_infec_est = "$-\\log(\\nu)$",
    rho_death_est = "$\\log\\left(\\frac{\\rho}{1 - \\rho}\\right)$",
    phi_death_est = "$-log(\\phi) / 2$",
    alpha_init_est = "$\\tilde{\\alpha}_1$",
    kappa_est = "$-\\log(\\kappa) / 2$",
    ifr_init_est = "$\\tilde{\\eta}_1$",
    sigma_R0_est = "$\\log(\\sigma_{R_0})$",
    sigma_ifr_est = "$\\log(\\sigma_{IFR})$",
    sigma_alpha_est = "$\\log(\\sigma_{\\alpha})$",
    "$\\frac{\\exp(S_0)}{\\exp(S_0) + 1}$" = "$\\frac{\\exp(S_0)}{\\exp(S_0) + 1}$",
    "$\\frac{\\exp(\\tilde{I}_0)}{\\exp(\\tilde{I}_0) + 1}$" = "$\\frac{\\exp(\\tilde{I}_0)}{\\exp(\\tilde{I}_0) + 1}$",
    "$\\frac{\\exp(\\tilde{R}_0)}{\\exp(\\tilde{R}_0) + 1}$" = "$\\frac{\\exp(\\tilde{R}_0)}{\\exp(\\tilde{R}_0) + 1}$",
    "$\\frac{\\exp(\\tilde{D}_0)}{\\exp(\\tilde{D}_0) + 1}$" = "$\\frac{\\exp(\\tilde{D}_0)}{\\exp(\\tilde{D}_0) + 1}$")

main_results_mcmc_summary <-
  MCMCsummary(main_results_mcmc_list) %>%
  rownames_to_column("Variable") %>%
  as_tibble() %>%
  mutate(Variable = variable_latex_converter[Variable],
         n.eff = as.integer(n.eff)) %>%
  select(Variable, "$\\hat{R}$" = Rhat, ESS = n.eff)

simulated_results_mcmc_summary <-
  MCMCsummary(simulated_results_mcmc_list) %>%
  rownames_to_column("Variable") %>%
  as_tibble() %>%
  mutate(Variable = variable_latex_converter[Variable],
         n.eff = as.integer(n.eff)) %>%
  select(Variable, "$\\hat{R}$" = Rhat, ESS = n.eff)


main_results_mcmc_summary %>%
  xtable(caption = "Convergence diagnostics for the main model fit to the Orange County data set.",
         label = "table:main_results_mcmc_summary") %>%
  print(file = "~/Documents/oc_covid19_stemr_manuscript/tables/main_results_mcmc_summary.tex",
        sanitize.text.function = function(x){x}, floating = T, comment = F, include.rownames = F,
  )

simulated_results_mcmc_summary %>%
  xtable(caption = "Convergence diagnostics for the main model fit to the simulated data set.",
         label = "table:simulated_results_mcmc_summary") %>%
  print(file = "~/Documents/oc_covid19_stemr_manuscript/tables/simulated_results_mcmc_summary.tex",
        sanitize.text.function = function(x){x}, floating = T, comment = F, include.rownames = F,
        )

save_plot(filename = "~/Documents/oc_covid19_stemr_manuscript/figures/main_log_lik_plot.pdf",
          plot = main_log_lik_plot,
          ncol = 2,
          nrow = 5,
          base_height =  11 / 5)

save_plot(filename = "~/Documents/oc_covid19_stemr_manuscript/figures/simulated_log_lik_plot.pdf",
          plot = simulated_log_lik_plot,
          ncol = 2,
          nrow = 5,
          base_height =  11 / 5)

save_plot(filename = "~/Documents/oc_covid19_stemr_manuscript/figures/main_log_post_plot.pdf",
          plot = main_log_post_plot,
          ncol = 2,
          nrow = 5,
          base_height =  11 / 5)

save_plot(filename = "~/Documents/oc_covid19_stemr_manuscript/figures/simulated_log_post_plot.pdf",
          plot = simulated_log_post_plot,
          ncol = 2,
          nrow = 5,
          base_height =  11 / 5)
