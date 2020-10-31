library(tidyverse)
library(tidybayes)
library(stemr)
library(rstan)

model_objects <- read_rds("~/Documents/uci_covid_modeling/code/SEIeIpRD/oc/2020-08-12_2020-09-16/model_objects.rds")
# fixed_init_multi_chain_stem_fit <- read_rds("fixed_init_sim/fixed_init_multi_chain_stem_fit.rds")
fixed_init_multi_chain_stem_fit <- read_rds("fixed_init_sim/fixed_init_multi_chain_stem_fit_plus_one.rds")

init_states <- c(S = 3012564L, E = 14866L, Ie = 20907L, Ip = 20590L, R = 0L, D = 0L)
model_objects$popsize <- sum(init_states)

model_objects$frac_carrs <- 1 - init_states[["S"]] / model_objects$popsize
model_objects$frac_carrs_infec <- (init_states[["Ie"]] + init_states[["Ip"]]) / (init_states[["E"]] + init_states[["Ie"]] + init_states[["Ip"]])
model_objects$frac_infec_early <- init_states[["Ie"]] / (init_states[["Ie"]] + init_states[["Ip"]])

i <- 10

fixed_pars <- fixed_init_multi_chain_stem_fit$stem_fit_list[[1]]$results$posterior$parameter_samples_est[i,] %>%
  enframe() %>%
  pivot_wider(names_from = name, values_from = value) %>%
  transmute(log_R0 = R0_est,
            latent_dur = exp(-dur_latent_est),
            early_dur = exp(-dur_early_est),
            prog_dur = exp(-dur_progress_est),
            IFR = expit(ifr_est),
            rho_death = expit(rho_death_est),
            sqrt_phi_inv_death = exp(phi_death_est),
            alpha_incid_0 = exp(alpha0_est),
            alpha_incid_1 = expit(alpha1_est),
            sqrt_kappa_inv_incid = exp(kappa_est)) %>%
  pivot_longer(everything()) %>%
  deframe() %>%
  as.list()

stan_result <- stan(file = "~/Documents/stemr_oc/fixed_init_sim/SEIeIpRD.stan",
                data = model_objects,
                iter = 1,
                algorithm = "Fixed_param",
                init = list(fixed_pars),
                seed = 0,
                chains = 1)

unlist(extract(stan_result, "data_log_lik"))
fixed_init_multi_chain_stem_fit$stem_fit_list[[1]]$results$posterior$data_log_lik[i]

unlist(extract(stan_result, "params_log_prior"))
fixed_init_multi_chain_stem_fit$stem_fit_list[[1]]$results$posterior$params_log_prior[i]
