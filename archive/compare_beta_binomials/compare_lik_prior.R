library(tidyverse)
library(tidybayes)
library(stemr)
library(rstan)

fixed_init_multi_chain_stem_fit <- read_rds("compare_beta_binomials_fixed_tests/fixed_init_multi_chain_stem_fit.rds")
dat <- read_rds("~/Documents/uci_covid_modeling/code/SEIeIpRD/oc/2020-08-12_2020-09-16/model_objects.rds")$lumped_ochca_covid
test_override <- as.integer(mean(dat$new_tests))
dat$new_tests <- test_override

model_objects <- read_rds("~/Documents/uci_covid_modeling/code/SEIeIpRD/oc/2020-08-12_2020-09-16/model_objects.rds")
model_objects$tests[seq_along(dat$new_tests)] <- test_override
model_objects$x_i[1,13:24] <- test_override
model_objects$x_i_pp[1, 29:40] <- model_objects
model_objects$lumped_ochca_covid$new_tests <- test_override
model_objects$lumped_ochca_covid_forecast$new_tests[1:12] <- test_override


dat$new_tests <- as.integer(mean(dat$new_tests))
# fixed_init_multi_chain_stem_fit <- read_rds("fixed_init_sim/fixed_init_multi_chain_stem_fit.rds")
fixed_init_multi_chain_stem_fit <- read_rds("fixed_init_sim/fixed_init_multi_chain_stem_fit.rds")

init_states <- c(S = 3012564L, E = 14866L, Ie = 20907L, Ip = 20590L, R = 0L, D = 0L)
model_objects$popsize <- sum(init_states)

model_objects$frac_carrs <- 1 - init_states[["S"]] / model_objects$popsize
model_objects$frac_carrs_infec <- (init_states[["Ie"]] + init_states[["Ip"]]) / (init_states[["E"]] + init_states[["Ie"]] + init_states[["Ip"]])
model_objects$frac_infec_early <- init_states[["Ie"]] / (init_states[["Ie"]] + init_states[["Ip"]])

i <- 12

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

stan_result <- stan(file = "compare_beta_binomials_fixed_tests/comp_bb_ft_SEIeIpRD.stan",
                data = model_objects,
                iter = 1,
                algorithm = "Fixed_param",
                init = list(fixed_pars),
                seed = 0,
                chains = 1)

unlist(rstan::extract(stan_result, "data_log_lik"))
fixed_init_multi_chain_stem_fit$stem_fit_list[[1]]$results$posterior$data_log_lik[i]

unlist(extract(stan_result, "params_log_prior"))
fixed_init_multi_chain_stem_fit$stem_fit_list[[1]]$results$posterior$params_log_prior[i]
