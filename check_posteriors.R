library(tidyverse)
library(tidybayes)
library(stemr)
library(coda)
theme_set(theme_minimal())
source('stemr_functions.R')
source("~/Documents/uci_covid_modeling/code/SEIeIpRD/plot_functions.R")

# stan_oc_post <- read_rds("/Users/damon/Documents/uci_covid_modeling/code/SEIeIpRD/oc/2020-08-12_2020-09-16/oc_post.rds")

# multi_chain_stem_fit <- read_rds("fixed_init_sim/multi_chain_stem_fit.rds")

multi_chain_stem_fit <- read_rds("fixed_init_sim/multi_chain_stem_fit_R0D0.rds")
# stan_oc_post <- read_rds("fixed_init_sim/oc_post.rds")
# stan_oc_model_objects <- read_rds("fixed_init_sim/model_objects.rds")

stan_oc_post <- read_rds("fixed_init_sim/oc_post_more_precise.rds")
stan_oc_model_objects <- read_rds("fixed_init_sim/model_objects.rds")


# Human Scale -------------------------------------------------------------
human_scale_mcmc_list <- extract_stem_parameter_posterior(multi_chain_stem_fit, transform = to_human_scale)

compare_posteriors <- stan_oc_post %>%
  spread_draws(alpha_incid_0,
               alpha_incid_1,
               early_dur,
               latent_dur,
               prog_dur,
               IFR,
               sqrt_kappa_inv_incid,
               sqrt_phi_inv_death,
               log_R0,
               rho_death) %>%
  select(alpha0 = alpha_incid_0,
         alpha1 = alpha_incid_1,
         dur_early = early_dur,
         dur_latent = latent_dur,
         dur_progress = prog_dur,
         ifr = IFR,
         kappa = sqrt_kappa_inv_incid,
         phi_death = sqrt_phi_inv_death,
         R0 = log_R0,
         rho_death = rho_death) %>%
  mutate(R0 = exp(R0)) %>%
  mutate(model = "stan") %>%
  bind_rows(tidy_draws(human_scale_mcmc_list) %>%
              select(-starts_with(".")) %>%
              mutate(model = "stemr"))


# Estimation Scale --------------------------------------------------------
estim_scale_mcmc_list <- extract_stem_parameter_posterior(multi_chain_stem_fit)

compare_posteriors <- stan_oc_post %>%
  spread_draws(alpha_incid_0,
               alpha_incid_1,
               early_dur,
               latent_dur,
               prog_dur,
               IFR,
               sqrt_kappa_inv_incid,
               sqrt_phi_inv_death,
               log_R0,
               rho_death) %>%
  select(alpha0 = alpha_incid_0,
         alpha1 = alpha_incid_1,
         dur_early = early_dur,
         dur_latent = latent_dur,
         dur_progress = prog_dur,
         ifr = IFR,
         kappa = sqrt_kappa_inv_incid,
         phi_death = sqrt_phi_inv_death,
         R0 = log_R0,
         rho_death = rho_death) %>%
  mutate(alpha0 = log(alpha0),
         alpha1 = logit(alpha1),
         dur_early = -log(dur_early),
         dur_latent = -log(dur_latent),
         dur_progress = -log(dur_progress),
         ifr = logit(ifr),
         kappa = log(kappa),
         phi_death = log(phi_death),
         R0 = R0,
         rho_death = logit(rho_death)) %>%
  mutate(model = "stan") %>%
  bind_rows(tidy_draws(estim_scale_mcmc_list) %>%
              rename_all(~str_remove(., "_est")) %>%
              select(-starts_with(".")) %>%
              mutate(model = "stemr"))


compare_posteriors %>%
  pivot_longer(-model) %>%
  group_by(model, name) %>%
  median_qi() %>%
  ggplot(aes(x = value, xmin = .lower, xmax = .upper, color = model)) +
  geom_pointinterval(position = position_dodge(width = .3)) +
  facet_wrap(. ~ name, scales = "free_x") +
  cowplot::theme_minimal_vgrid() +
  xlab(NULL) +
  ggtitle("Fixed Inits")


tidy_draws(extract_stem_parameter_posterior(multi_chain_stem_fit)) %>%
  rename_all(~str_remove(., "_est")) %>%
  pivot_longer(-c(.chain, .iteration, .draw)) %>%
  ggplot(aes(.iteration, value, color = as_factor(.chain))) +
  geom_line() +
  facet_grid(name ~ .chain, scales = "free_y") +
  theme(legend.position = "none")

# add prior intervals to these plots

# Maybe just a discrepency in the priors
# 1. Latent trajectories (prevalence is problably most important)
# 2. fix initial conditions
# 3. scatterplot matrix of initial conditions
# Maybe I comparment is too high
