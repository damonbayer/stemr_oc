library(stemr)
library(tidyverse)
source('stemr_functions.R')
source("/Users/damon/Documents/uci_covid_modeling/code/SEIeIpRD/plot_functions.R")
source("/Users/damon/Documents/stemr_oc/stemr_functions.R")
oc_prior <-  read_rds("fixed_init_sim/oc_prior.rds")
model_objects <-  read_rds("fixed_init_sim/model_objects.rds")
phi_death_inds <- 1:model_objects$n_obs_pp * 10
kappa_incid_inds  = 1:model_objects$n_obs_pp * 10 - 1
log_mu_death_inds = 1:model_objects$n_obs_pp * 10 - 2
rho_incid_inds    = 1:model_objects$n_obs_pp * 10 - 3

popsize <- model_objects$popsize

names(oc_prior)[str_detect(names(oc_prior), "rho")]

stan_results <- prep_epi_curves(oc_prior, model_objects) %>%
  filter(g_text %in% c("rho incidence", "kappa incid")) %>%
  ungroup() %>%
  pivot_wider(names_from = g_text, values_from = epi_curves) %>%
  rename(rho_incid = `rho incidence`, kappa_incid = `kappa incid`) %>%
  full_join(spread_draws(oc_prior, alpha_incid_0, alpha_incid_1) %>% select(.draw, starts_with("alpha"))) %>%
  mutate(bb_alpha = expit(alpha_incid_0 + alpha_incid_1 * logit(rho_incid)) * kappa_incid,
         bb_beta = (1 - expit(alpha_incid_0 + alpha_incid_1 * logit(rho_incid))) * kappa_incid)


beta_binomial_lpmf(
  x_i[1:n_obs_times] |
    x_i[(1 + n_obs_times):(2 * n_obs_times)],
  kappa_incid * rho_incid,
  kappa_incid * (1 - rho_incid)) +
  neg_binomial_2_log_lpmf(
    x_i[(1 + 2 * n_obs_times):(3 * n_obs_times)] |
      log_mu_death, phi_death));

neg_binomial_2_log_rng(
  log(rho_death) + log_popsize + log(epi_curves[log_mu_death_inds[j]]),
  epi_curves[phi_death_inds[j]]);

# inv_logit(alpha_incid_0 + alpha_incid_1 * logit(epi_curves[rho_incid_inds[j]])) * epi_curves[kappa_incid_inds[j]]
# (1 - inv_logit(alpha_incid_0 + alpha_incid_1 * logit(epi_curves[rho_incid_inds[j]]))) * epi_curves[kappa_incid_inds[j]])



# -------------------------------------------------------------------------

stan_params <- oc_prior %>%
  spread_draws(alpha_incid_0,
               alpha_incid_1,
               early_dur,
               IFR,
               latent_dur,
               prog_dur,
               prog_dur,
               sqrt_kappa_inv_incid,
               sqrt_phi_inv_death,
               log_R0,
               rho_death) %>%
  transmute(R0_est = log_R0,
            dur_latent_est = -log(latent_dur),
            dur_early_est = -log(early_dur),
            dur_progress_est = -log(prog_dur),
            ifr_est = logit(IFR),
            rho_death_est = logit(rho_death),
            phi_death_est = log(sqrt_phi_inv_death),
            alpha0_est = log(alpha_incid_0),
            alpha1_est = logit(alpha_incid_1),
            kappa_est = log(sqrt_kappa_inv_incid)) %>%
  as.matrix() %>%
  apply(1, from_estimation_scale) %>%
  t() %>%
  split(., row(.)) %>%
  map(~`names<-`(.,c("beta", "gamma", "nu_early", "mu_rec", "mu_death", "rho_death",
                     "phi_death", "alpha0", "alpha1", "kappa")))



stemr_prior <- simulate_stem(stem_object = stem_object, method = "ode", nsim = 4000,
              simulation_parameters = stan_params)

write_rds(stemr_prior, "compare_beta_binomials/stemr_prior.rds")
stemr_prior <- read_rds("compare_beta_binomials/stemr_prior.rds")



stemr_results <- imap_dfr(stemr_prior$paths, ~as_tibble(.x) %>% select(time, Ie2Ip) %>% mutate(.draw = .y)) %>%
  left_join(stan_params %>%
              bind_rows() %>%
              select(alpha0, alpha1, kappa) %>%
              mutate(.draw = row_number())) %>%
  filter(time != 0) %>%
  mutate(time = model_objects$dates[time * 7 / 3]) %>%
  mutate(bb_alpha = kappa * (exp(alpha0) * (Ie2Ip / model_objects$popsize) ^ alpha1) / (exp(alpha0) * (Ie2Ip / model_objects$popsize) ^ alpha1 + (1 - Ie2Ip/model_objects$popsize) ^ alpha1),
         bb_beta = kappa * ((1 - Ie2Ip/model_objects$popsize) ^ alpha1) / (exp(alpha0) * (Ie2Ip / model_objects$popsize) ^ alpha1 + (1 - Ie2Ip/model_objects$popsize) ^ alpha1))


tmp <-
left_join(select(stan_results, .draw, t, starts_with("bb")) %>%
  pivot_longer(cols = starts_with("bb"), names_prefix = "bb_", names_to = "param", values_to = "stan"),
select(stemr_results, .draw, t = time, starts_with("bb")) %>%
  pivot_longer(cols = starts_with("bb"), names_prefix = "bb_", names_to = "param", values_to = "stemr")) %>%
  drop_na() %>%
  mutate(diff = stemr - stan)


tmp %>%
  group_by(t, param) %>%
  summarize(mean(abs(diff)),
            var(abs(diff))) %>%
  View()


group_by(param) %>%
  summarize(`mean diff` = mean(diff),
            `sd diff` = sd(diff))


"kappa * (exp(alpha0) * (Ie2Ip / popsize) ^ alpha1) / (exp(alpha0) * (Ie2Ip / popsize) ^ alpha1 + (1 - Ie2Ip/popsize) ^ alpha1)",
"kappa * ((1 - Ie2Ip/popsize) ^ alpha1) / (exp(alpha0) * (Ie2Ip / popsize) ^ alpha1 + (1 - Ie2Ip/popsize) ^ alpha1)"), # distribution pars, here overdispersion and mean
