library(tidyverse)
library(tidybayes)
library(stemr)
library(coda)
library(scales)
library(fs)
library(extraDistr)
theme_set(theme_minimal())
source('stemr_functions.R')
source("~/Documents/uci_covid_modeling/code/SEIeIpRD/plot_functions.R")

# stan_oc_post <- read_rds("~/Documents/uci_covid_modeling/code/SEIeIpRD/oc/2020-08-12_2020-09-16/oc_post.rds")
# multi_chain_stem_fit <- read_rds("fixed_init_sim/multi_chain_stem_fit.rds")

# multi_chain_stem_fit <- read_rds("fixed_init_sim/multi_chain_stem_fit_R0D0.rds")
# stan_oc_post <- read_rds("fixed_init_sim/oc_post.rds")
# stan_oc_model_objects <- read_rds("fixed_init_sim/model_objects.rds")
# stan_oc_prior <- read_rds("fixed_init_sim/oc_prior.rds")


# multi_chain_stem_fit <- read_rds("fixed_init_sim/multi_chain_stem_fit_R0D0_bugfix.rds")
# stan_oc_post <- read_rds("fixed_init_sim/oc_post_bugfix.rds")
# stan_oc_model_objects <- read_rds("fixed_init_sim/model_objects.rds")
# stan_oc_prior <- read_rds("fixed_init_sim/oc_prior_bugfix.rds")
# stemr_prior <- extract_stem_parameter_posterior(read_rds("priors_only_sim/priors_only_multi_chain_stem_fit_2020_10_16_10_56_13.rds")) %>%
#   tidy_draws() %>%
#   select(-starts_with(".")) %>%
#   mutate(model = "stemr", type = "prior")

multi_chain_stem_fit <- read_rds("mvnss_sim/mvnss_multi_chain_stem_fit_2020_10_16_12_25_59.rds")
stan_oc_post <- read_rds("fixed_init_sim/oc_post_bugfix.rds")
stan_oc_model_objects <- read_rds("fixed_init_sim/model_objects.rds")
stan_oc_prior <- read_rds("fixed_init_sim/oc_prior_bugfix.rds")
stemr_prior <- extract_stem_parameter_posterior(read_rds("priors_only_sim/priors_only_multi_chain_stem_fit_2020_10_16_10_56_13.rds")) %>%
  tidy_draws() %>%
  select(-starts_with(".")) %>%
  mutate(model = "stemr", type = "prior")

# prior_post_estim_scale_plot ---------------------------------------------
stemr_post <- extract_stem_parameter_posterior(multi_chain_stem_fit) %>%
  tidy_draws() %>%
  select(-starts_with(".")) %>%
  mutate(model = "stemr", type = "post")

stan_post <- stan_oc_post %>%
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
  mutate(model = "stan", type = "post")


stan_prior <- stan_oc_prior %>%
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
            dur_latent_est = log(latent_dur),
            dur_early_est = log(early_dur),
            dur_progress_est = log(prog_dur),
            ifr_est = logit(IFR),
            rho_death_est = logit(rho_death),
            phi_death_est = log(sqrt_phi_inv_death),
            alpha0_est = log(alpha_incid_0),
            alpha1_est = logit(alpha_incid_1),
            kappa_est = log(sqrt_kappa_inv_incid)) %>%
  mutate(model = "stan", type = "prior")

priors_only_plot <- bind_rows(stemr_prior, stemr_post, stan_prior, stan_post) %>%
  filter(type == "prior") %>%
  unite(type, model, col = "source") %>%
  pivot_longer(-source) %>%
  group_by(source, name) %>%
  median_qi(.width = c(0.5, 0.95)) %>%
  ggplot(aes(x = value, xmin = .lower, xmax = .upper, color = source)) +
  geom_pointinterval(position = position_dodge(width = .3)) +
  facet_wrap(. ~ name, scales = "free_x") +
  cowplot::theme_minimal_vgrid() +
  xlab(NULL) +
  scale_y_continuous(labels = NULL) +
  ggtitle("Fixed Inits") +
  theme(legend.position = c(.6, .2)) +
  guides(colour = guide_legend(reverse = T))

prior_post_estim_scale_plot <-
  bind_rows(stemr_prior, stemr_post, stan_prior, stan_post) %>%
  unite(type, model, col = "source") %>%
  pivot_longer(-source) %>%
  group_by(source, name) %>%
  median_qi(.width = c(0.5, 0.95)) %>%
  ggplot(aes(x = value, xmin = .lower, xmax = .upper, color = source)) +
  geom_pointinterval(position = position_dodge(width = .3)) +
  facet_wrap(. ~ name, scales = "free_x") +
  cowplot::theme_minimal_vgrid() +
  xlab(NULL) +
  scale_y_continuous(labels = NULL) +
  ggtitle("Fixed Inits") +
  theme(legend.position = c(.6, .2)) +
  guides(colour = guide_legend(reverse = T))


# Stemr trace plot --------------------------------------------------------
stemr_trace_plot <- extract_stem_parameter_posterior(multi_chain_stem_fit) %>%
  tidy_draws() %>%
  pivot_longer(-c(.chain, .iteration, .draw)) %>%
  ggplot(aes(.iteration, value, color = as_factor(.chain))) +
  geom_line() +
  facet_grid(name ~ .chain, scales = "free_y") +
  theme(legend.position = "none")



# Prevalence Curves -------------------------------------------------------
stemr_epi_curves <- extract_epi_curves(multi_chain_stem_fit, curve_type = "p", tidy = T) %>%
  # filter(name %in% c("S", "E", "Ie", "Ip")) %>%
  select(time, name, value) %>%
  mutate(model = "stemr")


stan_epi_curves <-
prep_epi_curves(stan_oc_post, stan_oc_model_objects) %>%
  pivot_wider( values_from = epi_curves, names_from = g_text) %>%
  select(.draw, t, S, E, IE, IP, D) %>%
  mutate(R = 1 - (S + E + IE + IP + D)) %>%
  pivot_longer(-c(.draw, t)) %>%
  rename(time = t) %>%
  ungroup() %>%
  mutate(time = as.numeric(time - min(time) + 3) / 7,
         name = str_to_title(name)) %>%
  filter(time %in% unique(stemr_epi_curves$time)) %>%
  # filter(name %in% unique(stemr_epi_curves$name)) %>%
  select(time, name, value) %>%
  bind_rows(spread_draws(stan_oc_post, `init_fracs[1]`, `init_fracs[2]`, `init_fracs[3]`) %>%
              mutate(S = 1 - `init_fracs[1]`,
                     E = exp(log(`init_fracs[1]`) + log(1 - `init_fracs[2]`)),
                     Ie = exp(log(`init_fracs[1]`) + log(`init_fracs[2]`) + log(`init_fracs[3]`)),
                     Ip = exp(log(`init_fracs[1]`) + log(`init_fracs[2]`) + log(1-`init_fracs[3]`)),
                     R = 0,
                     D = 0) %>%
              select(-starts_with("init"), -starts_with(".")) %>%
              pivot_longer(everything()) %>%
              mutate(time = 0), .) %>%
  mutate(value = value * stan_oc_model_objects$popsize) %>%
  mutate(name = fct_inorder(name),
         model = "stan")

prevalence_curves_plot <-
  bind_rows(stemr_epi_curves, stan_epi_curves) %>%
  group_by(name, model, time) %>%
  median_qi() %>%
  ggplot(aes(time, value, ymin = .lower, ymax = .upper, fill = model, color = model, group = model)) +
  facet_wrap(.~name, scale = "free_y") +
  # facet_grid(name ~ model, scale = "free_y") +
  geom_lineribbon(alpha = 0.5) +
  ylab("People") +
  scale_y_continuous(labels = comma) +
  ggtitle("Prevalence Curves for Stan and Stemr") +
  theme(legend.position = "bottom")


# Save --------------------------------------------------------------------
folder_name <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
dir_create(path = path("images", folder_name))

dpi <- 144

ggsave(filename = "prior_post_estim_scale_plot.png",
       plot = prior_post_estim_scale_plot, device = "png",
       path =  path("images", folder_name), width = 1000 / dpi, height = 800 / dpi, dpi = dpi)


 ggsave(filename = "prevalence_curves_plot.png",
       plot = prevalence_curves_plot, device = "png",
       path =  path("images", folder_name), width = 800 / dpi, height = 800 / dpi, dpi = dpi)
