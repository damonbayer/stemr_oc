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

multi_chain_stem_fit <- read_rds("fixed_init_sim/multi_chain_stem_fit_R0D0.rds")
stan_oc_post <- read_rds("fixed_init_sim/oc_post.rds")
stan_oc_model_objects <- read_rds("fixed_init_sim/model_objects.rds")
stan_oc_prior <- read_rds("fixed_init_sim/oc_prior.rds")


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

n <- 8000

stemr_prior <- tibble(log_R0 = rnorm(n, -0.25, 0.7),
                      log_dur_latent = rnorm(n, 0, 0.22),
                      log_dur_early = rnorm(n, 0, 0.22),
                      log_dur_progress = rnorm(n, 0, 0.22),
                      ifr = rbeta(n, 1.5, 200),
                      rho_death = rbeta(n, 8, 2),
                      sqrt_phi_death_inv = rexp(n, 1),
                      alpha0 = rtnorm(n, 4, 2, a = 0),
                      alpha1 = rbeta(n, 3, 1),
                      sqrt_kappa_inv = rexp(n, 1)) %>%
  transmute(R0_est = log_R0,
            dur_latent_est = log_dur_latent,
            dur_early_est = log_dur_early,
            dur_progress_est = log_dur_progress,
            ifr_est = logit(ifr),
            rho_death_est = logit(rho_death),
            phi_death_est = log(sqrt_phi_death_inv),
            alpha0_est = alpha0,
            alpha1_est = logit(alpha1),
            kappa_est = log(sqrt_kappa_inv)) %>%
  mutate(model = "stemr", type = "prior")

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
            alpha0_est = alpha_incid_0,
            alpha1_est = logit(alpha_incid_1),
            kappa_est = log(sqrt_kappa_inv_incid)) %>%
  mutate(model = "stan", type = "prior")


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
  filter(name %in% c("S", "E", "Ie", "Ip")) %>%
  select(time, name, value) %>%
  mutate(model = "stemr")

stan_epi_curves <- prep_epi_curves(stan_oc_post, stan_oc_model_objects) %>%
  rename(name = g_text,
         value = epi_curves,
         time = t) %>%
  ungroup() %>%
  mutate(time = as.numeric(time - min(time) + 3) / 7,
         name = str_to_title(name)) %>%
  filter(time %in% unique(stemr_epi_curves$time)) %>%
  filter(name %in% unique(stemr_epi_curves$name)) %>%
  mutate(value = value * stan_oc_model_objects$popsize) %>%
  select(time, name, value) %>%
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
