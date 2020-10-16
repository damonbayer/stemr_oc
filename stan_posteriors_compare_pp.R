source('~/Documents/uci_covid_modeling/code/SEIeIpRD/plot_functions.R')
# oc_post <- read_rds("fixed_init_sim/oc_post.rds")
oc_post <- read_rds("fixed_init_sim/oc_prior.rds")
model_objects <- read_rds("fixed_init_sim/model_objects.rds")

oc_post_stan_params <- oc_post %>%
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

stemr_pp <-
  simulate_stem(stem_object = stem_object, method = "ode", nsim = 4000,
                simulation_parameters = oc_post_stan_params) %>%
  `$`(dataset) %>%
  map_dfr(as_tibble) %>%
  pivot_longer(-time) %>%
  group_by(time, name) %>%
  median_qi()


stan_pp <- prep_pp(model_objects, oc_post)$posterior %>%
  filter(usage == "train",
         name %in% c("cases", "deaths")) %>%
  mutate(time = as.numeric((t - min(t) + 3) / 7)) %>%
  select(time, name, value) %>%
  group_by(time, name) %>%
  median_qi()

bind_rows(mutate(stemr_pp, model = "stemr"),
          mutate(stan_pp, model = "stan")) %>%
ggplot(aes(time, value, ymin = .lower, ymax = .upper, fill = model)) +
  geom_lineribbon(alpha = 0.5) +
# facet_grid(name ~ model, scales = "free_y")
  facet_wrap(.~name, scales = "free_y") +
  guides(fill=guide_legend("pp simulation from")) +
  cowplot::theme_minimal_grid() +
  theme(legend.position = "bottom") +
  # ggtitle("Predictives based on Stan Posteriors", subtitle = "95% CI") +
  ggtitle("Predictives based on Stan Priors", subtitle = "95% CI") +
  scale_y_continuous(labels = scales::comma)

