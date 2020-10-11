n <- 8000
stemr_oc_prior <- tibble(R0 = rnorm(n, -0.25, 0.7), # log(R0)
                       dur_latent = rnorm(n, 0, 0.22), # -log(dur_latent)
                       dur_early = rnorm(n, 0, 0.22), # -log(dur_early)
                       dur_progress = rnorm(n, 0, 0.22), # -log(dur_progress)
                       ifr = rbeta(n, 1.5, 200), # logit(ifr)
                       rho_death = rbeta(n, 8, 2), # logit(rho_death)
                       phi_death = rexp(n, 1), # -0.5 * log(phi_death)
                       alpha0 = rnorm(n, 4, 2), # log(alpha0)
                       alpha1 = rbeta(n, 3, 1), # logit(alpha1)
                       kappa = rexp(n, 1)) # -0.5 * log(kappa)

stan_oc_model_objects <- read_rds("fixed_init_sim/model_objects.rds")

stan_oc_prior <- read_rds("fixed_init_sim/oc_prior.rds") %>%
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
  mutate(R0 = exp(R0),
         kappa = 1 / kappa^2,
         phi_death = 1 / phi_death^2) %>%
  mutate(model = "stan")


stemr_oc_prior %>%
  mutate(R0 = exp(R0),
         dur_early = exp(-dur_early),
         dur_latent = exp(-dur_latent),
         dur_progress = exp(-dur_progress),
         kappa = 1 / kappa^2,
         phi_death = 1 / phi_death^2) %>%
  mutate(model = "stemr") %>%
bind_rows(stan_oc_prior) %>%
pivot_longer(-model) %>%
  group_by(model, name) %>%
  median_qi(.width = c(0.5, 0.8, 0.95)) %>%
  ggplot(aes(x = value, xmin = .lower, xmax = .upper, color = model)) +
  geom_pointinterval(position = position_dodge(width = .3)) +
  facet_wrap(. ~ name, scales = "free_x") +
  cowplot::theme_minimal_vgrid() +
  xlab(NULL) +
  ggtitle("priors")
