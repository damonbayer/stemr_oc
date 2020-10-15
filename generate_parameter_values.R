n <- 8000
tibble(log_R0 = rnorm(n, -0.25, 0.7),
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
  as.matrix() %>%
  apply(1, from_estimation_scale) %>%
  t() %>%
  apply(2, median) %>%
  round(4) %>%
  dput()


as.numeric(stemr_prior[1,])

t(apply(as.matrix(stemr_prior), 1, from_estimation_scale))

from_estimation_scale()
