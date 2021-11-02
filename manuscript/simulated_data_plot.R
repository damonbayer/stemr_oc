library(tidyverse)
library(tidybayes)
library(fs)
library(cowplot)
library(scales)
library(stemr)
library(coda)
library(extraDistr)
library(lubridate)
library(latex2exp)

# City --------------------------------------------------------------------
model_name_conversion <- c(
  `2021-02-28` = "Orange County",
  double_initial_counts = "Sensitivity 1",
  higher_ifr = "Sensitivity 4",
  huntington_beach = "Huntington Beach",
  irvine = "Irvine",
  longer_infectious_period = "Sensitivity 2",
  lower_R0 = "Sensitivity 3",
  lower_alpha = "Sensitivity 5",
  santa_ana = "Santa Ana",
  simulated_data = "Simulated Data"
)

n_sim <- 2000
original_priors <- tibble(
  "$S_0$" = expit(rnorm(n_sim, 6.02, 0.55)),
  "$\\tilde{I}_0$" = expit(rnorm(n_sim, 0.62, 0.03)),
  "$\\tilde{R}_0$" = expit(rnorm(n_sim, -7.8977, 9.85e-05)),
  "$\\tilde{D}_0$" = expit(rnorm(n_sim, -7.89730, 9.85e-06)),
  "$1 / \\gamma$" = exp(rnorm(n_sim, -0.85, 0.25)),
  "$1 / \\nu" = exp(rnorm(n_sim, 0, 0.225)),
  "$\\rho$" = rbeta(n_sim, 200, 20),
  "$1 / \\sqrt{\\phi}$" = rexp(n_sim, 1),
  "$1 / \\sqrt{\\kappa}$" = rexp(n_sim, 1),
  "$\\exp{\\tilde{R}_{01}}$" = exp(rnorm(n_sim, 0, 0.225)),
  "$\\expit{\\tilde{\\eta}_1}$" = rbeta(n_sim, 30, 10000),
  "$\\exp{\\tilde{\\alpha}_1}$" = rtnorm(n_sim, 4, 0.5, a = 0),
  "$\\sigma_{R_0}$" = rtnorm(n_sim, 0.05, 0.01, a = 0),
  "$\\sigma_{IFR}$" = rtnorm(n_sim, 0.08, 0.01, a = 0),
  "$\\sigma_{\\alpha}$" = rtnorm(n_sim, 0.05, 0.01, a = 0))


stemr_time_constant_params <-
  c("final_models/results/stemr_all_params/simulated_data_stemr_all_params.rds" = path("final_models/results/stemr_all_params/simulated_data_stemr_all_params.rds")) %>%
  map(~read_rds(.) %>%
        filter(time == 0) %>%
        select(starts_with("."), gamma, dur_infec, rho_death, phi_death, kappa,
               sigma_R0, sigma_ifr, sigma_alpha,
               S, E, I, R, D,
               R0_init, alpha_init, ifr_init)) %>%
  `names<-`(., names(.) %>% path_file() %>% path_ext_remove() %>% str_sub(end = -18)) %>%
  `names<-`(., model_name_conversion[names(.)]) %>%
  imap_dfr(~mutate(.x, model = .y)) %>%
  transmute(model = model,
            "$S_0$" = S / (S + E + I + R + D),
            "$\\tilde{I}_0$" = I / (E + I + R + D),
            "$\\tilde{R}_0$" = R / (E + R + D),
            "$\\tilde{D}_0$" = D / (E + D),
            "$1 / \\gamma$" = 1 / gamma,
            "$1 / \\nu" = dur_infec,
            "$\\rho$" = rho_death,
            "$1 / \\sqrt{\\phi}$" = exp(-log(phi_death)/2),
            "$1 / \\sqrt{\\kappa}$" = exp(-log(kappa)/2),
            "$\\exp{\\tilde{R}_{01}}$" = R0_init,
            "$\\expit{\\tilde{\\eta}_1}$" = ifr_init,
            "$\\exp{\\tilde{\\alpha}_1}$" = alpha_init,
            "$\\sigma_{R_0}$" = sigma_R0,
            "$\\sigma_{IFR}$" = sigma_ifr,
            "$\\sigma_{\\alpha}$" = sigma_alpha) %>%
  pivot_longer(-model) %>%
  mutate(Distribution = "Posterior") %>%
  bind_rows(bind_rows(
    mutate(original_priors,
           model = "Orange County"),
    mutate(original_priors,
           model = "Huntington Beach"),
    mutate(original_priors,
           model = "Irvine"),
    mutate(original_priors,
           model = "Santa Ana"),
    mutate(original_priors,
           model = "Sensitivity 1",
           "$S_0$" = expit(rnorm(n_sim, 5.33, 0.12)),
           "$\\tilde{I}_0$" = expit(rnorm(n_sim, 0.62, 0.06)),
           "$\\tilde{R}_0$" = expit(rnorm(n_sim, -8.59, 2.38e-04)),
           "$\\tilde{D}_0$" = expit(rnorm(n_sim, -8.59, 2.38e-04))),
    mutate(original_priors,
           model = "Sensitivity 2",
           "$1 / \\nu" = exp(rnorm(n_sim, -0.34, 0.20))),
    mutate(original_priors,
           model = "Sensitivity 3",
           "$\\exp{\\tilde{R}_{01}}$" = exp(rnorm(n_sim, -0.69, 0.40))),
    mutate(original_priors,
           model = "Sensitivity 4",
           "$\\expit{\\tilde{\\eta}_1}$" = rbeta(n_sim, 60, 9000)),
    mutate(original_priors,
           model = "Sensitivity 5",
           "$\\exp{\\tilde{\\alpha}_1}$" = rtnorm(n_sim, 2, 0.5, a = 0)),
    mutate(original_priors,
           model = "Simulated Data")) %>%
      pivot_longer(-model) %>%
      mutate(Distribution = "Prior")
  ) %>%
  mutate(name = as_factor(name)) %>%
  mutate(name = fct_relabel(name, ~as.character(TeX(.)))) %>%
  filter(model == "Simulated Data")


simulation_parameters <-
  c(
    R0_init = 1.2, gamma = 3.1, dur_infec = 0.68, rho_death = 0.92,
    phi_death = 2800, alpha_init = 4.1, kappa = 1.9e+08, ifr_init = 0.0022,
    sigma_R0 = 0.054, sigma_ifr = 0.086, sigma_alpha = 0.072
  )

init_dist <- c(S = 3170103, I = 3594, R = 1, D = 1, E = 1993)

true_values <-
  enframe(c(simulation_parameters, init_dist)) %>%
  pivot_wider(everything()) %>%
  transmute("$S_0$" = S / (S + E + I + R + D),
            "$\\tilde{I}_0$" = I / (E + I + R + D),
            "$\\tilde{R}_0$" = R / (E + R + D),
            "$\\tilde{D}_0$" = D / (E + D),
            "$1 / \\gamma$" = 1 / gamma,
            "$1 / \\nu" = dur_infec,
            "$\\rho$" = rho_death,
            "$1 / \\sqrt{\\phi}$" = exp(-log(phi_death)/2),
            "$1 / \\sqrt{\\kappa}$" = exp(-log(kappa)/2),
            "$\\exp{\\tilde{R}_{01}}$" = R0_init,
            "$\\expit{\\tilde{\\eta}_1}$" = ifr_init,
            "$\\exp{\\tilde{\\alpha}_1}$" = alpha_init,
            "$\\sigma_{R_0}$" = sigma_R0,
            "$\\sigma_{IFR}$" = sigma_ifr,
            "$\\sigma_{\\alpha}$" = sigma_alpha) %>%
  pivot_longer(everything()) %>%
  mutate(name = as_factor(name)) %>%
  mutate(name = fct_relabel(name, ~as.character(TeX(.))))

time_constant_plot <-
  stemr_time_constant_params %>%
  mutate(model = model %>%
           fct_inorder() %>%
           fct_rev()) %>%
  ggplot(aes(value, Distribution, group = Distribution, fill = Distribution, color = Distribution)) +
  facet_wrap(. ~ name, scales = "free", labeller = label_parsed, ncol = 3) +
  stat_halfeye(normalize = "xy", alpha = 0.5) +
  geom_vline(data = true_values, mapping = aes(xintercept = value), linetype = "dashed") +
  labs(x = NULL, y = NULL) +
  cowplot::theme_minimal_vgrid() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(reverse = TRUE),
         fill = guide_legend(reverse = TRUE)) +
  ggtitle("Prior and Posterior Distributions for Simulated Data Analysis",
          subtitle = "Dashed line indicates true value used in simulation")



save_plot(filename = "~/Documents/oc_covid19_stemr_manuscript/figures/simulated_data_time_constant_plot.pdf",
          plot = time_constant_plot,
          ncol = 3,
          nrow = 5,
          base_asp = 1.5,
          # base_asp = (8.5 * 3) / (11 * 5),
          base_height = 2.5)
