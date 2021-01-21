library(stemr)
library(tidyverse)
library(tidybayes)
library(scales)
library(coda)
library(cowplot)
library(extraDistr)
source('~/Documents/stemr_oc/stemr_functions.R')
multi_chain_stem_fit <- list()
multi_chain_stem_fit$n_iterations <- 3000
# multi_chain_stem_fit$stem_fit_list <- readRDS("~/Documents/stemr_oc/res_2020_11_02_12_56_25.rds")
# multi_chain_stem_fit$stem_fit_list <- readRDS("~/Documents/stemr_oc/res_2020_11_03_10_44_51.rds") # the incorrect tail
# multi_chain_stem_fit$stem_fit_list <- readRDS("~/Documents/stemr_oc/res_2020_11_03_12_14_06.rds") # just the tail
# multi_chain_stem_fit$stem_fit_list <- read_rds("res_2020_11_05_16_18_06.rds") # correct data, double initial counts
# multi_chain_stem_fit$stem_fit_list <- read_rds("res_2020_11_05_17_31_55.rds") # correct data, triple initial counts
# multi_chain_stem_fit$stem_fit_list <- read_rds("res_2020_11_06_13_53_57.rds") # correct data, quadruple initial counts
# multi_chain_stem_fit$stem_fit_list <- read_rds("res_2020_11_10_13_53_47.rds") # fixed sigma 0.05 triple initial counts
# multi_chain_stem_fit$stem_fit_list <- read_rds("res_2020_11_10_14_44_57.rds") # fixed sigma 0.05 quadruple initial counts
# multi_chain_stem_fit$stem_fit_list <- read_rds("res_2020_11_10_15_56_53.rds") # fixed sigma 0.1 quadruple initial counts
# multi_chain_stem_fit$stem_fit_list <- read_rds("res_2020_11_10_18_04_44.rds") # fixed sigma 0.2 quadruple initial counts
multi_chain_stem_fit$stem_fit_list <- read_rds("res_original_prior_2020_11_12_17_17_19.rds")  # adjusted start date fixed sigma 0.2 original progress prior
# multi_chain_stem_fit$stem_fit_list <- read_rds("res_modified_prior_2020_11_12_17_15_49.rds") # adjusted start date fixed sigma 0.2 longer progress prior
multi_chain_stem_fit$stem_fit_list <- read_rds("res_2020_11_16_17_28_33.rds") # adjusted start date fixed sigma 0.2 longer progress prior, time varying ifr
n_chains <- 4
thinning_interval <- 100

popsize <- 3068927L
log_popsize <- log(popsize)


# Epi Curves --------------------------------------------------------------
epi_curves_p <- extract_epi_curves(multi_chain_stem_fit = multi_chain_stem_fit, curve_type = "p", tidy = T)

# Separate By Chain
epi_curves_p %>%
  filter(time %in% c(0,1)) %>%
  group_by(.chain, time, name) %>%
  select(value) %>%
  median_qi(.width = c(0.5, 0.8, 0.95)) %>%
  ggplot(aes(time, value, ymin = .lower, ymax = .upper, group = name)) +
  geom_lineribbon() +
  facet_grid(name ~ .chain, scales = "free_y") +
  scale_fill_brewer() +
  scale_y_continuous(labels = comma) +
  cowplot::theme_minimal_hgrid() +
  theme(legend.position = "none")

epi_curves_p %>%
  filter(time %in% c(0,1)) %>%
  group_by(time, name) %>%
  select(value) %>%
  median_qi(.width = c(0.5, 0.8, 0.95)) %>%
  ggplot(aes(time, value, ymin = .lower, ymax = .upper, group = name)) +
  geom_lineribbon() +
  facet_grid(. ~ name, scales = "free_y") +
  scale_fill_brewer() +
  scale_y_continuous(labels = comma) +
  cowplot::theme_minimal_hgrid() +
  theme(legend.position = "none")


# Beta_t ------------------------------------------------------------------
tmp_beta_t <- imap_dfr(multi_chain_stem_fit$stem_fit_list, ~{
  beta_t <- .x$results$posterior$tparam_samples
  nu_early <- .x$results$posterior$parameter_samples_nat[,"nu_early"]
  mu_death <- .x$results$posterior$parameter_samples_nat[,"mu_death"]
  mu_rec <- .x$results$posterior$parameter_samples_nat[,"mu_rec"]

  exp(log(beta_t) + log_popsize + rep(log(1 / nu_early + 0.8 / (mu_rec + mu_death)), each = dim(beta_t)[1])) %>%
    gather_array(value = value, time, var, .iteration) %>%
    as_tibble() %>%
    select(-var) %>%
    mutate(time = time - 1) %>%
    mutate(.chain = .y)
}) %>%
    mutate(.draw = tidybayes:::draw_from_chain_and_iteration_(.chain, .iteration)) %>%
    left_join(epi_curves_p %>%
                pivot_wider(names_from = name, values_from = value) %>%
                mutate(prop_s = S / popsize) %>%
                select(time, .iteration, .chain, .draw, prop_s)) %>%
    mutate(Rt = value * prop_s)

epi_estim_output_3_15_10_15 <- readRDS("/Users/damon/Downloads/epi_estim_output_3_15_10_15.rds")
library(EpiEstim)

Rt_plot <- ggplot() +
  geom_lineribbon(data = tmp_beta_t %>%
                    select(.chain, time, Rt) %>%
                    group_by(.chain, time) %>%
                    # select(time, Rt) %>%
                    # group_by(time) %>%
                    median_qi(.width = c(0.5, 0.8, 0.95)),
                  mapping = aes(time, Rt, ymin = .lower, ymax = .upper)) +
  facet_wrap(. ~ .chain) +
  scale_fill_brewer() +
  cowplot::theme_minimal_hgrid() +
  ggtitle("Stemr vs epiestim?")

Rt_plot

cowplot::save_plot(filename = "~/Desktop/stemr_vs_epiestim.pdf", ncol = 1, nrow = 1, plot = Rt_plot)


# Posterior Trace Plots ---------------------------------------------------

# Posterior trace plots
tidy_draws(extract_stem_parameter_posterior(multi_chain_stem_fit = multi_chain_stem_fit, transform = "estimation")) %>%
  pivot_longer(-starts_with(".")) %>%
  mutate(.chain = as_factor(.chain)) %>%
  ggplot(aes(.iteration, value, color = .chain, group = .chain)) +
  geom_line() +
  facet_wrap(. ~ name, scales = "free_y") +
  theme_cowplot() +
  theme(legend.position = "none")



multi_chain_stem_fit$stem_fit_list %>%
  map("results") %>%
  map("posterior") %>%
  map("data_log_lik") %>%
  enframe(name = ".chain", value = "data_log_lik") %>%
  left_join(multi_chain_stem_fit$stem_fit_list %>%
              map("results") %>%
              map("posterior") %>%
              map("data_log_lik") %>%
              enframe(name = ".chain", value = "params_log_prior")) %>%
  unnest(c("data_log_lik", "params_log_prior")) %>%
  mutate(lp = data_log_lik + params_log_prior) %>%
  group_by(.chain) %>%
  mutate(.iteration = row_number()) %>%
  mutate(.draw = tidybayes:::draw_from_chain_and_iteration_(.chain, .iteration)) %>%
  select(starts_with("."), everything()) %>%
  pivot_longer(-starts_with(".")) %>%
  mutate(.chain = as_factor(.chain)) %>%
  ggplot(aes(.iteration, value, color = .chain, group = .chain)) +
  geom_line() +
  facet_wrap(. ~ name, scales = "free_y") +
  theme_cowplot() +
  theme(legend.position = "none")


# Posterior Predictive ----------------------------------------------------
stem_object$dynamics$parameters
stemr_params <- map(multi_chain_stem_fit$stem_fit_list,
                            function(stem_fit) {
                              cbind(
                                stem_fit$results$posterior$parameter_samples_nat,
                                stem_fit$results$posterior$initdist_samples) %>%
                                split(., row(.)) %>%
                                # map(~`names<-`(.,c("log_R0_init", "gamma", "nu_early", "mu_rec", "mu_death", "rho_death", "phi_death", "alpha0", "alpha1", "kappa", "sigma", "S_0", "E_0", "Ie_0", "Ip_0", "R_0", "D_0")))
                              map(~`names<-`(.,c("log_R0_init", "gamma", "nu_early", "mu_rec", "mu_death", "rho_death", "phi_death", "alpha0", "alpha1", "kappa", "S_0", "E_0", "Ie_0", "Ip_0", "R_0", "D_0")))
                            }) %>%
  unlist(recursive = F)

# post processing is a bit confusing since things aren't named


# this should be a list of lists of matrices
# first index is simulation
# second index is tparam
# This sucks
stemr_tparams <- map(multi_chain_stem_fit$stem_fit_list,
                     function(stem_fit) {
                       tmp <- stem_fit$results$posterior$tparam_draws %>%
                         map(~split_along_dim(., n = 2))

                       i <- 1:length(tmp)
                       j <- 1:length(tmp[[1]])
                       lapply(j, function(j) lapply(i, function(i) tmp[[i]][[j]]))
                       }) %>%
  unlist(recursive = F)


# tmp_pp <- simulate_stem(stem_object = stem_object,
#                      nsim = 12000,
#                      simulation_parameters = stemr_params,
#                      tparam_draws = stemr_tparams,
#                      method = "ode", paths = F)
#
# tmp_pp2 <- stemr::simulate_stem(stem_object = stem_object,
#                                nsim = 12000,
#                                simulation_parameters = stemr_params,
#                                tparam_draws = stemr_tparams,
#                                method = "ode", paths = F)
#
# tmp_pp$natural_paths[[1]] == tmp_pp2$natural_paths[[1]]
# tmp_pp$paths[[1]] == tmp_pp2$paths[[1]]


pp_data <- map(1:12000, ~{
  stem_object$dynamics$parameters <- head(stemr_params[[.]], -6)
  stem_object$dynamics$initializer[[1]]$init_states <- tail(stemr_params[[.]], 6)

  simulate_stem(stem_object = stem_object,
                nsim = 1,
                tparam_draws = list(stemr_tparams[[.]]),
                method = "ode", paths = F)
}) %>%
  map("datasets") %>%
  map(pluck(1)) %>%
  do.call(rbind, .) %>%
  as_tibble() %>%
  left_join(select(dat, time, tests)) %>%
  mutate(pos = cases / tests) %>%
  select(time, deaths, pos) %>%
  pivot_longer(-time) %>%
  group_by(time, name) %>%
  median_qi(.width = c(0.5, 0.8, 0.95))


stemr_plot <- ggplot() +
  geom_lineribbon(data = pp_data, mapping = aes(time, value, ymin = .lower, ymax = .upper)) +
  geom_point(data = dat %>%
               mutate(pos = cases / tests) %>%
               select(time, deaths, pos) %>%
               pivot_longer(-time), mapping = aes(time, value)) +
  facet_wrap(. ~ name, scales = "free_y") +
  scale_fill_brewer() +
  theme_bw() +
  ggtitle("Posterior Predictive & Real Data (by stemr)") +
  theme(legend.position = "none")


stemr_plot


# Manual Posterior Predictive ---------------------------------------------
library(coda)
pp_data2 <- imap_dfr(multi_chain_stem_fit$stem_fit_list,  ~{
  .x$results$posterior$latent_paths %>%
    gather_array(value, row, name, .iteration) %>%
    as_tibble() %>%
    mutate(name = dimnames(.x$results$posterior$latent_paths)[[2]][name],
           time = .x$results$posterior$latent_paths[,1,1][row]) %>%
    filter(name != "time") %>%
    select(.iteration, time, name, value, -row) %>%
    pivot_wider(names_from = name, values_from = value) %>%
    mutate(.chain = .y)
}) %>%
  mutate(.draw = tidybayes:::draw_from_chain_and_iteration_(.chain, .iteration)) %>%
  select(.chain, .iteration, .draw, everything()) %>%
  left_join(tidy_draws(extract_stem_parameter_posterior(multi_chain_stem_fit = multi_chain_stem_fit, transform = "natural"))) %>%
  filter(time != 0) %>%
  left_join(select(dat, time, tests)) %>%
  mutate(tests = tests,
            alpha = kappa * (exp(alpha0) * (Ie2Ip / popsize) ^ alpha1) / (exp(alpha0) * (Ie2Ip / popsize) ^ alpha1 + (1 - Ie2Ip/popsize) ^ alpha1),
            beta = kappa * ((1 - Ie2Ip/popsize) ^ alpha1) / (exp(alpha0) * (Ie2Ip / popsize) ^ alpha1 + (1 - Ie2Ip/popsize) ^ alpha1),
            size = phi_death,
            mu = rho_death * Ip2D) %>%
  mutate(
    cases = pmap_dbl(
      list(tests = as.list(tests),
           alpha = as.list(alpha),
           beta = as.list(beta)),
      function(tests, alpha, beta) rbbinom(n = 1, size = tests, alpha = alpha, beta = beta)),
    deaths = pmap_dbl(
      list(mu = as.list(mu),
           size = as.list(size)),
      function(mu, size) rnbinom(n = 1, size = size, mu = mu))
    ) %>%
  select(starts_with("."), time, tests, cases, deaths) %>%
  mutate(pos = cases / tests) %>%
  select(.chain, time, deaths, pos) %>%
  pivot_longer(-c(.chain, time))


handmade_plot <- ggplot() +
  geom_lineribbon(data = pp_data2 %>%
                    group_by(.chain, time, name) %>%
                    median_qi(.width = c(0.5, 0.8, 0.95)),
                  mapping = aes(time, value, ymin = .lower, ymax = .upper)) +
  geom_point(data = dat %>%
               mutate(pos = cases / tests) %>%
               select(time, deaths, pos) %>%
               pivot_longer(-time), mapping = aes(time, value)) +
  facet_grid(name ~ .chain, scales = "free_y") +
  scale_fill_brewer() +
  theme_bw() +
  ggtitle("Posterior Predictive & Real Data (by hand)") +
  theme(legend.position = "none")

cowplot::save_plot(plot = handmade_plot + ggtitle("Posterior Predictive & Real Data (by hand)", subtitle = "longer prior"), filename = "~/Desktop/new_prior.pdf", nrow = 2, ncol = 4)
cowplot::save_plot(plot = handmade_plot + ggtitle("Posterior Predictive & Real Data (by hand)", subtitle = "original prior"), filename = "~/Desktop/original_prior.pdf", nrow = 2, ncol = 4)
pp_dat_tail_only <- cowplot::plot_grid(stemr_plot, handmade_plot, nrow = 2, ncol = 1, align = "hv")
cowplot::save_plot(plot = pp_dat_tail_only, filename = "~/Desktop/pp_dat_tail_only.pdf", nrow = 2, ncol = 2)

pp_dat_plot <- cowplot::plot_grid(stemr_plot, handmade_plot, nrow = 2, ncol = 1, align = "hv")
cowplot::save_plot(plot = pp_dat_plot, filename = "~/Desktop/pp_dat.pdf", nrow = 2, ncol = 2)
