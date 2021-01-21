library(stemr)
library(tidyverse)
library(tidybayes)
library(scales)
library(coda)
library(cowplot)
library(extraDistr)
source('stemr_functions.R')
theme_set(cowplot::theme_minimal_grid())
multi_chain_stem_fit <- read_rds("res_2021_01_13_15_09_32.rds")
multi_chain_stem_fit$n_iterations <- 3000

popsize <- multi_chain_stem_fit$stem_fit_list[[1]]$dynamics$popsize
log_popsize <- log(popsize)

epi_curves_p <- extract_epi_curves(multi_chain_stem_fit = multi_chain_stem_fit, curve_type = "p", tidy = T)

# Separate By Chain
epi_curves_plot <- epi_curves_p %>%
  # filter(time %in% c(0,1)) %>%
  # group_by(.chain, time, name) %>%
  group_by(time, name) %>%
  select(value) %>%
  median_qi(.width = c(0.5, 0.8, 0.95)) %>%
  ggplot(aes(time, value, ymin = .lower, ymax = .upper, group = name)) +
  geom_lineribbon() +
  # facet_grid(name ~ .chain, scales = "free_y") +
  facet_wrap(. ~ name, scales = "free_y") +
  scale_fill_brewer() +
  scale_y_continuous(labels = comma) +
  cowplot::theme_minimal_hgrid() +
  theme(legend.position = "none") +
  ylab("People")

epi_curves_plot


# Beta_t ------------------------------------------------------------------
tmp_beta_t <- imap_dfr(multi_chain_stem_fit$stem_fit_list, ~{
  beta_t <- .x$results$posterior$tparam_samples[,1,,drop = F]
  dur_infec <- exp(-.x$results$posterior$parameter_samples_est[,"dur_infec_est"])

  exp(log(beta_t) + log_popsize + rep(log(dur_infec), each = dim(beta_t)[1])) %>%
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
  theme(legend.position = "none")

Rt_plot

# Posterior Trace Plots ---------------------------------------------------

# Posterior trace plots
posterior_trace_plot <-
  tidy_draws(extract_stem_parameter_posterior(multi_chain_stem_fit = multi_chain_stem_fit, transform = "estimation")) %>%
  pivot_longer(-starts_with(".")) %>%
  mutate(.chain = as_factor(.chain)) %>%
  ggplot(aes(.iteration, value, color = .chain, group = .chain)) +
  geom_line() +
  facet_grid(name ~ .chain, scales = "free_y") +
  theme_cowplot() +
  theme(legend.position = "none")

posterior_trace_plot

posterior_lik_trace_plot <-
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
  facet_grid(name ~ .chain, scales = "free_y") +
  theme_cowplot() +
  theme(legend.position = "none")

posterior_lik_trace_plot

# Posterior Predictive ----------------------------------------------------
stemr_params <- map(multi_chain_stem_fit$stem_fit_list,
                    function(stem_fit) {
                      cbind(
                        stem_fit$results$posterior$parameter_samples_nat,
                        stem_fit$results$posterior$initdist_samples) %>%
                        split_along_dim(1)
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


# Manual Posterior Predictive ---------------------------------------------
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
         alpha = kappa * (exp(alpha0) * (E2I / popsize) ^ alpha1) / (exp(alpha0) * (E2I / popsize) ^ alpha1 + (1 - E2I / popsize) ^ alpha1),
         beta = kappa * ((1 - E2I / popsize) ^ alpha1) / (exp(alpha0) * (E2I / popsize) ^ alpha1 + (1 - E2I / popsize) ^ alpha1),
         size = phi_death,
         mu = rho_death * I2D) %>%
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



ggplot() +
  # geom_lineribbon(data = pp_data2 %>%
  #                   filter(name == "deaths") %>%
  #                   select(time, value) %>%
  #                   group_by(time) %>%
  #                   median_qi(.width = c(0.5, 0.8, 0.95)) %>%
  #                   left_join(select(multi_chain_stem_fit$data, time, date = end_date)),
  #                 mapping = aes(date, value, ymin = .lower, ymax = .upper)) +
  geom_point(data = multi_chain_stem_fit$data,
             mapping = aes(end_date, deaths)) +
  scale_fill_brewer() +
  scale_y_continuous(name = "Deaths", labels = comma) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %e") +
  theme(legend.position = "none")

ggplot() +
  geom_point(data = multi_chain_stem_fit$data,
             mapping = aes(x = end_date, y = deaths))

ggplot(multi_chain_stem_fit$data, aes(end_date, deaths)) + geom_point()
