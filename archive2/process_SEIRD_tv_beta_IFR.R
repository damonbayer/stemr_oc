library(stemr)
library(tidyverse)
library(tidybayes)
library(scales)
library(coda)
library(cowplot)
library(extraDistr)
source('~/Documents/stemr_oc/stemr_functions.R')

# oc_data <- read_rds("data/oc_data_no_snf.rds")
oc_data <- read_rds("data/oc_data.rds")

lump_oc_data <-
  function(oc_data,
           time_interval_in_days,
           first_day = "0000-01-01",
           last_day = "9999-12-31") {
    oc_data %>%
      filter(date >= lubridate::ymd(first_day),
             date <= lubridate::ymd(last_day)) %>%
      group_by(lump = as.integer(floor((max(date) - date) / time_interval_in_days))) %>%
      filter(n() == time_interval_in_days) %>%
      dplyr::summarize(start_date = min(date),
                       end_date = max(date),
                       cases = sum(cases),
                       tests = sum(tests),
                       deaths = sum(deaths)) %>%
      dplyr::select(-lump) %>%
      arrange(start_date)
  }

lumped_oc_data <- lump_oc_data(oc_data, time_interval_in_days = 7, first_day = "2020-03-30", last_day = "2020-10-11") %>%
  left_join(select(oc_data, end_date = date, lagged_pos_ifr))

dat <- lumped_oc_data %>%
  mutate(time = as.numeric((end_date - min(start_date) + 1L) / 7)) %>%
  select(time, cases, tests, deaths, lagged_pos_ifr)

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
# multi_chain_stem_fit$stem_fit_list <- read_rds("res_original_prior_2020_11_12_17_17_19.rds")  # adjusted start date fixed sigma 0.2 original progress prior
# multi_chain_stem_fit$stem_fit_list <- read_rds("res_modified_prior_2020_11_12_17_15_49.rds") # adjusted start date fixed sigma 0.2 longer progress prior
# multi_chain_stem_fit$stem_fit_list <- read_rds("res_2020_11_16_17_28_33.rds") # adjusted start date fixed sigma 0.2 longer progress prior, time varying ifr
# multi_chain_stem_fit$stem_fit_list <- read_rds("res_2020_11_18_17_07_30.rds") # adjusted start date fixed sigma 0.2 longer progress prior, time varying ifr, no nursing homes
# multi_chain_stem_fit$stem_fit_list <- read_rds("res_2020_11_19_16_44_46.rds") # adjusted start date fixed sigma 0.2 longer progress prior, time varying ifr, reduced initial infections (10x), nursing homes included
# multi_chain_stem_fit$stem_fit_list <- read_rds("res_2020_11_20_15_08_20.rds") # adjusted start date fixed sigma 0.2 longer progress prior, time varying ifr, reduced initial infections (7x), nursing homes included
multi_chain_stem_fit$stem_fit_list <- read_rds("res_2020_11_23_16_40_41.rds")
n_chains <- 4
thinning_interval <- 100

popsize <- 3068927L
log_popsize <- log(popsize)


# Epi Curves --------------------------------------------------------------
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
  dur_infec <- .x$results$posterior$parameter_samples_nat[,"dur_infec"]

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

# IFR_t -------------------------------------------------------------------
IFRt_plot <- imap_dfr(multi_chain_stem_fit$stem_fit_list, ~{
  .x$results$posterior$tparam_samples[,2,,drop=F] %>%
    gather_array(value = value, time, var, .iteration) %>%
    as_tibble() %>%
    select(-var) %>%
    mutate(time = time - 1) %>%
    mutate(.chain = .y)}) %>%
  group_by(time, .chain) %>%
  median_qi(.width = c(0.5, 0.8, 0.95)) %>%
  ggplot(aes(time, value, ymin = .lower, ymax = .upper)) +
  facet_wrap(. ~ .chain) +
  geom_lineribbon() +
  scale_fill_brewer() +
  cowplot::theme_minimal_hgrid() +
  theme(legend.position = "none") +
  geom_line(data = dat %>% mutate(time = time - 1), mapping = aes(time, lagged_pos_ifr, ymin = NULL, ymax = NULL), color = "red") +
  ylab("IFR")

IFRt_plot

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


handmade_pp_plot <- ggplot() +
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

handmade_pp_plot

cowplot::save_plot(plot = handmade_plot + ggtitle("Posterior Predictive & Real Data (by hand)", subtitle = "longer prior"), filename = "~/Desktop/new_prior.pdf", nrow = 2, ncol = 4)
cowplot::save_plot(plot = handmade_plot + ggtitle("Posterior Predictive & Real Data (by hand)", subtitle = "original prior"), filename = "~/Desktop/original_prior.pdf", nrow = 2, ncol = 4)
pp_dat_tail_only <- cowplot::plot_grid(stemr_plot, handmade_plot, nrow = 2, ncol = 1, align = "hv")
cowplot::save_plot(plot = pp_dat_tail_only, filename = "~/Desktop/pp_dat_tail_only.pdf", nrow = 2, ncol = 2)

pp_dat_plot <- cowplot::plot_grid(stemr_plot, handmade_plot, nrow = 2, ncol = 1, align = "hv")
cowplot::save_plot(plot = pp_dat_plot, filename = "~/Desktop/pp_dat.pdf", nrow = 2, ncol = 2)
