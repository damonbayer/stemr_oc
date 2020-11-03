library(tidyverse)
library(tidybayes)
library(scales)
source('~/Documents/stemr_oc/stemr_functions.R')
multi_chain_stem_fit <- list()
multi_chain_stem_fit$n_iterations <- 3000
multi_chain_stem_fit$stem_fit_list <- readRDS("~/Documents/stemr_oc/res_2020_11_02_12_56_25.rds")
n_chains <- 4
thinning_interval <- 100

popsize <- 3068927L
log_popsize <- log(popsize)


# Epi Curves --------------------------------------------------------------
epi_curves_p <- extract_epi_curves(multi_chain_stem_fit = multi_chain_stem_fit, curve_type = "p", tidy = T) %>%
  group_by(time, name) %>%
  select(value) %>%
  median_qi(.width = c(0.5, 0.8, 0.95))

epi_curves_p %>%
  ggplot(aes(time, value, ymin = .lower, ymax = .upper, group = name)) +
  geom_lineribbon() +
  facet_wrap(. ~ name, scales = "free_y") +
  scale_fill_brewer() +
  scale_y_continuous(labels = comma) +
  cowplot::theme_minimal_hgrid()
# kinda sus


# Beta_t ------------------------------------------------------------------
beta_t <- multi_chain_stem_fit$stem_fit_list[[1]]$results$posterior$tparam_samples
nu_early <- multi_chain_stem_fit$stem_fit_list[[1]]$results$posterior$parameter_samples_nat[,"nu_early"]
mu_death <- multi_chain_stem_fit$stem_fit_list[[1]]$results$posterior$parameter_samples_nat[,"mu_death"]
mu_rec <- multi_chain_stem_fit$stem_fit_list[[1]]$results$posterior$parameter_samples_nat[,"mu_rec"]

tmp <- exp(log(beta_t) + log_popsize + rep(log(1 / nu_early + 0.8 / (mu_rec + mu_death)), each = 31)) %>%
  gather_array(value = value, time, var, .draw) %>%
  as_tibble() %>%
  select(-var) %>%
  mutate(time = time - 1) %>%
  group_by(time) %>%
  median_qi(.width = c(0.5, 0.8, 0.95))

tmp %>%
  mutate(time = min(lumped_oc_data$start_date) + 7 * time) %>%
  ggplot(aes(time, value, ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  scale_fill_brewer() +
  # cowplot::theme_minimal_hgrid() +
  ylab("Rt or Re? idk")

# Really sus


# Posterior Predictive ----------------------------------------------------
stem_object$dynamics$parameters
stemr_params <- map(multi_chain_stem_fit$stem_fit_list,
                            function(stem_fit) {
                              cbind(
                                stem_fit$results$posterior$parameter_samples_nat,
                                stem_fit$results$posterior$initdist_samples) %>%
                                split(., row(.)) %>%
                                map(~`names<-`(.,c("log_R0_init", "gamma", "nu_early", "mu_rec", "mu_death", "rho_death", "phi_death", "alpha0", "alpha1", "kappa", "sigma", "S_0", "E_0", "Ie_0", "Ip_0", "R_0", "D_0")))
                            }) %>%
  unlist(recursive = F)

# post processing is a bit confusing since things aren't named


multi_chain_stem_fit$stem_fit_list[[1]]$results$posterior$tparam_draws
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

debugonce(simulate_stem)
stem_object$dynamics$fixed_inits <- T

tmp_pp <- simulate_stem(stem_object = stem_object,
                     nsim = 12000,
                     simulation_parameters = stemr_params,
                     tparam_draws = stemr_tparams,
                     method = "ode", paths = F)

tmp_pp2 <- stemr::simulate_stem(stem_object = stem_object,
                               nsim = 12000,
                               simulation_parameters = stemr_params,
                               tparam_draws = stemr_tparams,
                               method = "ode", paths = F)

tmp_pp$natural_paths[[1]] == tmp_pp2$natural_paths[[1]]
tmp_pp$paths[[1]] == tmp_pp2$paths[[1]]



tmp <- map(1:12000, ~{
  stem_object$dynamics$parameters <- head(stemr_params[[1]], -6)
  stem_object$dynamics$initializer[[1]]$init_states <- tail(stemr_params[[1]], 6)

  simulate_stem(stem_object = stem_object,
                nsim = 1,
                tparam_draws = list(stemr_tparams[[1]]),
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



ggplot() +
  geom_lineribbon(data = tmp, mapping = aes(time, value, ymin = .lower, ymax = .upper)) +
  geom_point(data = dat %>%
               mutate(pos = cases / tests) %>%
               select(time, deaths, pos) %>%
               pivot_longer(-time), mapping = aes(time, value)) +
  facet_wrap(. ~ name, scales = "free_y") +
  scale_fill_brewer() +
  theme_bw() +
  ggtitle("Posterior Predictive & Real Data") +
  theme(legend.position = "none")

tmp %>% map("datasets") %>% map(pluck(1)) %>% do.call(rbind, .) %>% as_tibble()

# stem_object$dynamics$parameters <- head(stemr_params[[1]], -6)
# # stem_object$dynamics$initdist_params <- tail(stemr_params[[1]], 6)
# stem_object$dynamics$initializer[[1]]$init_states <- tail(stemr_params[[1]], 6)
#
# tmptmp <- simulate_stem(stem_object = stem_object,
#               nsim = 1,
#               # simulation_parameters = stemr_params,
#               tparam_draws = list(stemr_tparams[[1]]),
#               method = "ode", paths = F)
# tmptmp$paths
#
# tmptmp$natural_paths

multi_chain_stem_fit$stem_fit_list[[1]]$results$posterior$latent_paths[,,1]
pp_data <- do.call(rbind, tmp_pp$datasets) %>%
  as_tibble() %>%
  left_join(select(dat, time, tests)) %>%
  mutate(pos = cases / tests) %>%
  select(time, deaths, pos) %>%
  pivot_longer(-time) %>%
  group_by(time, name) %>%
  median_qi(.width = c(0.5, 0.8, 0.95))



ggplot() +
  geom_lineribbon(data = pp_data, mapping = aes(time, value, ymin = .lower, ymax = .upper)) +
  geom_point(data = dat %>%
               mutate(pos = cases / tests) %>%
               select(time, deaths, pos) %>%
               pivot_longer(-time), mapping = aes(time, value)) +
  facet_wrap(. ~ name, scales = "free_y") +
  scale_fill_brewer() +
  theme_bw() +
  ggtitle("Posterior Predictive & Real Data") +
  theme(legend.position = "none")



multi_chain_stem_fit$stem_fit_list[[1]]$results$posterior$parameter_samples_nat[1,]
