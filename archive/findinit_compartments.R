library(tidyverse)
library(stemr)
library(extraDistr)
library(foreach)
library(doRNG)
library(doParallel)
library(scales)
theme_set(cowplot::theme_cowplot())

# Read Data ---------------------------------------------------------------
# dat <- read_csv("oc_data_subset.csv")
# dat <- read_csv("~/Documents/stemr/.development_files/oc_covid19/oc_data_subset.csv")
dat <- read_rds("~/Documents/uci_covid_modeling/code/SEIeIpRD/oc/2020-08-12_2020-09-16/model_objects.rds")$lumped_ochca_covid
popsize <- 3175692L

dat %>%
  select(-start_date) %>%
  pivot_longer(-end_date) %>%
  ggplot(aes(x=end_date, y = value)) +
  geom_point() +
  geom_line() +
  facet_wrap(. ~ name, scales = "free_y") +
  labs(x = "Week", y = "Count", title = "Observed Incidence")

# Process old fit for initialization --------------------------------------
oc_post <- read_rds("~/Documents/uci_covid_modeling/code/SEIeIpRD/oc/2020-08-07_2020-09-11/oc_post.rds")
model_objects <- read_rds("~/Documents/uci_covid_modeling/code/SEIeIpRD/oc/2020-08-07_2020-09-11/model_objects.rds")
source("~/Documents/uci_covid_modeling/code/SEIeIpRD/plot_functions.R")

target <- prep_epi_curves(oc_post, model_objects) %>%
  filter(t == min(dat$start_date)) %>%
  ungroup() %>%
  select(-t) %>%
  filter(str_detect(g_text, "\\s", negate = T)) %>%
  pivot_wider(values_from = epi_curves, names_from = g_text) %>%
  mutate(across(-.draw, ~{round(. * model_objects$popsize)})) %>%
  mutate(R = model_objects$popsize - c(S + E + IE + IP + D)) %>%
  left_join(spread_draws(oc_post, IFR) %>% select(.draw, IFR)) %>%
  mutate(D = round(D + IFR * (popsize - model_objects$popsize)),
         R = round(R + (1 - IFR) * (popsize - model_objects$popsize))) %>%
  rename_all(str_to_title) %>%
  select(S, E, Ie, Ip, R, D)

init_states <- round(colMeans(target))

# Maximize likelihood
C <- optimize(f = function(C) sum(ddirmnom(x = target, size = popsize, alpha = init_states / C, log = T)), lower = 1, upper = 10000, maximum = T)$maximum

stemr_prior_samples <- rdirmnom(n = 8000, size = popsize, alpha = init_states / C) %>% `colnames<-`(names(init_states)) %>% as_tibble()

# Compare to Beta Prior ---------------------------------------------------
oc_prior <- read_rds("~/Documents/uci_covid_modeling/code/SEIeIpRD/oc/2020-08-12_2020-09-16/oc_prior.rds")

stan_prior_samples <- spread_draws(model = oc_prior, frac_carrs, frac_carrs_infec, frac_infec_early) %>%
  transmute(EIeIp = frac_carrs * model_objects$popsize,
            S = (1 - frac_carrs) * model_objects$popsize,
            IeIP = frac_carrs_infec * EIeIp,
            E = (1 - frac_carrs_infec) * EIeIp,
            Ie = frac_infec_early * IeIP,
            Ip = (1 - frac_infec_early) * IeIP) %>%
  select(S, E, Ie, Ip) %>%
  round() %>%
  as_tibble()

oc_post <- read_rds("~/Documents/uci_covid_modeling/code/SEIeIpRD/oc/2020-08-12_2020-09-16/oc_post.rds")

stan_post_samples <- spread_draws(model = oc_post, frac_carrs, frac_carrs_infec, frac_infec_early) %>%
  transmute(EIeIp = frac_carrs * model_objects$popsize,
         S = (1 - frac_carrs) * model_objects$popsize,
         IeIP = frac_carrs_infec * EIeIp,
         E = (1 - frac_carrs_infec) * EIeIp,
         Ie = frac_infec_early * IeIP,
         Ip = (1 - frac_infec_early) * IeIP) %>%
  select(S, E, Ie, Ip) %>%
  round() %>%
  as_tibble()

# And actual stemr --------------------------------------------------------

multi_chain_stem_fit <- read_rds("multi_chain_stem_fit.rds")

stemr_post_samples <- map_dfr(multi_chain_stem_fit$stem_fit_list, ~as_tibble(.$results$posterior$initdist_samples)) %>%
  rename_all(~str_remove(., "_0"))

all_samples <- bind_rows(mutate(stemr_post_samples, dist = "post", model = "stemr"),
                         mutate(stemr_prior_samples, dist = "prior", model = "stemr"),
                         mutate(stan_post_samples, dist = "post", model = "stan"),
                         mutate(stan_prior_samples, dist = "prior", model = "stan")) %>%
  pivot_longer(-c(dist, model)) %>%
  mutate(name = fct_inorder(name))


all_samples %>%
  drop_na() %>%
  group_by(dist, model, name) %>%
  median_qi(.width = c(0.5, 0.8, 0.95)) %>%
  unite(model, dist, col = "source", remove = F, sep = "\n") %>%
  mutate(source = fct_reorder2(source, model, dist)) %>%
  ggplot(aes(source, value, ymin = .lower, ymax =.upper, color = source)) +
  facet_wrap(. ~ name, scales = "free_y") +
  geom_pointinterval() +
  cowplot::theme_minimal_hgrid() +
  ggtitle("Initial Conditions for Stan and Stemr") +
  theme(legend.position = "none") +
  scale_y_continuous(labels = comma) +
  ylab("People")


all_samples %>%
  drop_na() %>%
  ggplot(aes(x = model, y = value, fill = model)) +
  # facet_grid(dist ~ name, scales = "free") +
  facet_grid(name ~ dist, scales = "free") +
  ggdist::stat_eye(normalize = "xy") +
  theme_bw() +
  cowplot::theme_minimal_hgrid() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Initial Conditions for Stan and Stemr") +
  theme(legend.position = "bottom")


all_samples %>%
  drop_na() %>%
  group_by(dist, model, name) %>%
  median_qi(.width = c(0.5, 0.8, 0.95)) %>%
  ggplot(aes(value, model, xmin = .lower, xmax =.upper, color = dist)) +
  facet_grid(dist ~ name, scales = "free") +
  geom_pointinterval() +
  cowplot::theme_minimal_vgrid()

stemr_post_samples %>% mutate(.iteration = row_number()) %>%
  pivot_longer(-.iteration) %>%
  ggplot(aes(.iteration, value)) +
  facet_wrap(. ~ name, scales = "free_y") +
  geom_line()
