library(tidyverse)
library(tidybayes)
library(scales)
library(stemr)
source("~/Documents/uci_covid_modeling/code/SEIeIpRD/plot_functions.R")
source('stemr_functions.R')
# multi_chain_stem_fit <- read_rds("multi_chain_stem_fit.rds")
# stan_oc_post <- read_rds("~/Documents/uci_covid_modeling/code/SEIeIpRD/oc/2020-08-12_2020-09-16/oc_post.rds")
# stan_oc_model_objects <- read_rds("~/Documents/uci_covid_modeling/code/SEIeIpRD/oc/2020-08-12_2020-09-16/model_objects.rds")
theme_set(cowplot::theme_cowplot())
# multi_chain_stem_fit <- read_rds("fixed_init_sim/multi_chain_stem_fit.rds")
multi_chain_stem_fit <- read_rds("fixed_init_sim/multi_chain_stem_fit_R0D0.rds")
stan_oc_post <- read_rds("fixed_init_sim/oc_post.rds")
stan_oc_model_objects <- read_rds("fixed_init_sim/model_objects.rds")


# Prevalence --------------------------------------------------------------
stemr_epi_curves <- extract_epi_curves(multi_chain_stem_fit, curve_type = "p", tidy = T) %>%
  filter(name %in% c("S", "E", "Ie", "Ip")) %>%
  select(time, name, value)

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
  mutate(name = fct_inorder(name))

all_epi_curves <- bind_rows(
  mutate(stemr_epi_curves, model = "stemr"),
  mutate(stan_epi_curves, model = "stan"))

all_epi_curves %>%
  group_by(name, model, time) %>%
  median_qi() %>%
  ggplot(aes(time, value, ymin = .lower, ymax = .upper, fill = model, color = model, group = model)) +
  facet_wrap(.~name, scale = "free_y") +
  # facet_grid(name ~ model, scale = "free_y") +
  geom_lineribbon(alpha = 0.5) +
  ylab("People") +
  scale_y_continuous(labels = comma) +
  ggtitle("Prevalence Curves for Stan and Stemr")



# Incidence ---------------------------------------------------------------
stemr_i_curves <- extract_epi_curves(multi_chain_stem_fit, curve_type = "i", tidy = T)
stan_i_curves <- prep_epi_curves(stan_oc_post, stan_oc_model_objects)

unique(stemr_i_curves$name)

stan_i_curves %>%
  mutate(g_text = str_to_title(g_text)) %>%
  pivot_wider(names_from = g_text, values_from = epi_curves,
              id_cols = c(.draw, t, g_text)) %>%
  mutate(R = 1 - (S + E + Ie + Ip + D)) %>%
  select(.draw, t, S, E, Ie, Ip, R, D)
  transmute(S2E = S - lag(S),
            E2Ie = E - lag(E),
            E2Ie2 = lag(Ie) - Ie)
