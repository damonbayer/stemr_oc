library(tidyverse)
library(tidybayes)
library(fs)
library(cowplot)
library(scales)
library(stemr)
library(coda)
library(extraDistr)
library(lubridate)
source('stemr_functions.R')
dat <- read_rds("final_models/data/dat.rds")
all_forecast_results <- read_rds("final_models/forecast_results/forecast_results.rds") %>%
  filter(date <= max(dat$end_date)) %>%
  mutate(weeks_ahead = as.numeric(date - model_name, "weeks")) %>%
  mutate(weeks_ahead_label = str_c(weeks_ahead, if_else(weeks_ahead == 1, "Week", "Weeks"), "Ahead Forecast", sep = " "))
multi_chain_stem_fit <- read_rds("final_models/results/2021-02-28.rds")

popsize <- multi_chain_stem_fit$stem_fit_list[[1]]$dynamics$popsize
ci_width <- c(0.95)

epi_curves_p <- extract_epi_curves(multi_chain_stem_fit = multi_chain_stem_fit, curve_type = "p", tidy = T)
stemr_all_params <- read_rds("final_models/results/stemr_all_params/2021-02-28_stemr_all_params.rds")
pp_data <- read_rds("final_models/results/pp_data/2021-02-28_pp_data.rds")

pp_data_intervals <-
  pp_data %>%
  select(.chain, date, name, value) %>%
  group_by(.chain, date, name) %>%
  median_qi(.width = ci_width)

all_pp_intervals <-
  bind_rows(
    pp_data_intervals %>%
      mutate(type = "Posterior Predictive",
             weeks_ahead = 0,
             weeks_ahead_label = "Posterior Predictive")#,
    # all_forecast_results %>%
    #   mutate(type = "Posterior Forecast") %>%
    #   select(-model_name)
    ) %>%
  ungroup() %>%
  filter(weeks_ahead %in% c(0, 1, 4)) %>%
  mutate(type = fct_inorder(type),
         weeks_ahead_label = fct_inorder(weeks_ahead_label),
         .width = .width %>% percent() %>% fct_inorder()) %>%
  arrange(date, type, .width)


all_dat <-
  all_pp_intervals %>%
  select(date, weeks_ahead_label) %>%
  distinct() %>%
  left_join(dat %>%
              mutate(pos = cases / tests) %>%
              select(date = end_date, pos, deaths) %>%
              pivot_longer(-date))

deaths_forecast_plot_chain <-
  ggplot() +
  geom_lineribbon(
    data = all_pp_intervals %>%
      mutate(.chain = as_factor(.chain)) %>%
      filter(name == "deaths"),
    mapping = aes(date, value, ymin = .lower, ymax = .upper, group = .chain, fill = .chain),
    alpha = 0.1,
    show.legend = T) +
  geom_point(size = 0.5, data = all_dat %>% filter(name == "deaths"),
             mapping = aes(date, value)) +
  scale_y_continuous(name = "Deaths", labels = comma) +
  scale_x_date(name = "Date", date_breaks = "3 months", date_labels = "%b '%y") +
  scale_fill_discrete(name = "Chain") +
  theme_minimal_grid() +
  theme(legend.position = c(.05, .8),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white")) +
ggtitle("Posterior Predictive Deaths",
        subtitle = "95% Confidence Intervals")


pos_forecast_plot_chain <-
  ggplot() +
  geom_lineribbon(
    data = all_pp_intervals %>%
      mutate(.chain = as_factor(.chain)) %>%
      filter(name == "pos"),
    mapping = aes(date, value, ymin = .lower, ymax = .upper, group = .chain, fill = .chain),
    alpha = 0.1,
    show.legend = T) +
  geom_point(size = 0.5, data = all_dat %>% filter(name == "pos"),
             mapping = aes(date, value)) +
  scale_y_continuous(name = "Testing Positivity", labels = function(.) scales::percent(., accuracy = 1)) +
  scale_x_date(name = "Date", date_breaks = "3 months", date_labels = "%b '%y") +
  scale_fill_discrete(name = "Chain") +
  theme_minimal_grid() +
  theme(legend.position = "none") +
  ggtitle("Posterior Predictive Testing Positivity",
          subtitle = "95% Confidence Intervals")



pp_forecasts_plot_chain <- plot_grid(deaths_forecast_plot_chain, pos_forecast_plot_chain, align = "hv", nrow = 2, ncol = 1)

save_plot(pp_forecasts_plot_chain,
          filename = "~/Documents/oc_covid19_stemr_manuscript/figures/pp_forecast_plot_chain.pdf",
          nrow = 2,
          ncol = 1,
          base_asp = 2)

