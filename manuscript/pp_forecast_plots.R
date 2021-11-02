library(tidyverse)
library(tidybayes)
library(fs)
library(cowplot)
library(scales)
library(stemr)
library(coda)
library(extraDistr)
source('stemr_functions.R')
all_forecast_results <- read_rds("final_models/forecast_results/forecast_results.rds")
dat <- read_rds("final_models/data/dat.rds")
multi_chain_stem_fit <- read_rds("final_models/results/2021-02-28.rds")

popsize <- multi_chain_stem_fit$stem_fit_list[[1]]$dynamics$popsize
ci_width <- c(0.5, 0.8, 0.95)

epi_curves_p <- extract_epi_curves(multi_chain_stem_fit = multi_chain_stem_fit, curve_type = "p", tidy = T)
stemr_all_params <- read_rds("final_models/results/stemr_all_params/2021-02-28_stemr_all_params.rds")
pp_data <- read_rds("final_models/results/pp_data/2021-02-28_pp_data.rds")

pp_data_intervals <-
  pp_data %>%
  group_by(date, name) %>%
  median_qi(.width = ci_width)

weeks_ahead_key <-
  all_forecast_results %>%
  select(date, model_name) %>%
  distinct() %>%
  mutate(weeks_ahead = as.numeric(date - model_name, "weeks")) %>%
  select(date, weeks_ahead) %>%
  mutate(weeks_ahead_label = str_c(weeks_ahead, if_else(weeks_ahead == 1, "week", "weeks"), sep = " "))

all_pp_intervals <-
  bind_rows(
    pp_data_intervals %>%
      mutate(type = "Posterior Predictive",
             weeks_ahead = 0,
             weeks_ahead_label = "0 weeks"),
    all_forecast_results %>%
      mutate(type = "Posterior Forecast") %>%
      left_join(weeks_ahead_key) %>%
      select(-model_name)) %>%
  ungroup() %>%
  mutate(type = fct_inorder(type),
         .width = .width %>% percent() %>% fct_inorder()) %>%
  arrange(date, type, .width)

tmp <- dat %>%
  mutate(pos = cases / tests) %>%
  select(date = end_date, pos, deaths) %>%
  pivot_longer(-date)

all_dat <-
  rbind(
    mutate(tmp,
           weeks_ahead_label = "0 weeks",
           type = "Posterior Predictive"),
    left_join(tmp, weeks_ahead_key) %>%
      mutate(type = "Posterior Forecast") %>%
      select(-weeks_ahead) %>%
      drop_na()) %>%
  mutate(type = fct_inorder(type))

deaths_pp_forecast_plot <-
  ggplot() +
  facet_wrap(. ~ type,
             nrow = 2,
             scales = "free_y",
             labeller = as_labeller(function(.) paste0(., " Deaths"))) +
  geom_interval(
    data = all_pp_intervals %>%
      filter(name == "deaths"),
    mapping = aes(date, value, ymin = .lower, ymax = .upper),
    show_point = F,
    shape = "_",
    point_color = "black",
    fatten_point = 1,
    show.legend = T) +
  geom_point(data = all_dat %>% filter(name == "deaths"),
             mapping = aes(date, value, shape = weeks_ahead_label)) +
  scale_y_continuous(name = "Deaths", labels = comma) +
  scale_x_date(name = "Date", date_breaks = "3 months", date_labels = "%b '%y") +
  scale_shape_discrete(name = "Forecast\nHorizon") +
  # scale_color_brewer(name = "Credibility", labels = function(.) percent(as.numeric(.))) +
  scale_color_brewer(name = "Credibility") +
  theme_minimal_grid() +
  # ggtitle("Posterior Predictive Deaths & Forecast") +
  theme(legend.position = c(.05, .25),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white")) +
  guides(shape = guide_legend(oreder = 1, ncol = 2),
         color = guide_legend(order = 2, reverse = T))



pos_pp_forecast_plot <-
  ggplot() +
  facet_wrap(. ~ type,
             nrow = 2,
             scales = "free_y",
             labeller = as_labeller(function(.) paste0(., " Testing Positivity"))) +
  geom_interval(
    data = all_pp_intervals %>%
      filter(name == "pos"),
    mapping = aes(date, value, ymin = .lower, ymax = .upper),
    show_point = F,
    shape = "_",
    point_color = "black",
    fatten_point = 1,
    show.legend = T) +
  geom_point(data = all_dat %>% filter(name == "pos"),
             mapping = aes(date, value, shape = weeks_ahead_label)) +
  scale_y_continuous(name = "Testing Positivity", labels = function(.) scales::percent(., accuracy = 1)) +
  scale_x_date(name = "Date", date_breaks = "3 months", date_labels = "%b '%y") +
  scale_shape_discrete(name = "Forecast\nHorizon") +
  # scale_color_brewer(name = "Credibility", labels = function(.) percent(as.numeric(.))) +
  scale_color_brewer(name = "Credibility") +
  theme_minimal_grid() +
  # ggtitle("Posterior Predictive Testing Positivity & Forecast") +
  theme(legend.position = "none")


pp_forecast_plot <-
  cowplot::plot_grid(deaths_pp_forecast_plot, pos_pp_forecast_plot, align = "hv", axis = "l", ncol = 2) +
  patchwork::plot_annotation(title = "Posterior Predictive & Forecasts")

save_plot(filename = "~/Documents/oc_covid19_stemr_manuscript/figures/pp_forecast_plot.pdf",
          plot = pp_forecast_plot,
          ncol = 2, nrow = 2, base_asp = 7/5)


write_rds(x = deaths_pp_forecast_plot, "ISBA_2021/plots_rds/deaths_pp_forecast_plot.rds")
write_rds(x = pos_pp_forecast_plot +
            theme(legend.position = c(.05, .25),
                  legend.direction = "horizontal",
                  legend.background = element_rect(fill = "white")) +
            guides(shape = guide_legend(oreder = 1, ncol = 2),
                   color = guide_legend(order = 2, reverse = T)),
          "ISBA_2021/plots_rds/pos_pp_forecast_plot.rds")
