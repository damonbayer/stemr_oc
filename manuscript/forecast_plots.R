library(tidyverse)
library(tidybayes)
library(fs)
library(cowplot)
library(scales)
all_forecast_results <- read_rds("final_models/forecast_results/forecast_results.rds")
dat <- read_rds("final_models/data/dat.rds")

weeks_ahead_key <-
  all_forecast_results %>%
  select(date, model_name) %>%
  distinct() %>%
  mutate(weeks_ahead = as.numeric(date - model_name, "weeks")) %>%
  select(date, weeks_ahead) %>%
  mutate(weeks_ahead_label = str_c(weeks_ahead, if_else(weeks_ahead == 1, "week", "weeks"), sep = " "))


deaths_forecast_plot <-
  ggplot() +
  geom_interval(
    data = all_forecast_results %>%
      filter(name == "deaths"),
    mapping = aes(date, value, ymin = .lower, ymax = .upper),
    show_point = T,
    shape = "_",
    point_color = "black",
    fatten_point = 1,
    show.legend = F) +
  geom_point(
    data = dat %>%
      mutate(pos = cases / tests) %>%
      select(date = end_date, pos, deaths) %>%
      pivot_longer(-date) %>%
      filter(name == "deaths") %>%
      left_join(weeks_ahead_key) %>%
      drop_na(),
    mapping = aes(date, value, shape = weeks_ahead_label)) +
  scale_shape(name = "Forecast\nHorizon") +
  scale_y_continuous(name = "Deaths", labels = comma) +
  scale_x_date(name = "Date") +
  scale_color_brewer() +
  theme_minimal_hgrid() +
  ggtitle("Posterior Deaths Forecasts")

pos_forecast_plot <-
  ggplot() +
  geom_interval(
    data = all_forecast_results %>%
      filter(name == "pos"),
    mapping = aes(date, value, ymin = .lower, ymax = .upper),
    show_point = T,
    shape = "_",
    point_color = "black",
    fatten_point = 1,
    show.legend = F) +
  geom_point(
    data = dat %>%
      mutate(pos = cases / tests) %>%
      select(date = end_date, pos, deaths) %>%
      pivot_longer(-date) %>%
      filter(name == "pos") %>%
      left_join(weeks_ahead_key) %>%
      drop_na(),
    mapping = aes(date, value, shape = weeks_ahead_label)) +
  scale_shape(name = "Forecast\nHorizon") +
  scale_y_continuous(name = "Testing Positivity", labels = percent) +
  scale_x_date(name = "Date") +
  scale_color_brewer() +
  theme_minimal_hgrid() +
  ggtitle("Posterior Testing Positivity Forecasts")

save_plot("~/Documents/oc_covid19_stemr_manuscript/figures/deaths_forecast_plot.pdf", deaths_forecast_plot)
save_plot("~/Documents/oc_covid19_stemr_manuscript/figures/pos_forecast_plot.pdf", pos_forecast_plot)
