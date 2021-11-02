library(tidyverse)
library(tidybayes)
library(scales)
library(cowplot)
library(fs)
library(stemr)
library(coda)
library(lubridate)
library(glue)
source('stemr_functions.R')
theme_set(cowplot::theme_minimal_grid())

# Plot Options
ci_width <- c(0.5, 0.8, 0.95)
date_breaks <- "2 month"
date_breaks3 <- "3 months"
date_labels <- "%b '%y"


d_i_p_curves_city <-
  stemr_all_time_varying_params_city %>%
  select(date,
         model,
         `Cumulative Incidence` = cumulative_incidence,
         Prevalence = prevalence,
         `Cumulative Deaths` = cumulative_latent_deaths) %>%
  pivot_longer(-c(date, model)) %>%
  group_by(date, model, name) %>%
  median_qi(.width = ci_width)

d_i_data_city <-
  stemr_all_time_varying_params_city %>%
  select(date,
         model,
         `Cumulative Incidence` = cumulative_cases,
         `Cumulative Deaths` = cumulative_deaths) %>%
  distinct() %>%
  pivot_longer(-c(date, model))

d_i_p_curves_sensitivity <-
  stemr_all_time_varying_params_sensitivity %>%
    select(date,
           model,
           `Cumulative Incidence` = cumulative_incidence,
           Prevalence = prevalence,
           `Cumulative Deaths` = cumulative_latent_deaths) %>%
    pivot_longer(-c(date, model)) %>%
    group_by(date, model, name) %>%
    median_qi(.width = ci_width)

d_i_data_sensitivity <-
  stemr_all_time_varying_params_sensitivity %>%
  select(date,
         model,
         `Cumulative Incidence` = cumulative_cases,
         `Cumulative Deaths` = cumulative_deaths) %>%
  distinct() %>%
  pivot_longer(-c(date, model))


d_i_p_city_plot <-
  ggplot() +
  facet_wrap(model ~ name, scales = "free_y", ncol = 3) +
  geom_lineribbon(data = d_i_p_curves_city,
                  mapping = aes(date, value, ymin = .lower, ymax = .upper),
                  color = "steelblue4",
                  show.legend = F) +
  geom_line(data = d_i_data_city,
            mapping = aes(date, value)) +
  geom_point(data = d_i_data_city,
             mapping = aes(date, value)) +
  scale_fill_brewer(name = "Credibility Level") +
  scale_y_continuous(name = "People", labels = comma) +
  scale_x_date(name = "Date", date_breaks = date_breaks3, date_labels = date_labels)


d_i_p_sensitivity_plot <-
  ggplot() +
  facet_wrap(model ~ name, scales = "free_y", ncol = 3) +
  geom_lineribbon(data = d_i_p_curves_sensitivity,
                  mapping = aes(date, value, ymin = .lower, ymax = .upper),
                  color = "steelblue4",
                  show.legend = F) +
  geom_line(data = d_i_data_sensitivity,
            mapping = aes(date, value)) +
  geom_point(data = d_i_data_sensitivity,
             mapping = aes(date, value)) +
  scale_fill_brewer(name = "Credibility Level") +
  scale_y_continuous(name = "People", labels = comma) +
  scale_x_date(name = "Date", date_breaks = "4 months", date_labels = date_labels)




save_plot(filename = "~/Documents/oc_covid19_stemr_manuscript/figures/d_i_p_sensitivity_plot.pdf",
          plot = d_i_p_sensitivity_plot,
          ncol = 3,
          nrow = 6,
          base_height = 1.75,
          base_asp = (8.5 * 6) / (11 * 3))


# save_plot(filename = path("~/Documents/stemr_oc/manuscript/city_sensitivity_plots/", "d_i_p_city_plot", ext = "pdf"),
#           plot = d_i_p_city_plot,
#           ncol = 3, nrow = 4,
#           # base_asp = (16 * 4) / (9 * 3))
#           base_asp = (8.5 * 4) / (11 * 3))
#
# save_plot(filename = path("~/Documents/stemr_oc/manuscript/city_sensitivity_plots/", "d_i_p_sensitivity_plot", ext = "pdf"),
#           plot = d_i_p_sensitivity_plot,
#           ncol = 3, nrow = 6,
#           # base_asp = (16 * 6) / (9 * 3))
#           base_asp = (8.5 * 6) / (11 * 3))
