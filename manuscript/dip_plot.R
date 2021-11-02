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
date_labels <- "%b %y"

multi_chain_stem_fit <- read_rds("final_models/results/2021-02-28.rds")
time_date_conversion <- select(multi_chain_stem_fit$data, time, date = end_date) %>%
  add_row(., time = 0, date = min(.$date) - 7, .before = T)
popsize <- multi_chain_stem_fit$stem_fit_list[[1]]$dynamics$popsize
log_popsize <- log(popsize)

# Gather Parameters -------------------------------------------------------
epi_curves_p <- extract_epi_curves(multi_chain_stem_fit = multi_chain_stem_fit, curve_type = "p", tidy = F)
stemr_all_params <- read_rds("final_models/results/stemr_all_params/2021-02-28_stemr_all_params.rds")

# DIP Plot ----------------------------------------------------------------
d_i_data <-
  multi_chain_stem_fit$data %>%
  mutate(`Cumulative Incidence` = cumsum(cases),
         `Cumulative Deaths` = cumsum(deaths)) %>%
  select(date = end_date, starts_with("Cumulative")) %>%
  pivot_longer(-date)

d_i_p_curves <- epi_curves_p %>%
  mutate(`Cumulative Incidence` = popsize - S,
         `Prevalence` = E + I,
         `Cumulative Deaths` = D) %>%
  left_join(time_date_conversion) %>%
  select(-c(time, S, E, I, R, D)) %>%
  select(-starts_with(".")) %>%
  pivot_longer(-date) %>%
  group_by(date, name) %>%
  median_qi(.width = ci_width) %>%
  mutate(.width = percent(.width))

cumulative_deaths_plot <-
  ggplot() +
  # geom_interval(data = d_i_p_curves %>%
  #                 filter(name == "Cumulative Deaths"),
  #               mapping = aes(date, value, ymin = .lower, ymax = .upper),
  #               show_point = T,
  #               shape = "_",
  #               point_color = "black",
  #               fatten_point = 1) +
  geom_lineribbon(data = d_i_p_curves %>%
                    filter(name == "Cumulative Deaths"),
                  mapping = aes(date, value, ymin = .lower, ymax = .upper),
                  color = "steelblue4") +
  geom_line(data = d_i_data %>%
              filter(name == "Cumulative Deaths"),
            mapping = aes(date, value)) +
  geom_point(data = d_i_data %>%
               filter(name == "Cumulative Deaths"),
             mapping = aes(date, value)) +
  annotate(geom = "curve",
           x = lubridate::ymd("2020-09-01"),
           xend = lubridate::ymd("2021-02-07") - days(4),
           y = 3500,
           yend = d_i_p_curves %>%
             filter(date == "2021-02-07",
                    name == "Cumulative Deaths",
                    .width == "95%") %>%
             pull(value),
           curvature = -0.25,
           arrow = arrow(length = unit(4, "mm"))) +
  annotate(geom = "text",
           x = lubridate::ymd("2020-09-01"),
           y = 3500,
           label = "Posterior",
           hjust = "center",
           vjust = "top") +
  annotate(geom = "curve",
           x = lubridate::ymd("2020-12-31"),
           xend = lubridate::ymd("2021-02-07"),
           y = 500,
           yend = d_i_data %>%
             filter(date == "2021-02-07",
                    name == "Cumulative Deaths") %>%
             pull(value) %>%
             `-`(100),
           curvature = 0.25,
           arrow = arrow(length = unit(4, "mm"))) +
  annotate(geom = "text",
           x = lubridate::ymd("2020-12-31"),
           y = 500,
           label = "Observed Data", hjust = "right", vjust = "center") +
  scale_fill_brewer(name = "Credibility Level") +
  scale_color_brewer() +
  scale_y_continuous(name = "Deaths", labels = comma) +
  scale_x_date(name = "Date", date_breaks = date_breaks3, date_labels = date_labels) +
  theme(legend.position = "none") +
  ggtitle("Posterior Latent and\nObserved Deaths")


seroprev_interval <-
  d_i_p_curves %>%
  filter(name == "Cumulative Incidence",
         date == "2020-08-16",
         .width == "95%") %>%
  ungroup() %>%
  select(value, .lower, .upper) %>%
  pivot_longer(everything()) %>%
  mutate(value = scales::percent(value / popsize)) %>%
  deframe()


cumulative_incidence_plot <-
  ggplot() +
  # geom_interval(data = d_i_p_curves %>%
  #                 filter(name == "Cumulative Incidence"),
  #               mapping = aes(date, value, ymin = .lower, ymax = .upper),
  #               show_point = T,
  #               shape = "_",
  #               point_color = "black",
  #               fatten_point = 1) +
  geom_lineribbon(data = d_i_p_curves %>%
                  filter(name == "Cumulative Incidence"),
                  mapping = aes(date, value, ymin = .lower, ymax = .upper),
                color = "steelblue4") +
  geom_line(data = d_i_data %>%
              filter(name == "Cumulative Incidence"),
            mapping = aes(date, value)) +
  geom_point(data = d_i_data %>%
               filter(name == "Cumulative Incidence"),
             mapping = aes(date, value)) +
  annotate(geom = "curve",
           x = lubridate::ymd("2020-07-01"),
           xend = lubridate::ymd("2020-08-16"),
           y = 1250000,
           yend = d_i_p_curves %>%
             filter(date == "2020-08-16",
                    name == "Cumulative Incidence",
                    .width == "95%") %>%
             pull(.upper),
           curvature = 0.25,
           # arrow = arrow = arrow(length = unit(4, "mm"))) +
           arrow = arrow(length = unit(4, "mm"))) +
  annotate(geom = "text",
           x = lubridate::ymd("2020-07-01"),
           y = 1250000 + 50000,
           label = glue("Posterior\nSeroprevalence:\n{seroprev_interval[['value']]}\n({seroprev_interval[['.lower']]}, {seroprev_interval[['.upper']]})"),
           hjust = "center",
           vjust = "bottom") +
  scale_fill_brewer(name = "Credibility Level") +
  scale_color_brewer() +
  scale_y_continuous(name = "Incidence", labels = comma) +
  scale_x_date(name = "Date", date_breaks = date_breaks3, date_labels = date_labels) +
  theme(legend.position = "none") +
  ggtitle("Posterior Latent and Observed\nCumulative Incidence")

prevalence_plot <-
  ggplot() +
  # geom_interval(
  #   data = d_i_p_curves %>%
  #     filter(name == "Prevalence"),
  #   mapping = aes(date, value, ymin = .lower, ymax = .upper),
  #   show_point = T,
  #   shape = "_",
  #   point_color = "black",
  #   fatten_point = 1) +
  geom_lineribbon(
    data = d_i_p_curves %>%
      filter(name == "Prevalence"),
    mapping = aes(date, value, ymin = .lower, ymax = .upper),
    color = "steelblue4") +
  scale_fill_brewer(name = "Credibility Level") +
  scale_color_brewer() +
  scale_y_continuous(name = "Prevalence", labels = comma) +
  scale_x_date(name = "Date", date_breaks = date_breaks3, date_labels = date_labels) +
  theme(legend.position = c(.1, .7),
        legend.background = element_rect(fill = "white") ) +
  ggtitle("Posterior Latent\nPrevalence") +
  guides(color = guide_legend(title = "Credibility", reverse = T))

dip_plot <- plot_grid(cumulative_deaths_plot,
                      cumulative_incidence_plot,
                      prevalence_plot,
                      align = "hv", ncol = 3, nrow = 1)

save_plot(filename = "~/Documents/oc_covid19_stemr_manuscript/figures/dip_plot.pdf",
          plot = dip_plot,
          ncol = 3,
          nrow = 1,
          base_asp = 1.25)

write_rds(dip_plot, "ISBA_2021/plots_rds/dip_plot.rds")
