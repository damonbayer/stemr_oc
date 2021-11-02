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
ci_width <- 0.95
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

d_i_p_curves_chain <-
  epi_curves_p %>%
  mutate(`Cumulative Incidence` = popsize - S,
         `Prevalence` = E + I,
         `Cumulative Deaths` = D) %>%
  left_join(time_date_conversion) %>%
  select(-c(time, S, E, I, R, D)) %>%
  select(-.iteration, -.draw) %>%
  pivot_longer(-c(date, .chain)) %>%
  group_by(date, name, .chain) %>%
  median_qi(.width = 0.95) %>%
  ungroup() %>%
  mutate(.chain = as_factor(.chain))


cumulative_deaths_plot_chain <-
  ggplot() +
  geom_lineribbon(data = d_i_p_curves_chain %>%
                    filter(name == "Cumulative Deaths"),
                  mapping = aes(date, value, ymin = .lower, ymax = .upper, fill = .chain),
                  alpha = 0.1) +
  # scale_fill_brewer(name = "Credibility Level") +
  scale_color_brewer() +
  scale_y_continuous(name = "Deaths", labels = comma) +
  scale_x_date(name = "Date", date_breaks = date_breaks3, date_labels = date_labels) +
  # theme(legend.position = "none") +
  ggtitle("Posterior Latent Deaths",
          subtitle = "95% Credible Intervals") +
  scale_fill_discrete("Chain")

cumulative_incidence_plot_chain <-
  ggplot() +
  geom_lineribbon(data = d_i_p_curves_chain %>%
                    filter(name == "Cumulative Incidence"),
                  mapping = aes(date, value, ymin = .lower, ymax = .upper, fill = .chain),
                  alpha = 0.1) +
  # scale_fill_brewer(name = "Credibility Level") +
  scale_color_brewer() +
  scale_y_continuous(name = "Cumulative Incidence", labels = comma) +
  scale_x_date(name = "Date", date_breaks = date_breaks3, date_labels = date_labels) +
  theme(legend.position = "none") +
  ggtitle("Posterior Cumulative Incidence",
          subtitle = "95% Credible Interval") +
  scale_fill_discrete("Chain")

prevalence_plot_chain <-
  ggplot() +
  geom_lineribbon(data = d_i_p_curves_chain %>%
                    filter(name == "Prevalence"),
                  mapping = aes(date, value, ymin = .lower, ymax = .upper, fill = .chain),
                  alpha = 0.1) +
  scale_color_brewer() +
  scale_y_continuous(name = "Prevalence", labels = comma) +
  scale_x_date(name = "Date", date_breaks = date_breaks3, date_labels = date_labels) +
  theme(legend.position = "none") +
  ggtitle("Posterior Cumulative Incidence",
          subtitle = "95% Credible Intervals") +
  scale_fill_discrete("Chain")


dip_plot_chain <- plot_grid(cumulative_deaths_plot_chain + theme(legend.position = c(0.1, 0.5)) + guides(fill=guide_legend(ncol=2)),
                            cumulative_incidence_plot_chain,
                            prevalence_plot_chain,
                            align = "hv", ncol = 3, nrow = 1)

save_plot(filename = "~/Documents/oc_covid19_stemr_manuscript/figures/dip_plot_chain.pdf",
          plot = dip_plot_chain,
          ncol = 3,
          nrow = 1,
          base_asp = 1.25)
