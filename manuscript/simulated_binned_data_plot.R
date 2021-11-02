library(tidyverse)
library(scales)
library(patchwork)
library(cowplot)
theme_set(theme_set(theme_minimal_grid()))


date_breaks <- "3 months"
date_labels <- "%b '%y"

dat <- read_csv("final_models/data/simulated_data.csv")
multi_chain_stem_fit <- read_rds("final_models/results/simulated_data.rds")

binned_data_plot <-
  ggplot(dat, aes(end_date, tests)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(name = "Tests", labels = comma) +
  scale_x_date(name = "Date", date_breaks = date_breaks, date_labels = date_labels) +
  ggplot(multi_chain_stem_fit$data, aes(end_date, cases)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(name = "Cases", labels = comma) +
  scale_x_date(name = "Date", date_breaks, date_labels = date_labels) +
  ggplot(multi_chain_stem_fit$data, aes(end_date, deaths)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(name = "Deaths", labels = comma) +
  scale_x_date(name = "Date", date_breaks, date_labels = date_labels) +
  ggplot(multi_chain_stem_fit$data, aes(end_date, cases / tests)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(name = "Testing Positivity",
                     labels = function(.) scales::percent(., accuracy = 1),
                     limits = c(0, NA)) +
  scale_x_date(name = "Date", date_breaks, date_labels = date_labels) +
  patchwork::plot_layout(ncol = 2, nrow = 2) +
  patchwork::plot_annotation(title = "Simulated Data")

save_plot(filename = "~/Documents/oc_covid19_stemr_manuscript/figures/simulated_binned_data_plot.pdf",
          plot = binned_data_plot,
          ncol = 2,
          nrow = 2,
          base_height = 2.5)
