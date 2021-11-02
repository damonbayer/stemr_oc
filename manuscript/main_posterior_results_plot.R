library(tidyverse)
library(tidybayes)
library(scales)
library(cowplot)
library(fs)
library(stemr)
library(coda)
source('stemr_functions.R')
theme_set(cowplot::theme_minimal_hgrid())

# Plot Options
ci_width <- c(0.5, 0.8, 0.95)
date_breaks <- "2 month"
date_breaks3 <- "3 months"
date_labels <- "%b '%y"

multi_chain_stem_fit <-
  read_rds("final_models/results/2021-02-28.rds")

popsize <- multi_chain_stem_fit$stem_fit_list[[1]]$dynamics$popsize
log_popsize <- log(popsize)

time_date_conversion <-
  select(multi_chain_stem_fit$data, time, date = end_date) %>%
  add_row(.,
          time = 0,
          date = min(.$date) - 7,
          .before = T)

stemr_all_params <-
  read_rds("final_models/results/stemr_all_params/2021-02-28_stemr_all_params.rds")

weekly_case_ratio <-
  extract_epi_curves(
    multi_chain_stem_fit = multi_chain_stem_fit,
    curve_type = "i",
    tidy = F
  ) %>%
  filter(time != 0) %>%
  select(starts_with("."), time, latent_cases = "E2I") %>%
  left_join(multi_chain_stem_fit$data %>% select(time, date = end_date, cases)) %>%
  mutate(case_detection_rate = cases / latent_cases,
         latent_case_ratio = latent_cases / cases)

initial_cases <-
  read_csv("final_models/data/oc_data.csv") %>%
  filter(date < "2020-03-30") %>%
  pull(cases) %>%
  sum()


cumulative_case_ratio <-
  extract_epi_curves(
    multi_chain_stem_fit = multi_chain_stem_fit,
    curve_type = "p",
    tidy = F
  ) %>%
  filter(time != 0) %>%
  mutate(cumulative_latent_cases = popsize - S) %>%
  select(time, starts_with("."), cumulative_latent_cases) %>%
  left_join(
    multi_chain_stem_fit$data %>%
      mutate(cumulative_cases = cumsum(cases) + initial_cases) %>%
      select(time, date = end_date, cumulative_cases)
  ) %>%
  mutate(
    cumulative_case_detection_rate = cumulative_cases / cumulative_latent_cases,
    cumalitive_latent_case_ratio = cumulative_latent_cases / cumulative_cases
  )

# R0 Plot -----------------------------------------------------------------
R0_plot <-
  stemr_all_params %>%
  mutate(R0 = exp(log(beta_t) + log_popsize + log(dur_infec)),
         prop_s = S / popsize) %>%
  select(R0, date = end_date) %>%
  group_by(date) %>%
  median_qi(.width = ci_width) %>%
  drop_na() %>%
  mutate(.width = scales::percent(.width)) %>%
  ggplot(aes(date, R0, ymin = .lower, ymax = .upper)) +
  # geom_interval(
  #   show_point = T,
  #   shape = "_",
  #   point_color = "black",
  #   fatten_point = 1
  # ) +
  geom_lineribbon(color = "steelblue4") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_brewer("Credibility") +
  # scale_color_brewer() +
  scale_y_continuous(
    name = expression(R[0]),
    breaks = seq(0, 3, by = 0.25),
    limits = c(.5, 2.05)
  ) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks,
               date_labels = date_labels) +
  # ggtitle(expression(bold(paste("Posterior ", R[0]))))
  ggtitle("Posterior Basic Reproductive Number")


# Rt Plot -----------------------------------------------------------------
Rt_plot <-
  stemr_all_params %>%
  mutate(
    R0 = exp(log(beta_t) + log_popsize + log(dur_infec)),
    prop_s = S / popsize,
    Rt = R0 * prop_s
  ) %>%
  select(Rt, date = end_date) %>%
  group_by(date) %>%
  median_qi(.width = ci_width) %>%
  drop_na() %>%
  ggplot(aes(date, Rt, ymin = .lower, ymax = .upper)) +
  geom_lineribbon(color = "steelblue4") +
  # geom_interval(
  #   show_point = T,
  #   shape = "_",
  #   point_color = "black",
  #   fatten_point = 1
  # ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_brewer() +
  scale_color_brewer() +
  scale_y_continuous(
    name = expression(R[t]),
    breaks = seq(0, 3, by = 0.25),
    limits = c(.5, 2.05)
  ) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks,
               date_labels = date_labels) +
  theme(legend.position = "none") +
  # ggtitle(expression(bold(paste("Posterior ", R[t]))))
  ggtitle("Posterior Effective Reproductive Number")



# IFR Plot ----------------------------------------------------------------
ifr_t_plot <-
  stemr_all_params %>%
  select(date = end_date, ifr_t) %>%
  group_by(date) %>%
  median_qi(.width = ci_width) %>%
  drop_na() %>%
  ggplot(aes(date, ifr_t, ymin = .lower, ymax = .upper)) +
  geom_lineribbon(color = "steelblue4") +
  # geom_interval(
  #   show_point = T,
  #   shape = "_",
  #   point_color = "black",
  #   fatten_point = 1
  # ) +
  scale_fill_brewer() +
  # scale_color_brewer() +
  theme(legend.position = "none") +
  scale_y_continuous(
    name = "IFR",
    labels = function(.)
      scales::percent(., accuracy = 0.1),
    # trans = scales::logit_trans(),
    breaks = seq(0, 0.01, by = 0.001),
    limits = c(0, NA)
  ) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks,
               date_labels = date_labels) +
  ggtitle("Posterior Infection Fatality Ratio")


# Alpha_t plot ------------------------------------------------------------
alpha_t_plot <-
  stemr_all_params %>%
  select(date = end_date, alpha0_t) %>%
  group_by(date) %>%
  median_qi(.width = ci_width) %>%
  drop_na() %>%
  ggplot(aes(date, alpha0_t, ymin = .lower, ymax = .upper)) +
  geom_lineribbon(color = "steelblue4") +
  # geom_interval(
  #   show_point = T,
  #   shape = "_",
  #   point_color = "black",
  #   fatten_point = 1
  # ) +
  scale_fill_brewer() +
  # scale_color_brewer() +
  theme(legend.position = "none") +
  scale_y_continuous(name = expression(alpha),
                     # scale_y_continuous(name = expression(alpha[0]),
                     # trans = scales::log_trans(),
                     breaks = 0:5) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks,
               date_labels = date_labels) +
  ggtitle("Posterior \u03B1")



# Weekly Latent:Case Ratio Plot -------------------------------------------
weekly_latent_case_ratio_plot <-
  weekly_case_ratio %>%
  select(date, latent_case_ratio) %>%
  group_by(date) %>%
  median_qi(.width = ci_width) %>%
  ggplot(aes(date, latent_case_ratio, ymin = .lower, ymax = .upper)) +
  geom_lineribbon(color = "steelblue4", show.legend = F) +
  # geom_interval(
  #   show_point = T,
  #   shape = "_",
  #   point_color = "black",
  #   fatten_point = 1
  # ) +
  scale_fill_brewer() +
  # scale_color_brewer() +
  theme(legend.position = "none") +
  scale_y_continuous(
    name = "Weekly Latent:Case Ratio",
    labels = function(.)
      str_c(scales::comma(., accuracy = 1), ":1"),
    breaks = seq(0, 40, 5),
    limits = c(1, 27)
  ) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks,
               date_labels = date_labels) +
  ggtitle("Posterior Weekly Latent:Case Ratio")


# Cumulative Latent:Case Ratio Plot ---------------------------------------
cummulative_latent_case_ratio_plot <-
  cumulative_case_ratio %>%
  select(date, cumalitive_latent_case_ratio) %>%
  group_by(date) %>%
  median_qi(.width = ci_width) %>%
  ggplot(aes(
    date,
    cumalitive_latent_case_ratio,
    ymin = .lower,
    ymax = .upper
  )) +
  geom_lineribbon(color = "steelblue4", show.legend = F) +
  # geom_interval(
  #   show_point = T,
  #   shape = "_",
  #   point_color = "black",
  #   fatten_point = 1
  # ) +
  scale_fill_brewer() +
  # scale_color_brewer() +
  theme(legend.position = "none") +
  scale_y_continuous(
    name = "Cumulative Latent:Case Ratio",
    labels = function(.)
      str_c(scales::comma(., accuracy = 1), ":1"),
    breaks = seq(0, 100, 5),
    limits = c(1, 27)
  ) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks,
               date_labels = date_labels) +
  ggtitle("Posterior Cumulative Latent:Case Ratio")



# Save Plot ---------------------------------------------------------------
main_posterior_results_plot <-
  plot_grid(
    R0_plot + theme(legend.position = c(0.05, 0.825)),
    Rt_plot,
    ifr_t_plot,
    alpha_t_plot,
    weekly_latent_case_ratio_plot,
    cummulative_latent_case_ratio_plot,
    align = "hv",
    nrow = 3,
    ncol = 2
  )

save_plot(
  filename = "~/Documents/oc_covid19_stemr_manuscript/figures/main_posterior_results_plot.pdf",
  plot = main_posterior_results_plot,
  ncol = 2,
  nrow = 3,
  device = cairo_pdf
)

write_rds(R0_plot, "ISBA_2021/plots_rds/R0_plot.rds")
write_rds(Rt_plot, "ISBA_2021/plots_rds/Rt_plot.rds")
write_rds(ifr_t_plot, "ISBA_2021/plots_rds/ifr_t_plot.rds")
write_rds(alpha_t_plot, "ISBA_2021/plots_rds/alpha_t_plot.rds")
write_rds(weekly_latent_case_ratio_plot, "ISBA_2021/plots_rds/weekly_latent_case_ratio_plot.rds")
write_rds(cummulative_latent_case_ratio_plot, "ISBA_2021/plots_rds/cummulative_latent_case_ratio_plot.rds")
