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
ci_width <- 0.95
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
R0_plot_chain <-
  stemr_all_params %>%
  mutate(R0 = exp(log(beta_t) + log_popsize + log(dur_infec)),
         prop_s = S / popsize) %>%
  select(.chain, R0, date = end_date) %>%
  group_by(date, .chain) %>%
  median_qi(.width = ci_width) %>%
  drop_na() %>%
  mutate(.width = scales::percent(.width),
         .chain = as_factor(.chain)) %>%
  ggplot(aes(date, R0, ymin = .lower, ymax = .upper, group = .chain, fill = .chain)) +
  geom_lineribbon(alpha = 0.1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_discrete("Chain") +
  scale_y_continuous(
    name = expression(R[0]),
    breaks = seq(0, 3, by = 0.25),
    limits = c(.5, 2.5)
  ) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks,
               date_labels = date_labels) +
  ggtitle("Posterior Basic Reproductive Number",
          subtitle = "95% Credible Intervals")


# Rt Plot -----------------------------------------------------------------
Rt_plot_chain <-
  stemr_all_params %>%
  mutate(
    R0 = exp(log(beta_t) + log_popsize + log(dur_infec)),
    prop_s = S / popsize,
    Rt = R0 * prop_s
  ) %>%
  select(.chain, Rt, date = end_date) %>%
  group_by(date, .chain) %>%
  median_qi(.width = ci_width) %>%
  drop_na() %>%
  mutate(.chain = as_factor(.chain)) %>%
  ggplot(aes(date, Rt, ymin = .lower, ymax = .upper, group = .chain, fill = .chain)) +
  geom_lineribbon(alpha = 0.1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_discrete("Chain") +
  scale_y_continuous(
    name = expression(R[t]),
    breaks = seq(0, 3, by = 0.25),
    limits = c(.5, 2.5)
  ) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks,
               date_labels = date_labels) +
  theme(legend.position = "none") +
  ggtitle("Posterior Effective Reproductive Number",
        subtitle = "95% Credible Intervals")



# IFR Plot ----------------------------------------------------------------
ifr_t_plot_chain <-
  stemr_all_params %>%
  select(.chain, date = end_date, ifr_t) %>%
  group_by(.chain, date) %>%
  median_qi(.width = ci_width) %>%
  drop_na() %>%
  mutate(.chain = as_factor(.chain)) %>%
  ggplot(aes(date, ifr_t, ymin = .lower, ymax = .upper, group = .chain, fill = .chain)) +
  geom_lineribbon(alpha = 0.1) +
  scale_fill_discrete() +
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
  ggtitle("Posterior Infection Fatality Ratio",
          subtitle = "95% Credible Intervals")


# Alpha_t plot ------------------------------------------------------------
alpha_t_plot_chain <-
  stemr_all_params %>%
  select(.chain, date = end_date, alpha0_t) %>%
  group_by(.chain, date) %>%
  median_qi(.width = ci_width) %>%
  drop_na() %>%
  mutate(.chain = as_factor(.chain)) %>%
  ggplot(aes(date, alpha0_t, ymin = .lower, ymax = .upper, group = .chain, fill = .chain)) +
  geom_lineribbon(alpha = 0.1) +
  scale_fill_discrete() +
  theme(legend.position = "none") +
  scale_y_continuous(name = expression(alpha),
                     # scale_y_continuous(name = expression(alpha[0]),
                     # trans = scales::log_trans(),
                     breaks = 0:5) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks,
               date_labels = date_labels) +
  ggtitle("Posterior \u03B1",
          subtitle = "95% Credible Intervals")



# Weekly Latent:Case Ratio Plot -------------------------------------------
weekly_latent_case_ratio_plot_chain <-
  weekly_case_ratio %>%
  select(.chain, date, latent_case_ratio) %>%
  group_by(.chain, date) %>%
  median_qi(.width = ci_width) %>%
  mutate(.chain = as_factor(.chain)) %>%
  ggplot(aes(date, latent_case_ratio, ymin = .lower, ymax = .upper, group = .chain, fill = .chain)) +
  geom_lineribbon(alpha = 0.1, show.legend = F) +
  scale_fill_discrete() +
  theme(legend.position = "none") +
  scale_y_continuous(
    name = "Weekly Latent:Case Ratio",
    labels = function(.)
      str_c(scales::comma(., accuracy = 1), ":1"),
    breaks = seq(0, 40, 5),
    limits = c(1, 33)
  ) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks,
               date_labels = date_labels) +
  ggtitle("Posterior Weekly Latent:Case Ratio",
          subtitle = "95% Credible Intervals")


# Cumulative Latent:Case Ratio Plot ---------------------------------------
cummulative_latent_case_ratio_plot_chain <-
  cumulative_case_ratio %>%
  select(.chain, date, cumalitive_latent_case_ratio) %>%
  group_by(.chain, date) %>%
  median_qi(.width = ci_width) %>%
  mutate(.chain = as_factor(.chain)) %>%
  ggplot(aes(
    date,
    cumalitive_latent_case_ratio,
    ymin = .lower,
    ymax = .upper, group = .chain, fill = .chain
  )) +
  geom_lineribbon(alpha = 0.1, show.legend = F) +
  scale_fill_discrete() +
  theme(legend.position = "none") +
  scale_y_continuous(
    name = "Cumulative Latent:Case Ratio",
    labels = function(.)
      str_c(scales::comma(., accuracy = 1), ":1"),
    breaks = seq(0, 100, 5),
    limits = c(1, 33)
  ) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks,
               date_labels = date_labels) +
  ggtitle("Posterior Cumulative Latent:Case Ratio",
          subtitle = "95% Credible Intervals")



# Save Plot ---------------------------------------------------------------
main_posterior_results_plot_chain <-
  plot_grid(
    R0_plot_chain + theme(legend.position = c(0.05, 0.825), legend.direction="horizontal"),
    Rt_plot_chain,
    ifr_t_plot_chain,
    alpha_t_plot_chain,
    weekly_latent_case_ratio_plot_chain,
    cummulative_latent_case_ratio_plot_chain,
    align = "hv",
    nrow = 3,
    ncol = 2
  )

save_plot(
  filename = "~/Documents/oc_covid19_stemr_manuscript/figures/main_posterior_results_plot_chain.pdf",
  plot = main_posterior_results_plot_chain,
  ncol = 2,
  nrow = 3,
  device = cairo_pdf
)
