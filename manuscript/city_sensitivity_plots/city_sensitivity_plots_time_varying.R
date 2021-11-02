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

model_name_conversion <- c(
  `2021-02-28` = "Orange County",
  double_initial_counts = "Sensitivity 1",
  higher_ifr = "Sensitivity 4",
  huntington_beach = "Huntington Beach",
  irvine = "Irvine",
  longer_infectious_period = "Sensitivity 2",
  lower_R0 = "Sensitivity 3",
  lower_alpha = "Sensitivity 5",
  santa_ana = "Santa Ana"
)

# initial_cases_county <-
#   read_csv("final_models/data/oc_data.csv") %>%
#   filter(date < "2020-03-30") %>%
#   pull(cases) %>%
#   sum()
#
# initial_cases_city <-
#   read_csv("~/Documents/uci_covid_modeling2/data/oc_city_data.csv") %>%
#   filter(date < "2020-03-30",
#          city %in% model_name_conversion) %>%
#   group_by(city) %>%
#   summarize(initial_cases = sum(cases)) %>%
#   deframe()

# Weekly Latent:Case Ratio Plot -------------------------------------------
# Cumulative Latent:Case Ratio Plot ---------------------------------------

stemr_all_time_varying_params <-
  dir_ls("final_models/results/stemr_all_params/") %>%
  map(~read_rds(.) %>%
        select(starts_with("."), date, beta_t, dur_infec, S, E, I, D, ifr_t, alpha0_t, E2I, cases, deaths)) %>%
  `names<-`(., names(.) %>% path_file() %>% path_ext_remove() %>% str_sub(end = -18)) %>%
  `names<-`(., model_name_conversion[names(.)]) %>%
  imap_dfr(~mutate(.x, model = .y)) %>%
  mutate(popsize = case_when(model == "Irvine" ~ 264975L,
                             model == "Santa Ana" ~ 198471L,
                             model == "Huntington Beach" ~ 358874L,
                             TRUE ~ 3175692L),
         initial_cases = case_when(model == "Irvine" ~ 55L,
                                   model == "Santa Ana" ~ 56L,
                                   model == "Huntington Beach" ~ 35L,
                                   TRUE ~ 760L)) %>%
  mutate(log_popsize = log(popsize)) %>%
  replace_na(list(cases = 0,
                  deaths = 0)) %>%
  group_by(.chain, .iteration, .draw, model) %>%
  mutate(cumulative_cases = cumsum(cases) + initial_cases,
         cumulative_deaths = cumsum(deaths)) %>%
  ungroup() %>%
  transmute(
    date = date,
    model = model,
    R0 = exp(log(beta_t) + log_popsize + log(dur_infec)),
    Rt = R0 * S / popsize,
    ifr_t = ifr_t,
    alpha0_t = alpha0_t,
    latent_case_ratio = E2I / cases,
    cumulative_incidence = popsize - S,
    cumalitive_latent_case_ratio = cumulative_incidence / cumulative_cases,
    prevalence = E + I,
    cumulative_latent_deaths = D,
    cumulative_cases,
    cumulative_deaths
    ) %>%
  mutate(model = fct_relevel(factor(model), "Orange County"))

stemr_all_time_varying_params_city <- stemr_all_time_varying_params %>%
  filter(str_starts(model, "Sensitivity", negate = T))

stemr_all_time_varying_params_sensitivity <- stemr_all_time_varying_params %>%
  mutate(model = fct_recode(model, "Original" = "Orange County")) %>%
  filter(str_starts(model, "Sensitivity", negate = F) | model == "Original")

rm(stemr_all_time_varying_params)


# R0 Plots ----------------------------------------------------------------
R0_city_plot <-
  stemr_all_time_varying_params_city %>%
  select(model, R0, date) %>%
  group_by(date, model) %>%
  median_qi(.width = ci_width) %>%
  drop_na() %>%
  ggplot(aes(date, R0, ymin = .lower, ymax = .upper)) +
  facet_wrap(. ~ model) +
  # geom_interval(
  #   show_point = T,
  #   shape = "_",
  #   point_color = "black",
  #   fatten_point = 1
  # ) +
  geom_lineribbon(
    color = "steelblue4") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_brewer() +
  scale_color_brewer() +
  scale_y_continuous(
    name = expression(R[0]),
    breaks = seq(0, 3, by = 0.25),
    limits = c(.5, 2.15)
  ) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks3,
               date_labels = date_labels) +
  theme(legend.position = "none") +
  ggtitle("Posterior Basic Reproductive Number")

R0_sensitivity_plot <-
  stemr_all_time_varying_params_sensitivity %>%
  select(model, R0, date) %>%
  group_by(date, model) %>%
  median_qi(.width = ci_width) %>%
  drop_na() %>%
  ggplot(aes(date, R0, ymin = .lower, ymax = .upper)) +
  facet_wrap(. ~ model) +
  # geom_interval(
  #   show_point = T,
  #   shape = "_",
  #   point_color = "black",
  #   fatten_point = 1
  # ) +
  geom_lineribbon(
    color = "steelblue4") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_brewer() +
  scale_color_brewer() +
  scale_y_continuous(
    name = expression(R[0]),
    breaks = seq(0, 3, by = 0.25),
    limits = c(.5, 2.25)
  ) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks3,
               date_labels = date_labels) +
  theme(legend.position = "none") +
  ggtitle("Posterior Basic Reproductive Number")

# Rt Plots ----------------------------------------------------------------
Rt_city_plot <-
  stemr_all_time_varying_params_city %>%
  select(date, model, Rt) %>%
  group_by(date, model) %>%
  median_qi(.width = ci_width) %>%
  drop_na() %>%
  ggplot(aes(date, Rt, ymin = .lower, ymax = .upper)) +
  facet_wrap(. ~ model) +
  # geom_lineribbon(color = "steelblue4") +
  # geom_interval(
  #   show_point = T,
  #   shape = "_",
  #   point_color = "black",
  #   fatten_point = 1
  # ) +
  geom_lineribbon(
    color = "steelblue4") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  # scale_fill_brewer() +
  scale_color_brewer() +
  scale_y_continuous(
    name = expression(R[t]),
    breaks = seq(0, 3, by = 0.25),
    limits = c(.5, 2.05)
  ) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks3,
               date_labels = date_labels) +
  theme(legend.position = "none") +
  # ggtitle(expression(bold(paste("Posterior ", R[t]))))
  ggtitle("Posterior Effective Reproductive Number")


Rt_sensitivity_plot <-
  stemr_all_time_varying_params_sensitivity %>%
  select(date, model, Rt) %>%
  group_by(date, model) %>%
  median_qi(.width = ci_width) %>%
  drop_na() %>%
  ggplot(aes(date, Rt, ymin = .lower, ymax = .upper)) +
  facet_wrap(. ~ model) +
  # geom_lineribbon(color = "steelblue4") +
  # geom_interval(
  #   show_point = T,
  #   shape = "_",
  #   point_color = "black",
  #   fatten_point = 1
  # ) +
  geom_lineribbon(
    color = "steelblue4") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_brewer() +
  # scale_color_brewer() +
  scale_y_continuous(
    name = expression(R[t]),
    breaks = seq(0, 3, by = 0.25),
    limits = c(.5, 2.05)
  ) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks3,
               date_labels = date_labels) +
  theme(legend.position = "none") +
  # ggtitle(expression(bold(paste("Posterior ", R[t]))))
  ggtitle("Posterior Effective Reproductive Number")



# IFR Plots ---------------------------------------------------------------
ifr_t_city_plot <-
  stemr_all_time_varying_params_city %>%
  select(date, model, ifr_t) %>%
  group_by(date, model) %>%
  median_qi(.width = ci_width) %>%
  drop_na() %>%
  ggplot(aes(date, ifr_t, ymin = .lower, ymax = .upper)) +
  facet_wrap(. ~ model) +
  # geom_lineribbon(color = "steelblue4") +
  # geom_interval(
  #   show_point = T,
  #   shape = "_",
  #   point_color = "black",
  #   fatten_point = 1
  # ) +
  geom_lineribbon(
    color = "steelblue4") +
  scale_fill_brewer() +
  # scale_color_brewer() +
  theme(legend.position = "none") +
  scale_y_continuous(
    name = "IFR",
    labels = function(.)
      scales::percent(., accuracy = 0.1),
    # trans = scales::logit_trans(),
    breaks = seq(0, 0.05, by = 0.005),
    limits = c(0, NA)
  ) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks3,
               date_labels = date_labels) +
  ggtitle("Posterior Infection Fatality Ratio")

ifr_t_sensitivity_plot <-
  stemr_all_time_varying_params_sensitivity %>%
  select(date, model, ifr_t) %>%
  group_by(date, model) %>%
  median_qi(.width = ci_width) %>%
  drop_na() %>%
  ggplot(aes(date, ifr_t, ymin = .lower, ymax = .upper)) +
  facet_wrap(. ~ model) +
  # geom_lineribbon(color = "steelblue4") +
  # geom_interval(
  #   show_point = T,
  #   shape = "_",
  #   point_color = "black",
  #   fatten_point = 1
  # ) +
  geom_lineribbon(
    color = "steelblue4") +
  scale_fill_brewer() +
  # scale_color_brewer() +
  theme(legend.position = "none") +
  scale_y_continuous(
    name = "IFR",
    labels = function(.)
      scales::percent(., accuracy = 0.1),
    # trans = scales::logit_trans(),
    breaks = seq(0, 0.05, by = 0.001),
    limits = c(0, NA)
  ) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks3,
               date_labels = date_labels) +
  ggtitle("Posterior Infection Fatality Ratio")




# alpha_t_plots -----------------------------------------------------------
alpha_t_city_plot <-
  stemr_all_time_varying_params_city %>%
  select(date, model, alpha0_t) %>%
  group_by(date, model) %>%
  median_qi(.width = ci_width) %>%
  drop_na() %>%
  ggplot(aes(date, alpha0_t, ymin = .lower, ymax = .upper)) +
  facet_wrap(. ~ model) +
  # geom_lineribbon(color = "steelblue4") +
  # geom_interval(
  #   show_point = T,
  #   shape = "_",
  #   point_color = "black",
  #   fatten_point = 1
  # ) +
  geom_lineribbon(
    color = "steelblue4") +
  scale_fill_brewer() +
  # scale_color_brewer() +
  theme(legend.position = "none") +
  scale_y_continuous(name = expression(alpha),
                     # scale_y_continuous(name = expression(alpha[0]),
                     # trans = scales::log_trans(),
                     breaks = 0:5) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks3,
               date_labels = date_labels) +
  ggtitle("Posterior \u03B1")


alpha_t_sensitivity_plot <-
  stemr_all_time_varying_params_sensitivity %>%
  select(date, model, alpha0_t) %>%
  group_by(date, model) %>%
  median_qi(.width = ci_width) %>%
  drop_na() %>%
  ggplot(aes(date, alpha0_t, ymin = .lower, ymax = .upper)) +
  facet_wrap(. ~ model) +
  # geom_lineribbon(color = "steelblue4") +
  # geom_interval(
  #   show_point = T,
  #   shape = "_",
  #   point_color = "black",
  #   fatten_point = 1
  # ) +
  geom_lineribbon(
    color = "steelblue4") +
  scale_fill_brewer() +
  # scale_color_brewer() +
  theme(legend.position = "none") +
  scale_y_continuous(name = expression(alpha),
                     # scale_y_continuous(name = expression(alpha[0]),
                     # trans = scales::log_trans(),
                     breaks = 0:5) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks3,
               date_labels = date_labels) +
  ggtitle("Posterior \u03B1")


# Weekly Latent:Case Ratio Plots ------------------------------------------
weekly_latent_case_ratio_city_plot <-
  stemr_all_time_varying_params_city %>%
  select(date, model, latent_case_ratio) %>%
  group_by(date, model) %>%
  median_qi(.width = ci_width) %>%
  drop_na() %>%
  ggplot(aes(date, latent_case_ratio, ymin = .lower, ymax = .upper)) +
  facet_wrap(. ~ model) +
  # geom_lineribbon(show.legend = F) +
  # geom_interval(
  #   show_point = T,
  #   shape = "_",
  #   point_color = "black",
  #   fatten_point = 1
  # ) +
  geom_lineribbon(
    color = "steelblue4") +
  scale_fill_brewer() +
  # scale_color_brewer() +
  theme(legend.position = "none") +
  scale_y_continuous(
    name = "Weekly Latent:Case Ratio",
    labels = function(.)
      str_c(scales::comma(., accuracy = 1), ":1"),
    breaks = seq(0, 240, 5),
    limits = c(0, 240)
  ) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks3,
               date_labels = date_labels) +
  ggtitle("Posterior Weekly Latent:Case Ratio")


weekly_latent_case_ratio_sensitivity_plot <-
  stemr_all_time_varying_params_sensitivity %>%
  select(date, model, latent_case_ratio) %>%
  group_by(date, model) %>%
  median_qi(.width = ci_width) %>%
  drop_na() %>%
  ggplot(aes(date, latent_case_ratio, ymin = .lower, ymax = .upper)) +
  facet_wrap(. ~ model) +
  # geom_lineribbon(show.legend = F) +
  # geom_interval(
  #   show_point = T,
  #   shape = "_",
  #   point_color = "black",
  #   fatten_point = 1
  # ) +
  geom_lineribbon(
    color = "steelblue4") +
  scale_fill_brewer() +
  # scale_color_brewer() +
  theme(legend.position = "none") +
  scale_y_continuous(
    name = "Weekly Latent:Case Ratio",
    labels = function(.)
      str_c(scales::comma(., accuracy = 1), ":1"),
    breaks = seq(0, 240, 5),
    limits = c(1, 40)
  ) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks3,
               date_labels = date_labels) +
  ggtitle("Posterior Weekly Latent:Case Ratio")

# Cumulative Latent:Case Ratio Plots --------------------------------------
cummulative_latent_case_ratio_city_plot <-
  stemr_all_time_varying_params_city %>%
  select(date, model, cumalitive_latent_case_ratio) %>%
  group_by(date, model) %>%
  median_qi(.width = ci_width) %>%
  drop_na() %>%
  ggplot(aes(
    date,
    cumalitive_latent_case_ratio,
    ymin = .lower,
    ymax = .upper
  )) +
  facet_wrap(. ~ model) +
  # geom_lineribbon(show.legend = F) +
  # geom_interval(
  #   show_point = T,
  #   shape = "_",
  #   point_color = "black",
  #   fatten_point = 1
  # ) +
  geom_lineribbon(
    color = "steelblue4") +
  scale_fill_brewer() +
  # scale_color_brewer() +
  theme(legend.position = "none") +
  scale_y_continuous(
    name = "Cumulative Latent:Case Ratio",
    labels = function(.)
      str_c(scales::comma(., accuracy = 1), ":1"),
    breaks = seq(0, 100, 5),
    limits = c(1, 80)
  ) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks3,
               date_labels = date_labels) +
  ggtitle("Posterior Cumulative Latent:Case Ratio")


cummulative_latent_case_ratio_sensitivity_plot <-
  stemr_all_time_varying_params_sensitivity %>%
  select(date, model, cumalitive_latent_case_ratio) %>%
  group_by(date, model) %>%
  median_qi(.width = ci_width) %>%
  drop_na() %>%
  ggplot(aes(
    date,
    cumalitive_latent_case_ratio,
    ymin = .lower,
    ymax = .upper
  )) +
  facet_wrap(. ~ model) +
  # geom_lineribbon(show.legend = F) +
  # geom_interval(
  #   show_point = T,
  #   shape = "_",
  #   point_color = "black",
  #   fatten_point = 1
  # ) +
  geom_lineribbon(
    color = "steelblue4") +
  scale_fill_brewer() +
  # scale_color_brewer() +
  theme(legend.position = "none") +
  scale_y_continuous(
    name = "Cumulative Latent:Case Ratio",
    labels = function(.)
      str_c(scales::comma(., accuracy = 1), ":1"),
    breaks = seq(0, 100, 5),
    limits = c(1, 30)
  ) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks3,
               date_labels = date_labels) +
  ggtitle("Posterior Cumulative Latent:Case Ratio")


# Save plots --------------------------------------------------------------
city_plot_names <- ls()[str_ends(ls(), "city_plot")]
sensitivity_plot_names <- ls()[str_ends(ls(), "sensitivity_plot")]

sensitivity_plot_names %>%
  map(~save_plot(plot = get(.), ncol = 3, nrow = 2,
                 filename = path("~/Documents/oc_covid19_stemr_manuscript/figures", . , ext = "pdf"),
                 device = cairo_pdf,
                 # base_asp = (8.5 * 2) / (11 * 3)))
                 base_asp = 1))

# city_plot_names %>%
#   map(~save_plot(plot = get(.), ncol = 2, nrow = 2,
#                  filename = path("~/Documents/stemr_oc/manuscript/city_sensitivity_plots/", . , ext = "pdf"),
#                  device = cairo_pdf,
#                  # base_asp = 16/9)),
#                  base_asp = 1))


# sensitivity_plot_names %>%
#   map(~save_plot(plot = get(.), ncol = 3, nrow = 2,
#                  filename = path("~/Documents/stemr_oc/manuscript/city_sensitivity_plots/", . , ext = "pdf"),
#                  device = cairo_pdf,
#                  # base_asp = (16 * 2) / (9 * 3)))
#                  base_asp = 2 / 3))
