library(extraDistr)
library(tidybayes)
library(broom)

nature_dat <-  tibble(age = c(2, 7, 12, 17, 22, 27, 32, 37, 42, 47, 52, 57, 62, 67, 72, 77, 85),
                      ifr = c(3e-05, 1e-05, 1e-05, 3e-05, 6e-05, 0.00013, 0.00024, 4e-04,
                              0.00075, 0.00121, 0.00207, 0.00323, 0.00456, 0.01075, 0.01674,
                              0.03203, 0.08292))


nature_dat_model <- glm(formula = ifr ~ age, family = gaussian(link = "logit"), data = nature_dat %>% mutate(age = age - 42))

oc_data <- read_rds("data/oc_data.rds")

lump_oc_data <-
  function(oc_data,
           time_interval_in_days,
           first_day = "0000-01-01",
           last_day = "9999-12-31") {
    oc_data %>%
      filter(date >= lubridate::ymd(first_day),
             date <= lubridate::ymd(last_day)) %>%
      group_by(lump = as.integer(floor((max(date) - date) / time_interval_in_days))) %>%
      filter(n() == time_interval_in_days) %>%
      dplyr::summarize(start_date = min(date),
                end_date = max(date),
                cases = sum(cases),
                tests = sum(tests),
                deaths = sum(deaths)) %>%
      dplyr::select(-lump) %>%
      arrange(start_date)
  }

lumped_oc_data <- lump_oc_data(oc_data, time_interval_in_days = 7, first_day = "2020-03-14", last_day = "2020-10-16") %>%
  left_join(select(oc_data, end_date = date, lagged_pos_age, lagged_pos_ifr))

dat <- lumped_oc_data %>%
  mutate(time = as.numeric((end_date - min(start_date) + 1L) / 7)) %>%
  select(time, cases, tests, deaths, lagged_pos_age, lagged_pos_ifr)


nature_dat_model %>%
  augment() %>%
  mutate(ifr = expit(`logit(ifr)`),
         fitted_ifr = expit(.fitted)) %>%
  select(age, ifr, fitted_ifr) %>%
  mutate(age = age + 42) %>%
  pivot_longer(-age) %>%
  ggplot(aes(age, value, group = name, color = name)) +
  geom_line() +
  geom_point() +
  scale_y_log10()


n <- 2000
eta1 <- rnorm(n, -7.3, 0.5)
eta2 <- rtnorm(n, 0.1102, 0.05, a = 0)



tmp <- tibble(age = seq(0, 85, 5)) %>%
  mutate(ifr = map(age, function(age) expit(eta1 + eta2 * (age - 42)))) %>%
  select(age, ifr) %>%
  unnest(ifr) %>%
  group_by(age) %>%
  median_qi(.width = c(0.5, 0.8, 0.9))


ggplot() +
  geom_lineribbon(data = tmp, mapping = aes(age, ifr, ymin = .lower, ymax = .upper)) +
  geom_point(data = nature_dat, mapping = aes(age, ifr), color = 2) +
  scale_fill_brewer() +
  cowplot::theme_minimal_grid() +
  scale_y_log10(labels = scales::percent) +
  theme(legend.position = "none") +
  ggtitle("IFR Prior Derived from Average Age")


ifr_prior_our_data <- dat %>%
  select(time, lagged_pos_age) %>%
  mutate(ifr = map(lagged_pos_age, function(lagged_pos_age) expit(eta1 + eta2 * (lagged_pos_age - 42)))) %>%
  select(time, ifr) %>%
  unnest(ifr) %>%
  group_by(time) %>%
  median_qi(.width = c(0.5, 0.8, 0.9))


age_plot <- ggplot(ifr_prior_our_data, aes(time, ifr, ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  scale_fill_brewer() +
  cowplot::theme_minimal_grid() +
  scale_y_log10(labels = scales::percent) +
  theme(legend.position = "none") +
  ggtitle("IFR Prior Derived from expected IFR of avg Age")

eta3 <- rnorm(n, sd = 0.4)

ifr_prior_our_data2 <- dat %>%
  select(time, lagged_pos_ifr) %>%
  mutate(ifr = map(lagged_pos_ifr, function(lagged_pos_ifr) expit(eta3 + logit(lagged_pos_ifr)))) %>%
  select(time, ifr) %>%
  unnest(ifr) %>%
  group_by(time) %>%
  median_qi(.width = c(0.5, 0.8, 0.9))

ifr_plot <- ggplot(ifr_prior_our_data2, aes(time, ifr, ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  scale_fill_brewer() +
  cowplot::theme_minimal_grid() +
  scale_y_log10(labels = scales::percent) +
  theme(legend.position = "none") +
  ggtitle("IFR Prior Derived from expected IFR")


cowplot::plot_grid(age_plot, ifr_plot, ncol = 1)
